"""
run_crossed_ensemble.py  --  4b: 600-member crossed ensemble for Helene.

Crosses the two DA-on perturbation arms:
    Forcing arm  : 30 met draws  (from run_perturbation_da_on.py Phase 2A)
    Hydro arm    : 20 state draws (from run_perturbation_da_on.py Phase 2B)

For each Helene issue time t0:
    Run 600 members = 30 forcing draws x 20 hydro-state draws.
    Each member (i, j) gets:
        - Initial state  = perturbed DA analysis state  (draw j of hydro arm)
        - Met forcing    = perturbed precip/PET sequence (draw i of forcing arm)

This gives the full uncertainty envelope combining both sources —
the operational probabilistic forecast picture.

Output (per catchment):
    <out-dir>/<cat-id>/<cat-id>_crossed_ensemble.parquet
    Columns: issue_time, lead_hour, member_0000..member_0599

Usage:
    python3 run_crossed_ensemble.py \\
        --cat-id cat-1016300 \\
        --forcing-dir  <NWM retro dir> \\
        --obs-dir      /mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance \\
        --cfe-dir      /mnt/disk2/suma_helen_poster/cfe_py \\
        --config-file  /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \\
        --out-dir      /mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble \\
        --test-forcing-dir1 <NWM op 2023-Feb2024> \\
        --test-forcing-dir2 <NWM op Feb2024-Sep2025>

Strategy:
    Re-runs DA analysis (same as Phase 1 of run_perturbation_da_on.py) to get
    the per-issue-time analysis state. Then generates N_FORCING sets of perturbed
    forcing sequences and N_HYDRO perturbed states. Runs all 600 combinations.

    To avoid re-running Phase 1, if the Phase 1 snapshot file exists from a
    prior run_perturbation_da_on.py run, it is loaded directly.

NOTE: All special characters kept ASCII to avoid encoding issues on GPU.
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

ALBEDO   = 0.20
ALPHA_PT = 1.26

N_FORCING = 30
N_HYDRO   = 20
N_MEMBERS = N_FORCING * N_HYDRO   # 600

TEST_START = "2024-08-24 00:00:00"   # 1 month spinup before Helene window
TEST_END   = "2024-10-31 23:00:00"

HELENE_START   = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END     = pd.Timestamp("2024-09-30 06:00:00")
FORECAST_HOURS = 18

PRECIP_SIGMA = 0.15
PET_SIGMA    = 0.10
STATE_FRAC   = 0.05

CAT_ID          = None
OBS_FILE        = None
CFE_CONFIG_FILE = None
TEST_FORCING_FILE = None
OUT_DIR         = None
HARDCODED_R     = None
DIRECT_VARIANCE = False   # when True: R = krig_var directly (no Vrugt formula)
bmi_cfe         = None


def priestley_taylor_pet(srad_wm2, T_kelvin, alpha=ALPHA_PT):
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    pet = alpha * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def load_test_forcing():
    df = pd.read_csv(TEST_FORCING_FILE)
    df["date"] = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    df["total_precipitation"] = df["APCP_surface"] * 3600.0
    df["potential_evaporation"] = priestley_taylor_pet(
        df["DSWRF_surface"].values, df["TMP_2maboveground"].values)
    return df


def load_obs():
    df = pd.read_csv(OBS_FILE)
    date_col = next((c for c in df.columns
                     if c.lower() in ("date", "time", "datetime", "timestamp")), None)
    q_col = next((c for c in df.columns
                  if "qkrig" in c.lower() and "var" not in c.lower()), None)
    var_col = next((c for c in df.columns if "var" in c.lower()), None)
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(date_col).reset_index(drop=True)
    date_strs = df[date_col].dt.strftime("%Y-%m-%d %H:%M:%S")
    obs_dict = dict(zip(date_strs, df[q_col]))
    var_vals = df[var_col] if var_col else pd.Series([0.07] * len(df))
    var_dict = dict(zip(date_strs, var_vals))
    return obs_dict, var_dict


def write_model_config(best_params):
    """Write the BMI config file once; reused for every model instantiation."""
    with open(CFE_CONFIG_FILE) as f:
        cfg = json.load(f)
    cfg["forcing_file"]          = TEST_FORCING_FILE
    cfg["soil_params"]["bb"]     = best_params["bb"]
    cfg["soil_params"]["smcmax"] = best_params["smcmax"]
    cfg["soil_params"]["satdk"]  = best_params["satdk"]
    cfg["slop"]                  = best_params["slop"]
    cfg["max_gw_storage"]        = best_params["max_gw_storage"]
    cfg["expon"]                 = best_params["expon"]
    cfg["Cgw"]                   = best_params["Cgw"]
    cfg["K_lf"]                  = best_params["K_lf"]
    cfg["K_nash"]                = best_params["K_nash"]
    cfg["partition_scheme"]      = ("Schaake" if best_params["scheme"] <= 0.5
                                    else "Xinanjiang")
    tmp_cfg = str(OUT_DIR / f"{CAT_ID}_bmi_config_crossed_stable.json")
    with open(tmp_cfg, "w") as f:
        json.dump(cfg, f)
    return tmp_cfg


def make_model(cfg_path):
    """Instantiate and initialize a CFE model from a pre-written config path."""
    m = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
    m.initialize()
    return m


def get_state(m):
    return {
        "soil_m": m.soil_reservoir["storage_m"],
        "gw_m":   m.gw_reservoir["storage_m"],
        "nash0":  float(m.nash_storage[0]),
        "nash1":  float(m.nash_storage[1]),
    }


def set_state(m, state):
    m.soil_reservoir["storage_m"] = state["soil_m"]
    m.gw_reservoir["storage_m"]   = state["gw_m"]
    m.nash_storage[0]              = state["nash0"]
    m.nash_storage[1]              = state["nash1"]


def perturb_state(state, rng, frac=STATE_FRAC):
    return {
        "soil_m": max(state["soil_m"] * (1.0 + frac * rng.standard_normal()), 1e-6),
        "gw_m":   max(state["gw_m"]   * (1.0 + frac * rng.standard_normal()), 1e-6),
        "nash0":  max(state["nash0"]  * (1.0 + frac * rng.standard_normal()), 0.0),
        "nash1":  max(state["nash1"]  * (1.0 + frac * rng.standard_normal()), 0.0),
    }


def run_model_step(m, precip_mmh, pet_mmh):
    m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                precip_mmh / 1000.0)
    m.set_value("water_potential_evaporation_flux",
                pet_mmh / 1000.0 / 3600.0)
    m.update()
    return m.get_value("land_surface_water__runoff_depth") * 1000.0


def enkf_update(state, q_sim, y_obs, krig_var, rng):
    if HARDCODED_R is not None:
        R = HARDCODED_R
    elif DIRECT_VARIANCE:
        R = max(krig_var, 1e-6)
    else:
        R = max((0.10 * max(y_obs, 0.0)) ** 2 + 0.001 * krig_var, 1e-6)
    R = max(R, 1e-6)
    P_yy = max(q_sim * 0.01, 1e-6)
    K = P_yy / (P_yy + R)
    scale = K * (y_obs - q_sim) / max(abs(q_sim), 1e-6)
    return {
        "soil_m": max(state["soil_m"] * (1.0 + scale), 1e-6),
        "gw_m":   max(state["gw_m"]   * (1.0 + scale), 1e-6),
        "nash0":  max(state["nash0"]  + state["nash0"] * scale, 0.0),
        "nash1":  max(state["nash1"]  + state["nash1"] * scale, 0.0),
    }


def run_phase1_da(cfg_path, df_test, obs_dict, var_dict, rng):
    """Run DA through test period; return snapshot dict {timestamp -> state}."""
    print(f"  Phase 1: DA analysis through test period...")
    model = make_model(cfg_path)
    snapshots = {}
    for _, row in df_test.iterrows():
        date_str = row["date"]
        P = float(row["total_precipitation"])
        E = float(row["potential_evaporation"])
        q_sim = run_model_step(model, P, E)
        y_obs = obs_dict.get(date_str, np.nan)
        if not np.isnan(y_obs):
            krig_var = var_dict.get(date_str, 0.07)
            state = get_state(model)
            updated = enkf_update(state, q_sim, y_obs, krig_var, rng)
            set_state(model, updated)
        date_ts = pd.Timestamp(date_str)
        if HELENE_START <= date_ts <= HELENE_END:
            snapshots[date_ts] = get_state(model)
    model.finalize()
    print(f"  Phase 1 done: {len(snapshots)} issue-time snapshots")
    return snapshots


def main():
    global CAT_ID, OBS_FILE, CFE_CONFIG_FILE, TEST_FORCING_FILE
    global OUT_DIR, HARDCODED_R, DIRECT_VARIANCE, bmi_cfe

    parser = argparse.ArgumentParser()
    parser.add_argument("--cat-id",            required=True)
    parser.add_argument("--forcing-dir",       required=True)
    parser.add_argument("--obs-dir",           required=True)
    parser.add_argument("--cfe-dir",           required=True)
    parser.add_argument("--config-file",       required=True)
    parser.add_argument("--out-dir",           required=True)
    parser.add_argument("--test-forcing-dir1", default=None)
    parser.add_argument("--test-forcing-dir2", default=None)
    parser.add_argument("--hardcoded-r",       type=float, default=None,
                        help="Fixed R value for EnKF (e.g. 0.07). "
                             "Omit to use dynamic Vrugt R: R=(0.10*y)^2 + 0.001*krig_var")
    parser.add_argument("--direct-variance",   action="store_true", default=False,
                        help="Use krig_var directly as R (no Vrugt scaling)")
    parser.add_argument("--rng-seed",          type=int, default=None,
                        help="RNG seed for reproducibility (default: hash of cat-id)")
    parser.add_argument("--n-forcing",         type=int, default=N_FORCING,
                        help="Number of met forcing draws (default 30)")
    parser.add_argument("--n-hydro",           type=int, default=N_HYDRO,
                        help="Number of hydro-state draws (default 20)")
    args = parser.parse_args()

    CAT_ID        = args.cat_id
    OBS_FILE      = os.path.join(args.obs_dir, f"{CAT_ID}.csv")
    CFE_CONFIG_FILE = args.config_file
    HARDCODED_R     = args.hardcoded_r
    DIRECT_VARIANCE = args.direct_variance
    OUT_DIR         = Path(args.out_dir) / CAT_ID
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    n_fa = args.n_forcing
    n_ha = args.n_hydro
    n_total = n_fa * n_ha
    print(f"[{CAT_ID}] Crossed ensemble: {n_fa} forcing x {n_ha} hydro = {n_total} members")

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    # Reuse combined forcing CSV from run_perturbation_da_on.py if it exists
    # (avoids writing a duplicate large file to disk)
    da_on_combined = (Path(args.out_dir).parent / "v2_perturbation_da_on" /
                      CAT_ID / f"{CAT_ID}_nwm_operational_combined.csv")
    crossed_combined = OUT_DIR / f"{CAT_ID}_nwm_operational_combined.csv"

    if da_on_combined.exists():
        TEST_FORCING_FILE = str(da_on_combined)
        print(f"  Reusing combined forcing: {TEST_FORCING_FILE}")
    elif crossed_combined.exists():
        TEST_FORCING_FILE = str(crossed_combined)
        print(f"  Reusing combined forcing: {TEST_FORCING_FILE}")
    elif args.test_forcing_dir1 and args.test_forcing_dir2:
        f1 = os.path.join(args.test_forcing_dir1, f"{CAT_ID}.csv")
        f2 = os.path.join(args.test_forcing_dir2, f"{CAT_ID}.csv")
        df1, df2 = pd.read_csv(f1), pd.read_csv(f2)
        combined = (pd.concat([df1, df2], ignore_index=True)
                    .drop_duplicates(subset="time").sort_values("time"))
        TEST_FORCING_FILE = str(crossed_combined)
        combined.to_csv(TEST_FORCING_FILE, index=False)
        print(f"  Wrote combined forcing: {TEST_FORCING_FILE}")
    else:
        print("Need --test-forcing-dir1 and --test-forcing-dir2"); return

    best_params_file = OUT_DIR / f"{CAT_ID}_best_params.json"
    if not best_params_file.exists():
        # Also check the DA results dir from run_perturbation_da_on.py
        alt = Path(args.out_dir).parent / "v2_perturbation_da_on" / CAT_ID / f"{CAT_ID}_best_params.json"
        if alt.exists():
            best_params_file = alt
        else:
            print(f"Missing best params at {best_params_file}"); return
    with open(best_params_file) as f:
        best_params = json.load(f)["best_parameters"]

    obs_dict, var_dict = load_obs()
    df_forcing = load_test_forcing()
    t_mask = (df_forcing["date"] >= TEST_START) & (df_forcing["date"] <= TEST_END)
    df_test = df_forcing[t_mask].reset_index(drop=True)

    # Write BMI config once — reused for all 90,600 model instantiations
    cfg_path = write_model_config(best_params)
    print(f"  BMI config written once: {cfg_path}")

    rng = np.random.default_rng(args.rng_seed if args.rng_seed is not None
                                else hash(CAT_ID) & 0x7fffffff)

    # Phase 1: get DA analysis snapshots
    # Check if a prior run already saved them (skip re-running DA if so)
    snap_cache = OUT_DIR / f"{CAT_ID}_da_snapshots.parquet"
    if snap_cache.exists():
        print(f"  Loading cached DA snapshots from {snap_cache}")
        df_snaps = pd.read_parquet(snap_cache)
        snapshots = {
            pd.Timestamp(row["issue_time"]): {
                "soil_m": row["soil_m"], "gw_m": row["gw_m"],
                "nash0": row["nash0"],   "nash1": row["nash1"],
            }
            for _, row in df_snaps.iterrows()
        }
    else:
        snapshots = run_phase1_da(cfg_path, df_test, obs_dict, var_dict, rng)

    issue_times = sorted(snapshots.keys())
    t0_to_idx = {t: i for i, t in enumerate(df_test["date"].apply(pd.Timestamp))}

    # ------------------------------------------------------------------ #
    # Crossed ensemble: for each issue_time, build n_fa forcing draws and  #
    # n_ha state draws, then run all n_fa x n_ha combinations.             #
    # ------------------------------------------------------------------ #
    print(f"[{CAT_ID}] Running crossed ensemble ({n_total} members x "
          f"{len(issue_times)} issue times)...")

    all_records = []
    all_hydro_draw_records  = []   # (issue_time, hydro_draw_j, states)
    all_forcing_draw_records = []  # (issue_time, forcing_draw_i, lead, P_pert, E_pert, scales)

    for t_idx, t0 in enumerate(issue_times):
        snap = snapshots[t0]
        t0_str = t0.strftime("%Y-%m-%d %H:%M:%S")
        if t0 not in t0_to_idx:
            continue
        start_idx = t0_to_idx[t0]

        # Pre-draw n_fa forcing perturbation sequences
        forcing_seqs = []
        for i in range(n_fa):
            seq = []
            for lead in range(1, FORECAST_HOURS + 1):
                fi = start_idx + lead
                if fi >= len(df_test):
                    seq.append((0.0, 0.0))
                    all_forcing_draw_records.append({
                        "issue_time": t0_str, "forcing_draw_i": i,
                        "lead_hour": lead,
                        "P_pert_mmh": 0.0, "E_pert_mmh": 0.0,
                        "precip_scale": np.nan, "pet_scale": np.nan,
                    })
                    continue
                row = df_test.iloc[fi]
                P = float(row["total_precipitation"])
                E = float(row["potential_evaporation"])
                mu_p = -0.5 * PRECIP_SIGMA ** 2
                P_pert = float(P * rng.lognormal(mu_p, PRECIP_SIGMA)) if P > 0 else 0.0
                E_pert = max(float(E * (1.0 + PET_SIGMA * rng.standard_normal())), 0.0)
                seq.append((P_pert, E_pert))
                all_forcing_draw_records.append({
                    "issue_time": t0_str, "forcing_draw_i": i,
                    "lead_hour": lead,
                    "P_pert_mmh": P_pert, "E_pert_mmh": E_pert,
                    "precip_scale": P_pert / P if P > 0 else np.nan,
                    "pet_scale":    E_pert / E if E > 0 else np.nan,
                })
            forcing_seqs.append(seq)

        # Pre-draw n_ha hydro-state perturbations
        state_draws = [snap if j == 0 else perturb_state(snap, rng)
                       for j in range(n_ha)]

        # Record each hydro draw's actual initial state
        for j, st in enumerate(state_draws):
            all_hydro_draw_records.append({
                "issue_time": t0_str, "hydro_draw_j": j,
                "soil_m": st["soil_m"], "gw_m": st["gw_m"],
                "nash0":  st["nash0"],  "nash1": st["nash1"],
            })

        # q_matrix shape: (FORECAST_HOURS, n_fa, n_ha)
        q_matrix = np.full((FORECAST_HOURS, n_fa, n_ha), np.nan)

        for i, f_seq in enumerate(forcing_seqs):
            for j, state_j in enumerate(state_draws):
                m = make_model(cfg_path)
                set_state(m, state_j)
                for lead in range(1, FORECAST_HOURS + 1):
                    P_pert, E_pert = f_seq[lead - 1]
                    q = run_model_step(m, P_pert, E_pert)
                    q_matrix[lead - 1, i, j] = q
                m.finalize()

        # Flatten 600-member matrix into rows
        for lead in range(1, FORECAST_HOURS + 1):
            row_data = {"issue_time": t0_str, "lead_hour": lead}
            mem_idx = 0
            for i in range(n_fa):
                for j in range(n_ha):
                    row_data[f"member_{mem_idx:04d}"] = q_matrix[lead - 1, i, j]
                    mem_idx += 1
            all_records.append(row_data)

        if (t_idx + 1) % 10 == 0 or (t_idx + 1) == len(issue_times):
            print(f"  {t_idx + 1}/{len(issue_times)} issue times complete")

    df_out = pd.DataFrame(all_records)
    out_path = OUT_DIR / f"{CAT_ID}_crossed_ensemble.parquet"
    df_out.to_parquet(out_path, index=False)
    print(f"Saved: {out_path}")
    print(f"  Shape: {df_out.shape}  "
          f"({df_out['issue_time'].nunique()} issue times x "
          f"{FORECAST_HOURS} leads x {n_total} members)")

    # ------------------------------------------------------------------ #
    # Provenance files — enable full ensemble traceability and restart     #
    # ------------------------------------------------------------------ #

    # 1. DA analysis snapshots — unperturbed state at each issue_time.
    #    Load any row and set_state() to restart from that init time.
    snap_records = [
        {"issue_time": t.strftime("%Y-%m-%d %H:%M:%S"), **s}
        for t, s in snapshots.items()
    ]
    snap_path = OUT_DIR / f"{CAT_ID}_da_snapshots.parquet"
    pd.DataFrame(snap_records).to_parquet(snap_path, index=False)
    print(f"Saved: {snap_path}  ({len(snap_records)} snapshots)")

    # 2. Member manifest — decoder ring: member_col -> (forcing_draw_i, hydro_draw_j).
    #    member_k uses forcing draw k//n_ha and hydro draw k%n_ha.
    manifest = [
        {"member": f"member_{k:04d}", "member_idx": k,
         "forcing_draw_i": k // n_ha, "hydro_draw_j": k % n_ha}
        for k in range(n_total)
    ]
    manifest_path = OUT_DIR / f"{CAT_ID}_member_manifest.csv"
    pd.DataFrame(manifest).to_csv(manifest_path, index=False)
    print(f"Saved: {manifest_path}  ({n_total} members)")

    # 3. Hydro draw states — actual perturbed initial state for each of
    #    n_ha draws at every issue_time. Links member -> initial conditions.
    hydro_path = OUT_DIR / f"{CAT_ID}_hydro_draw_states.parquet"
    pd.DataFrame(all_hydro_draw_records).to_parquet(hydro_path, index=False)
    print(f"Saved: {hydro_path}  ({len(all_hydro_draw_records)} rows)")

    # 4. Forcing draw sequences — actual P/PET perturbation values for each
    #    of n_fa draws at every (issue_time, lead_hour). Links member -> forcing.
    forcing_path = OUT_DIR / f"{CAT_ID}_forcing_draw_sequences.parquet"
    pd.DataFrame(all_forcing_draw_records).to_parquet(forcing_path, index=False)
    print(f"Saved: {forcing_path}  ({len(all_forcing_draw_records)} rows)")


if __name__ == "__main__":
    main()
