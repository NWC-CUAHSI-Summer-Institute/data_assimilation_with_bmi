"""
calibrate_catchment_cfe_nwm.py  --  F5, 20% gauge holdout

Single script combining:
  1. DA run (EnKF, re-kriged sigma^2 R formula) through the test period
     to obtain warm posterior states at each Helene issue time.
  2. 600-member crossed ensemble forecasts (18 h lead) at each issue time:
         member(i,j) = forcing draw i  x  hydro-state draw j
         i in [0, N_FORCING)   -- met perturbation (lognormal precip, Gaussian PET)
         j in [0, N_HYDRO)     -- initial state perturbation around DA posterior
     Total: 30 x 20 = 600 members per issue time.

No routing -- outputs are per-catchment streamflow in mm/h.

Outputs per catchment (in --out-dir/<cat-id>/):
    <cat>_crossed_ensemble.parquet     -- 600-member, 18-lead forecasts (mm/h)
    <cat>_da_snapshots.parquet         -- DA posterior state at each issue time
    <cat>_member_manifest.csv          -- decoder: member_col -> (forcing_i, hydro_j)
    <cat>_hydro_draw_states.parquet    -- actual perturbed initial states
    <cat>_forcing_draw_sequences.parquet -- actual P/PET perturbation values

R formula (F5): R(t) = krig_var(t)  (re-kriged per-hour variogram variance, direct)

Usage (single catchment):
    python3 calibrate_catchment_cfe_nwm.py --cat-id cat-1016300

Usage (batch over all 21 catchments):
    for cat in cat-1016279 ... cat-1016315; do
        python3 calibrate_catchment_cfe_nwm.py --cat-id $cat &
    done

All paths below are server defaults for F5 20% holdout.
Override any path via CLI flags.
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

# ── Ensemble design ───────────────────────────────────────────────────────────
N_FORCING = 30
N_HYDRO   = 20
N_MEMBERS = N_FORCING * N_HYDRO   # 600

FORECAST_HOURS = 18

# Perturbation magnitudes
PRECIP_SIGMA = 0.15   # lognormal sigma for precipitation
PET_SIGMA    = 0.10   # Gaussian sigma (fraction) for PET
STATE_FRAC   = 0.05   # fractional std for initial state draws

# ── Time windows ──────────────────────────────────────────────────────────────
SPINUP_START = "2023-02-01 00:00:00"   # match main DA spinup
SPINUP_END   = "2023-09-30 23:00:00"
TEST_START   = "2023-10-01 00:00:00"
TEST_END     = "2024-10-31 23:00:00"

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-30 06:00:00")

# ── CFE physics ───────────────────────────────────────────────────────────────
ALBEDO   = 0.20
ALPHA_PT = 1.26

# ── Server defaults (F5, 20% holdout) ────────────────────────────────────────
DEFAULT_OBS_DIR    = "/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig"
DEFAULT_PARAM_SRC  = "/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded"
DEFAULT_CFE_DIR    = "/mnt/disk2/suma_helen_poster/cfe_py"
DEFAULT_CONFIG     = "/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
DEFAULT_FORCING1   = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                      "output_03463300_nwmoperational/03463300/2023_2024_feb/forcings")
DEFAULT_FORCING2   = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                      "output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings")
DEFAULT_OUT_DIR    = "/mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble_f5_20pct"


# ── CFE helpers ───────────────────────────────────────────────────────────────

def priestley_taylor_pet(srad_wm2, T_kelvin, alpha=ALPHA_PT):
    T     = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam   = 2.501 - 0.002361 * T
    pet   = alpha * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def load_forcing(path):
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    df["total_precipitation"] = df["APCP_surface"] * 3600.0
    df["potential_evaporation"] = priestley_taylor_pet(
        df["DSWRF_surface"].values, df["TMP_2maboveground"].values)
    return df


def load_obs(obs_file):
    df = pd.read_csv(obs_file)
    date_col = next(c for c in df.columns
                    if c.lower() in ("date", "time", "datetime", "timestamp"))
    q_col    = next(c for c in df.columns
                    if "qkrig" in c.lower() and "var" not in c.lower())
    var_col  = next((c for c in df.columns if "var" in c.lower()), None)
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(date_col)
    dates    = df[date_col].dt.strftime("%Y-%m-%d %H:%M:%S")
    obs_dict = dict(zip(dates, df[q_col]))
    var_dict = dict(zip(dates, df[var_col])) if var_col else {}
    return obs_dict, var_dict


def write_config(cfg_file, forcing_file, best_params, out_path):
    with open(cfg_file) as f:
        cfg = json.load(f)
    cfg["forcing_file"]          = forcing_file
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
    with open(out_path, "w") as f:
        json.dump(cfg, f)
    return out_path


def make_model(bmi_cfe, cfg_path):
    m = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
    m.initialize()
    return m


def step(m, precip_mmh, pet_mmh):
    m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                precip_mmh / 1000.0)
    m.set_value("water_potential_evaporation_flux",
                pet_mmh / 1000.0 / 3600.0)
    m.update()
    return m.get_value("land_surface_water__runoff_depth") * 1000.0


def get_state(m):
    return {"soil_m": float(m.soil_reservoir["storage_m"]),
            "gw_m":   float(m.gw_reservoir["storage_m"]),
            "nash0":  float(m.nash_storage[0]),
            "nash1":  float(m.nash_storage[1])}


def set_state(m, s):
    m.soil_reservoir["storage_m"] = s["soil_m"]
    m.gw_reservoir["storage_m"]   = s["gw_m"]
    m.nash_storage[0]              = s["nash0"]
    m.nash_storage[1]              = s["nash1"]


def perturb_state(s, rng, frac=STATE_FRAC):
    return {"soil_m": max(s["soil_m"] * (1 + frac * rng.standard_normal()), 1e-6),
            "gw_m":   max(s["gw_m"]   * (1 + frac * rng.standard_normal()), 1e-6),
            "nash0":  max(s["nash0"]  + abs(s["nash0"]) * frac * rng.standard_normal(), 0.0),
            "nash1":  max(s["nash1"]  + abs(s["nash1"]) * frac * rng.standard_normal(), 0.0)}


def enkf_r(y_obs, krig_var):
    """F5 R formula: re-kriged variance used directly as observation error variance."""
    return max(float(krig_var), 1e-6)


def enkf_update(state, q_sim, y_obs, krig_var):
    R    = enkf_r(y_obs, krig_var)
    P_yy = max(q_sim * 0.01, 1e-6)
    K    = P_yy / (P_yy + R)
    scale = K * (y_obs - q_sim) / max(abs(q_sim), 1e-6)
    return {"soil_m": max(state["soil_m"] * (1 + scale), 1e-6),
            "gw_m":   max(state["gw_m"]   * (1 + scale), 1e-6),
            "nash0":  max(state["nash0"]  + state["nash0"] * scale, 0.0),
            "nash1":  max(state["nash1"]  + state["nash1"] * scale, 0.0)}


# ── Phase 1: DA spinup + test period to collect warm posterior states ─────────

def run_da_phase(bmi_cfe, cfg_path, df_forcing, obs_dict, var_dict):
    """Run EnKF through spinup + test period.

    Returns: {issue_time Timestamp -> state dict} for all Helene issue times.
    """
    print("  Phase 1: DA spinup + test period...")

    sp_mask   = (df_forcing["date"] >= SPINUP_START) & (df_forcing["date"] <= SPINUP_END)
    test_mask = (df_forcing["date"] >= TEST_START)   & (df_forcing["date"] <= TEST_END)
    df_sp   = df_forcing[sp_mask]
    df_test = df_forcing[test_mask].reset_index(drop=True)

    model = make_model(bmi_cfe, cfg_path)

    # Spinup (no DA, just warm up soil/GW states)
    for _, row in df_sp.iterrows():
        step(model, float(row["total_precipitation"]), float(row["potential_evaporation"]))

    # Test period with EnKF updates
    snapshots = {}
    for _, row in df_test.iterrows():
        date_str = row["date"]
        P  = float(row["total_precipitation"])
        E  = float(row["potential_evaporation"])
        q  = step(model, P, E)
        y  = obs_dict.get(date_str, np.nan)
        if not np.isnan(y):
            kv = float(var_dict.get(date_str, 1e-3))
            set_state(model, enkf_update(get_state(model), q, y, kv))
        ts = pd.Timestamp(date_str)
        if HELENE_START <= ts <= HELENE_END:
            snapshots[ts] = get_state(model)

    model.finalize()
    print(f"  Phase 1 done: {len(snapshots)} Helene issue-time snapshots")
    return snapshots, df_test


# ── Phase 2: 600-member crossed ensemble at each issue time ───────────────────

def run_crossed_phase(bmi_cfe, cfg_path, snapshots, df_test, rng, n_fa, n_ha):
    """For each Helene issue time: generate 600 forecast trajectories.

    Crossing: member(i,j) = forcing_draw_i x hydro_draw_j
      i in [0, n_fa)  -- met perturbation
      j in [0, n_ha)  -- initial state perturbation

    Returns four lists of records for saving to parquet/CSV.
    """
    n_total  = n_fa * n_ha
    date_to_idx = {row["date"]: idx for idx, row in df_test.iterrows()}
    issue_times = sorted(snapshots.keys())

    all_records        = []
    hydro_draw_records = []
    forcing_draw_records = []

    print(f"  Phase 2: {n_fa} x {n_ha} = {n_total} members "
          f"x {len(issue_times)} issue times x {FORECAST_HOURS} leads ...")

    for t_num, t0 in enumerate(issue_times):
        snap    = snapshots[t0]
        t0_str  = t0.strftime("%Y-%m-%d %H:%M:%S")
        if t0_str not in date_to_idx:
            continue
        start_i = date_to_idx[t0_str]

        # -- n_fa forcing perturbation sequences (one per lead hour) --------
        forcing_seqs = []
        for i in range(n_fa):
            seq = []
            for lead in range(1, FORECAST_HOURS + 1):
                fi = start_i + lead
                if fi >= len(df_test):
                    seq.append((0.0, 0.0))
                    forcing_draw_records.append({
                        "issue_time": t0_str, "forcing_draw_i": i,
                        "lead_hour": lead, "P_pert_mmh": 0.0, "E_pert_mmh": 0.0})
                    continue
                row  = df_test.iloc[fi]
                P    = float(row["total_precipitation"])
                E    = float(row["potential_evaporation"])
                mu_p = -0.5 * PRECIP_SIGMA ** 2
                P_p  = float(P * rng.lognormal(mu_p, PRECIP_SIGMA)) if P > 0 else 0.0
                E_p  = max(float(E * (1.0 + PET_SIGMA * rng.standard_normal())), 0.0)
                seq.append((P_p, E_p))
                forcing_draw_records.append({
                    "issue_time": t0_str, "forcing_draw_i": i,
                    "lead_hour": lead, "P_pert_mmh": P_p, "E_pert_mmh": E_p})
            forcing_seqs.append(seq)

        # -- n_ha hydro-state draws (draw 0 = unperturbed DA posterior) -----
        state_draws = [snap if j == 0 else perturb_state(snap, rng)
                       for j in range(n_ha)]
        for j, st in enumerate(state_draws):
            hydro_draw_records.append({
                "issue_time": t0_str, "hydro_draw_j": j, **st})

        # -- Run all 600 combinations ---------------------------------------
        q_matrix = np.full((FORECAST_HOURS, n_fa, n_ha), np.nan, dtype="float32")
        for i, f_seq in enumerate(forcing_seqs):
            for j, st_j in enumerate(state_draws):
                m = make_model(bmi_cfe, cfg_path)
                set_state(m, st_j)
                for lead in range(1, FORECAST_HOURS + 1):
                    P_p, E_p = f_seq[lead - 1]
                    q_matrix[lead - 1, i, j] = step(m, P_p, E_p)
                m.finalize()

        # -- Flatten to rows ------------------------------------------------
        for lead in range(1, FORECAST_HOURS + 1):
            row_data = {"issue_time": t0_str, "lead_hour": lead}
            for mem_idx, (i, j) in enumerate(
                    (i, j) for i in range(n_fa) for j in range(n_ha)):
                row_data[f"member_{mem_idx:04d}"] = float(q_matrix[lead - 1, i, j])
            all_records.append(row_data)

        if (t_num + 1) % 10 == 0 or (t_num + 1) == len(issue_times):
            print(f"  {t_num + 1}/{len(issue_times)} issue times done")

    return all_records, hydro_draw_records, forcing_draw_records


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="F5 20%% holdout — DA + 18h + 600-member crossed ensemble")
    parser.add_argument("--cat-id",       required=True)
    parser.add_argument("--obs-dir",      default=DEFAULT_OBS_DIR)
    parser.add_argument("--param-src",    default=DEFAULT_PARAM_SRC)
    parser.add_argument("--cfe-dir",      default=DEFAULT_CFE_DIR)
    parser.add_argument("--config-file",  default=DEFAULT_CONFIG)
    parser.add_argument("--forcing-dir1", default=DEFAULT_FORCING1)
    parser.add_argument("--forcing-dir2", default=DEFAULT_FORCING2)
    parser.add_argument("--out-dir",      default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-forcing",    type=int, default=N_FORCING)
    parser.add_argument("--n-hydro",      type=int, default=N_HYDRO)
    parser.add_argument("--rng-seed",     type=int, default=None)
    args = parser.parse_args()

    cat_id  = args.cat_id
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{cat_id}] F5 20%% holdout | "
          f"{args.n_forcing} forcing x {args.n_hydro} hydro = "
          f"{args.n_forcing * args.n_hydro} members | 18h lead")

    # ── Load CFE Python package ───────────────────────────────────────────────
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe

    # ── Combine forcing CSVs ──────────────────────────────────────────────────
    combined_path = out_dir / f"{cat_id}_nwm_operational_combined.csv"
    if not combined_path.exists():
        f1 = os.path.join(args.forcing_dir1, f"{cat_id}.csv")
        f2 = os.path.join(args.forcing_dir2, f"{cat_id}.csv")
        combined = (pd.concat([pd.read_csv(f1), pd.read_csv(f2)], ignore_index=True)
                    .drop_duplicates(subset="time").sort_values("time"))
        combined.to_csv(combined_path, index=False)
        print(f"  Combined forcing written: {combined_path}")
    else:
        print(f"  Reusing combined forcing: {combined_path}")

    df_forcing = load_forcing(str(combined_path))

    # ── Load calibrated parameters ────────────────────────────────────────────
    params_path = Path(args.param_src) / cat_id / f"{cat_id}_best_params.json"
    if not params_path.exists():
        raise FileNotFoundError(f"Best params not found: {params_path}")
    with open(params_path) as f:
        best_params = json.load(f)["best_parameters"]
    print(f"  Params loaded from: {params_path}")

    # ── Write BMI config (once, reused for all model instantiations) ──────────
    cfg_stable = str(out_dir / f"{cat_id}_bmi_config_crossed_stable.json")
    write_config(args.config_file, str(combined_path), best_params, cfg_stable)
    print(f"  BMI config written: {cfg_stable}")

    # ── Load Qkrig observations ───────────────────────────────────────────────
    obs_file = os.path.join(args.obs_dir, f"{cat_id}.csv")
    obs_dict, var_dict = load_obs(obs_file)
    print(f"  Obs loaded: {len(obs_dict)} timesteps")

    rng = np.random.default_rng(
        args.rng_seed if args.rng_seed is not None else hash(cat_id) & 0x7fffffff)

    # ── Phase 1: DA warm states ───────────────────────────────────────────────
    snap_cache = out_dir / f"{cat_id}_da_snapshots.parquet"
    if snap_cache.exists():
        print(f"  Loading cached DA snapshots: {snap_cache}")
        df_snaps = pd.read_parquet(snap_cache)
        snapshots = {pd.Timestamp(r["issue_time"]): {
                        "soil_m": r["soil_m"], "gw_m": r["gw_m"],
                        "nash0":  r["nash0"],  "nash1": r["nash1"]}
                     for _, r in df_snaps.iterrows()}
        test_mask = ((df_forcing["date"] >= TEST_START) &
                     (df_forcing["date"] <= TEST_END))
        df_test   = df_forcing[test_mask].reset_index(drop=True)
    else:
        snapshots, df_test = run_da_phase(bmi_cfe, cfg_stable, df_forcing,
                                          obs_dict, var_dict)

    # ── Phase 2: 600-member crossed ensemble ──────────────────────────────────
    all_records, hydro_recs, forcing_recs = run_crossed_phase(
        bmi_cfe, cfg_stable, snapshots, df_test, rng,
        args.n_forcing, args.n_hydro)

    # ── Save outputs ──────────────────────────────────────────────────────────
    n_total = args.n_forcing * args.n_hydro

    # 1. Main ensemble parquet
    df_out   = pd.DataFrame(all_records)
    ens_path = out_dir / f"{cat_id}_crossed_ensemble.parquet"
    df_out.to_parquet(ens_path, index=False)
    print(f"Saved: {ens_path}  shape={df_out.shape}")

    # 2. DA snapshots (warm states at each issue time)
    snap_recs = [{"issue_time": t.strftime("%Y-%m-%d %H:%M:%S"), **s}
                 for t, s in snapshots.items()]
    pd.DataFrame(snap_recs).to_parquet(snap_cache, index=False)
    print(f"Saved: {snap_cache}  ({len(snap_recs)} snapshots)")

    # 3. Member manifest (decoder: member_col -> forcing_i, hydro_j)
    manifest = [{"member": f"member_{k:04d}", "member_idx": k,
                 "forcing_draw_i": k // args.n_hydro,
                 "hydro_draw_j":   k  % args.n_hydro}
                for k in range(n_total)]
    mf_path = out_dir / f"{cat_id}_member_manifest.csv"
    pd.DataFrame(manifest).to_csv(mf_path, index=False)
    print(f"Saved: {mf_path}  ({n_total} members)")

    # 4. Hydro draw initial states
    hd_path = out_dir / f"{cat_id}_hydro_draw_states.parquet"
    pd.DataFrame(hydro_recs).to_parquet(hd_path, index=False)
    print(f"Saved: {hd_path}")

    # 5. Forcing draw sequences
    fd_path = out_dir / f"{cat_id}_forcing_draw_sequences.parquet"
    pd.DataFrame(forcing_recs).to_parquet(fd_path, index=False)
    print(f"Saved: {fd_path}")

    print(f"\n[{cat_id}] Done.")


if __name__ == "__main__":
    main()
