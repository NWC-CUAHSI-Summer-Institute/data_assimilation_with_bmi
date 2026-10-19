"""
run_perturbation_da_on.py  —  2a/2b: ensemble spread with DA ON.

Runs two separated perturbation arms for the Helene window (Sep 24-30),
with DA actively running during the analysis period:

    ARM A — forcing only (30 members):
        DA analyzes states using Qkrig at every hour (DA ON).
        At each issue_time: fork 30 members with ONLY met forcing perturbed
        (different precip/PET draw per member). Hydro states = DA analysis
        mean (same for all 30). Free-run 18 hours, no further DA.

    ARM B — hydro states only (20 members):
        Same DA analysis. At each issue_time: fork 20 members with ONLY
        hydro states perturbed (multiplicative noise on soil/GW/Nash from
        the DA analysis ensemble). Forcing = deterministic. Free-run 18 hours.

These two arms isolate "where does forecast spread come from WITH DA?"
— figures for the assimilation script section.

Output (per catchment, per arm):
    <out-dir>/<cat-id>/<cat-id>_da_forcing_arm.csv
    <out-dir>/<cat-id>/<cat-id>_da_hydro_arm.csv
    Columns: issue_time, lead_hour, member_00..member_29 (forcing arm)
             issue_time, lead_hour, member_00..member_19 (hydro arm)

Usage (run for each catchment, both arms):
    python3 run_perturbation_da_on.py \\
        --cat-id cat-1016300 \\
        --forcing-dir  <NWM retro dir> \\
        --obs-dir      /mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance \\
        --cfe-dir      /mnt/disk2/suma_helen_poster/cfe_py \\
        --config-file  /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \\
        --param-bounds /mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json \\
        --out-dir      /mnt/disk2/suma_helen_poster/da_results/v2_perturbation_da_on \\
        --test-forcing-dir1 <NWM op 2023-Feb2024> \\
        --test-forcing-dir2 <NWM op Feb2024-Sep2025> \\
        --hardcoded-r 0.07
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

ALBEDO    = 0.20
ALPHA_PT  = 1.26

N_FORCING = 30   # members in forcing arm
N_HYDRO   = 20   # members in hydro-state arm

# Analysis period — run DA through entire test period
TEST_START = "2024-08-24 00:00:00"   # 1 month spinup before Helene window
TEST_END   = "2024-10-31 23:00:00"

# Issue times for the Helene window (hourly, Sep 24-30)
HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-30 06:00:00")
FORECAST_HOURS = 18

# Perturbation magnitudes (same as production)
PRECIP_SIGMA  = 0.15
PET_SIGMA     = 0.10
STATE_FRAC    = 0.05   # hydro-state perturbation fraction

# Set at runtime
CAT_ID          = None
OBS_FILE        = None
CFE_CONFIG_FILE = None
TEST_FORCING_FILE = None
OUT_DIR         = None
HARDCODED_R     = None
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
    if date_col is None:
        raise ValueError(f"No date column in {OBS_FILE}. Columns: {df.columns.tolist()}")
    q_col = next((c for c in df.columns if "qkrig" in c.lower() and "var" not in c.lower()), None)
    if q_col is None:
        raise ValueError(f"No qkrig column in {OBS_FILE}. Columns: {df.columns.tolist()}")
    var_col = next((c for c in df.columns if "var" in c.lower()), None)
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(date_col).reset_index(drop=True)
    date_strs = df[date_col].dt.strftime("%Y-%m-%d %H:%M:%S")
    obs_dict = dict(zip(date_strs, df[q_col]))
    var_vals = df[var_col] if var_col else pd.Series([0.07] * len(df))
    var_dict = dict(zip(date_strs, var_vals))
    print(f"  Obs: date='{date_col}', q='{q_col}', var='{var_col}'")
    return obs_dict, var_dict


def make_model(best_params):
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
    tmp_cfg = str(OUT_DIR / f"{CAT_ID}_bmi_config_temp_da_on.json")
    with open(tmp_cfg, "w") as f:
        json.dump(cfg, f)
    m = bmi_cfe.BMI_CFE(cfg_file=tmp_cfg)
    m.initialize()
    return m


def get_state(m):
    return {
        "soil_m":  m.soil_reservoir["storage_m"],
        "gw_m":    m.gw_reservoir["storage_m"],
        "nash0":   float(m.nash_storage[0]),
        "nash1":   float(m.nash_storage[1]),
    }


def set_state(m, state):
    m.soil_reservoir["storage_m"] = state["soil_m"]
    m.gw_reservoir["storage_m"]   = state["gw_m"]
    m.nash_storage[0]              = state["nash0"]
    m.nash_storage[1]              = state["nash1"]


def perturb_state(state, rng, frac=STATE_FRAC):
    sm_max = 1.0   # placeholder — will be clipped in set_state; actual max enforced by model
    out = {
        "soil_m": max(state["soil_m"] * (1.0 + frac * rng.standard_normal()), 1e-6),
        "gw_m":   max(state["gw_m"]   * (1.0 + frac * rng.standard_normal()), 1e-6),
        "nash0":  max(state["nash0"]  * (1.0 + frac * rng.standard_normal()), 0.0),
        "nash1":  max(state["nash1"]  * (1.0 + frac * rng.standard_normal()), 0.0),
    }
    return out


def run_model_step(m, precip_mmh, pet_mmh):
    m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                precip_mmh / 1000.0)
    m.set_value("water_potential_evaporation_flux",
                pet_mmh / 1000.0 / 3600.0)
    m.update()
    return m.get_value("land_surface_water__runoff_depth") * 1000.0  # m/h -> mm/h


def enkf_update(state, q_sim, y_obs, krig_var, rng, n_members=1):
    """Simple scalar EnKF update on the mean analysis state."""
    if HARDCODED_R is not None:
        R = HARDCODED_R
    else:
        alpha = 0.10
        R = (alpha * max(y_obs, 0.0)) ** 2 + 0.001 * krig_var
    R = max(R, 1e-6)
    P_yy = max(q_sim * 0.01, 1e-6)   # rough prior variance on Q
    K = P_yy / (P_yy + R)
    innovation = y_obs - q_sim
    # Apply correction to each state proportionally (simple scalar gain)
    scale = K * innovation / max(abs(q_sim), 1e-6)
    updated = {
        "soil_m": max(state["soil_m"] * (1.0 + scale), 1e-6),
        "gw_m":   max(state["gw_m"]   * (1.0 + scale), 1e-6),
        "nash0":  max(state["nash0"]  + state["nash0"] * scale, 0.0),
        "nash1":  max(state["nash1"]  + state["nash1"] * scale, 0.0),
    }
    return updated


def main():
    global CAT_ID, OBS_FILE, CFE_CONFIG_FILE, TEST_FORCING_FILE
    global OUT_DIR, HARDCODED_R, bmi_cfe

    parser = argparse.ArgumentParser()
    parser.add_argument("--cat-id",            required=True)
    parser.add_argument("--forcing-dir",       required=True)
    parser.add_argument("--obs-dir",           required=True)
    parser.add_argument("--cfe-dir",           required=True)
    parser.add_argument("--config-file",       required=True)
    parser.add_argument("--param-bounds",      required=True)
    parser.add_argument("--out-dir",           required=True)
    parser.add_argument("--test-forcing-dir1", default=None)
    parser.add_argument("--test-forcing-dir2", default=None)
    parser.add_argument("--hardcoded-r",       type=float, default=None)
    args = parser.parse_args()

    CAT_ID          = args.cat_id
    OBS_FILE        = os.path.join(args.obs_dir, f"{CAT_ID}.csv")
    CFE_CONFIG_FILE = args.config_file
    HARDCODED_R     = args.hardcoded_r
    OUT_DIR         = Path(args.out_dir) / CAT_ID
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    # Build combined test forcing
    if args.test_forcing_dir1 and args.test_forcing_dir2:
        f1 = os.path.join(args.test_forcing_dir1, f"{CAT_ID}.csv")
        f2 = os.path.join(args.test_forcing_dir2, f"{CAT_ID}.csv")
        df1, df2 = pd.read_csv(f1), pd.read_csv(f2)
        combined = (pd.concat([df1, df2], ignore_index=True)
                    .drop_duplicates(subset="time").sort_values("time"))
        TEST_FORCING_FILE = str(OUT_DIR / f"{CAT_ID}_nwm_operational_combined.csv")
        combined.to_csv(TEST_FORCING_FILE, index=False)
    else:
        print("Need --test-forcing-dir1 and --test-forcing-dir2"); return

    best_params_file = OUT_DIR / f"{CAT_ID}_best_params.json"
    if not best_params_file.exists():
        print(f"Missing best params at {best_params_file}"); return
    with open(best_params_file) as f:
        best_params = json.load(f)["best_parameters"]

    obs_dict, var_dict = load_obs()
    df_forcing = load_test_forcing()
    t_mask = (df_forcing["date"] >= TEST_START) & (df_forcing["date"] <= TEST_END)
    df_test = df_forcing[t_mask].reset_index(drop=True)

    rng = np.random.default_rng(hash(CAT_ID) & 0x7fffffff)

    # ------------------------------------------------------------------ #
    # Phase 1: Run DA through test period, snapshot analysis state at     #
    # each Helene issue time                                               #
    # ------------------------------------------------------------------ #
    print(f"[{CAT_ID}] Phase 1: DA analysis through test period...")
    model = make_model(best_params)

    helene_snapshots = {}   # issue_time -> analysis state dict

    for h, row in df_test.iterrows():
        date_str = row["date"]
        P = float(row["total_precipitation"])
        E = float(row["potential_evaporation"])
        q_sim = run_model_step(model, P, E)

        date_ts = pd.Timestamp(date_str)

        # DA update if obs available
        y_obs = obs_dict.get(date_str, np.nan)
        if not np.isnan(y_obs):
            krig_var = var_dict.get(date_str, 0.07)
            state = get_state(model)
            updated = enkf_update(state, q_sim, y_obs, krig_var, rng)
            set_state(model, updated)

        # Snapshot if this is a Helene issue time
        if HELENE_START <= date_ts <= HELENE_END:
            helene_snapshots[date_ts] = get_state(model)

    model.finalize()
    print(f"  Snapshots saved: {len(helene_snapshots)} issue times")

    # ------------------------------------------------------------------ #
    # Phase 2A: Forcing arm — 30 members, fixed analyzed state, varying   #
    # met forcing                                                          #
    # ------------------------------------------------------------------ #
    print(f"[{CAT_ID}] Phase 2A: Forcing arm ({N_FORCING} members)...")
    forcing_records = []

    issue_times = sorted(helene_snapshots.keys())
    t0_to_idx = {t: i for i, t in enumerate(df_test["date"].apply(pd.Timestamp))}

    for t0 in issue_times:
        snap = helene_snapshots[t0]
        t0_str = t0.strftime("%Y-%m-%d %H:%M:%S")

        # Find index in df_test
        if t0 not in t0_to_idx:
            continue
        start_idx = t0_to_idx[t0]

        q_matrix = np.full((FORECAST_HOURS, N_FORCING), np.nan)

        for mem in range(N_FORCING):
            m = make_model(best_params)
            set_state(m, snap)   # all members start from same DA analysis state

            for lead in range(1, FORECAST_HOURS + 1):
                fi = start_idx + lead
                if fi >= len(df_test):
                    break
                row = df_test.iloc[fi]
                P = float(row["total_precipitation"])
                E = float(row["potential_evaporation"])
                # Perturb forcing only
                mu_p = -0.5 * PRECIP_SIGMA ** 2
                P_pert = float(P * rng.lognormal(mu_p, PRECIP_SIGMA)) if P > 0 else 0.0
                E_pert = max(float(E * (1.0 + PET_SIGMA * rng.standard_normal())), 0.0)
                q = run_model_step(m, P_pert, E_pert)
                q_matrix[lead - 1, mem] = q

            m.finalize()

        row_data = {"issue_time": t0_str}
        for lead in range(1, FORECAST_HOURS + 1):
            row_data["lead_hour"] = lead
            for mem in range(N_FORCING):
                row_data[f"member_{mem:02d}"] = q_matrix[lead - 1, mem]
            forcing_records.append(dict(row_data))

    df_forcing_arm = pd.DataFrame(forcing_records)
    out_a = OUT_DIR / f"{CAT_ID}_da_forcing_arm.csv"
    df_forcing_arm.to_csv(out_a, index=False)
    print(f"  Saved: {out_a}")

    # ------------------------------------------------------------------ #
    # Phase 2B: Hydro-state arm — 20 members, deterministic forcing,      #
    # perturbed analysis states                                            #
    # ------------------------------------------------------------------ #
    print(f"[{CAT_ID}] Phase 2B: Hydro-state arm ({N_HYDRO} members)...")
    hydro_records = []

    for t0 in issue_times:
        snap = helene_snapshots[t0]
        t0_str = t0.strftime("%Y-%m-%d %H:%M:%S")

        if t0 not in t0_to_idx:
            continue
        start_idx = t0_to_idx[t0]

        q_matrix = np.full((FORECAST_HOURS, N_HYDRO), np.nan)

        for mem in range(N_HYDRO):
            # Perturb the DA analysis state for this member
            if mem == 0:
                state_m = snap   # member 0 = unperturbed analysis mean
            else:
                state_m = perturb_state(snap, rng)

            m = make_model(best_params)
            set_state(m, state_m)

            for lead in range(1, FORECAST_HOURS + 1):
                fi = start_idx + lead
                if fi >= len(df_test):
                    break
                row = df_test.iloc[fi]
                P = float(row["total_precipitation"])
                E = float(row["potential_evaporation"])
                # Deterministic forcing
                q = run_model_step(m, P, E)
                q_matrix[lead - 1, mem] = q

            m.finalize()

        row_data = {"issue_time": t0_str}
        for lead in range(1, FORECAST_HOURS + 1):
            row_data["lead_hour"] = lead
            for mem in range(N_HYDRO):
                row_data[f"member_{mem:02d}"] = q_matrix[lead - 1, mem]
            hydro_records.append(dict(row_data))

    df_hydro_arm = pd.DataFrame(hydro_records)
    out_b = OUT_DIR / f"{CAT_ID}_da_hydro_arm.csv"
    df_hydro_arm.to_csv(out_b, index=False)
    print(f"  Saved: {out_b}")
    print(f"[{CAT_ID}] Done.")


if __name__ == "__main__":
    main()
