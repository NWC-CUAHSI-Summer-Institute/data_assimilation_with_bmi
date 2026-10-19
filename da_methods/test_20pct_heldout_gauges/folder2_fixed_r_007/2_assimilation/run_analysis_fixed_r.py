"""
run_analysis_fixed_r.py  --  DA analysis trajectory with constant R=0.07.

Runs CFE + EnKF forward through the full test period for one catchment,
using a fixed observation-error variance (R=0.07) instead of the dynamic
Vrugt formula.  Produces a _test_results.csv that run_route.py can read.

Output:
    <out-dir>/<cat-id>/<cat-id>_test_results.csv
    Columns: date, sim_mm_h, obs_mm_h, precip_mm_h

Usage:
    python3 run_analysis_fixed_r.py \\
        --cat-id cat-1016300 \\
        --forcing-file /mnt/disk2/.../cat-1016300_nwm_operational_combined.csv \\
        --obs-dir      /mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance \\
        --params-file  /mnt/disk2/.../cat-1016300_best_params.json \\
        --cfe-dir      /mnt/disk2/suma_helen_poster/cfe_py \\
        --config-file  /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \\
        --out-dir      /mnt/disk2/suma_helen_poster/da_results/v2_fixed_r007_analysis \\
        --hardcoded-r  0.07
"""

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

ALBEDO   = 0.20
ALPHA_PT = 1.26

TEST_START = "2024-08-24 00:00:00"
TEST_END   = "2024-10-31 23:00:00"

CAT_ID          = None
OBS_FILE        = None
CFE_CONFIG_FILE = None
TEST_FORCING_FILE = None
OUT_DIR         = None
HARDCODED_R     = 0.07
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
    tmp_cfg = str(OUT_DIR / f"{CAT_ID}_bmi_config_fixed_r.json")
    with open(tmp_cfg, "w") as f:
        json.dump(cfg, f)
    return tmp_cfg


def run_model_step(m, precip_mmh, pet_mmh):
    m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                precip_mmh / 1000.0)
    m.set_value("water_potential_evaporation_flux",
                pet_mmh / 1000.0 / 3600.0)
    m.update()
    return m.get_value("land_surface_water__runoff_depth") * 1000.0


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


def enkf_update(state, q_sim, y_obs, krig_var):
    R = HARDCODED_R if HARDCODED_R is not None else max(
        (0.10 * max(y_obs, 0.0)) ** 2 + 0.001 * krig_var, 1e-6)
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


def run_da_analysis(cfg_path, df_test, obs_dict, var_dict):
    """Run DA analysis; return per-timestep (date, sim_mm_h, obs_mm_h, precip_mm_h)."""
    model = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
    model.initialize()

    records = []
    for _, row in df_test.iterrows():
        date_str = row["date"]
        P = float(row["total_precipitation"])
        E = float(row["potential_evaporation"])

        q_sim = run_model_step(model, P, E)
        y_obs = obs_dict.get(date_str, np.nan)

        if not np.isnan(y_obs):
            krig_var = var_dict.get(date_str, 0.07)
            state = get_state(model)
            updated = enkf_update(state, q_sim, y_obs, krig_var)
            set_state(model, updated)

        records.append({
            "date":       date_str,
            "sim_mm_h":   q_sim,
            "obs_mm_h":   y_obs,
            "precip_mm_h": P,
        })

    model.finalize()
    return pd.DataFrame(records)


def main():
    global CAT_ID, OBS_FILE, CFE_CONFIG_FILE, TEST_FORCING_FILE
    global OUT_DIR, HARDCODED_R, bmi_cfe

    parser = argparse.ArgumentParser()
    parser.add_argument("--cat-id",       required=True)
    parser.add_argument("--forcing-file", required=True,
                        help="Path to {cat}_nwm_operational_combined.csv")
    parser.add_argument("--obs-dir",      required=True,
                        help="Dir containing {cat}.csv with Qkrig + variance")
    parser.add_argument("--params-file",  required=True,
                        help="Path to {cat}_best_params.json")
    parser.add_argument("--cfe-dir",      required=True,
                        help="Directory containing bmi_cfe.py")
    parser.add_argument("--config-file",  required=True,
                        help="BMI config JSON template")
    parser.add_argument("--out-dir",      required=True)
    parser.add_argument("--hardcoded-r",  type=float, default=0.07,
                        help="Fixed R value for EnKF (default 0.07)")
    args = parser.parse_args()

    CAT_ID          = args.cat_id
    TEST_FORCING_FILE = args.forcing_file
    OBS_FILE        = os.path.join(args.obs_dir, f"{CAT_ID}.csv")
    CFE_CONFIG_FILE = args.config_file
    HARDCODED_R     = args.hardcoded_r
    OUT_DIR         = Path(args.out_dir) / CAT_ID
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    out_csv = OUT_DIR / f"{CAT_ID}_test_results.csv"
    if out_csv.exists():
        print(f"[{CAT_ID}] already done — {out_csv}")
        return

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    with open(args.params_file) as f:
        raw = json.load(f)
    best_params = raw.get("best_parameters", raw)

    obs_dict, var_dict = load_obs()
    df_forcing = load_test_forcing()
    t_mask = (df_forcing["date"] >= TEST_START) & (df_forcing["date"] <= TEST_END)
    df_test = df_forcing[t_mask].reset_index(drop=True)

    print(f"[{CAT_ID}] R={HARDCODED_R}  steps={len(df_test)}")

    cfg_path = write_model_config(best_params)
    df_out = run_da_analysis(cfg_path, df_test, obs_dict, var_dict)

    df_out.to_csv(out_csv, index=False)
    print(f"[{CAT_ID}] saved {out_csv}  ({len(df_out)} rows)")


if __name__ == "__main__":
    main()
