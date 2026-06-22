#!/usr/bin/env python3
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless -- must be set before bmi_cfe imports pyplot
import matplotlib.pyplot as plt

# bmi_cfe lives in --cfe-dir; imported once in main() after sys.path is set.
bmi_cfe = None

# ---------------------------------------------------------------------------
# Experiment constants -- 20% held-out (fold-0) configuration
# ---------------------------------------------------------------------------
DEFAULT_OBS_DIR     = "/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig"
DEFAULT_PARAMS_DIR  = "/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded"
DEFAULT_CFE_DIR     = "/mnt/disk2/suma_helen_poster/cfe_py"
DEFAULT_CONFIG_FILE = "/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
DEFAULT_FORCING_DIR1 = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                        "output_03463300_nwmoperational/03463300/2023_2024_feb/forcings")
DEFAULT_FORCING_DIR2 = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                        "output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings")
DEFAULT_OUT_DIR     = "/mnt/disk2/1400_sites_helene/da_results_enkf600_rekrig_direct"

SPINUP_START = "2023-02-01 00:00:00"   # spread develops here, no DA
SPINUP_END   = "2023-09-30 23:00:00"
TEST_START   = "2023-10-01 00:00:00"   # assimilation window (includes Helene)
TEST_END     = "2024-10-31 23:00:00"
HELENE_START = pd.Timestamp("2024-09-24 00:00:00")   # member states saved here
HELENE_END   = pd.Timestamp("2024-09-30 06:00:00")

# Perturbation magnitudes -- WIDENED for HEFS-tier ensemble spread (deliberately
# larger than calibrate_catchment_cfe_da_v2.py / the NWC-CUAHSI reference). The
# dominant spread driver is the persistent per-member precip bias (PRECIP_BIAS_SIGMA).
PRECIP_SIGMA      = 0.50    # lognormal sigma on the HOURLY precip multiplier (mean 1; ~0.3x-3x)
PRECIP_BIAS_SIGMA = 0.40    # lognormal sigma on the PERSISTENT per-member precip bias (drawn once)
PET_SIGMA         = 0.20    # Gaussian sigma on PET multiplier
INIT_STATE_SIGMA  = 0.20    # Gaussian sigma on initial soil/GW storages
SOIL_NOISE        = 0.010   # per-hour process noise: soil (slow store)
GW_NOISE          = 0.008   # per-hour process noise: groundwater (slowest store)
NASH_NOISE        = 0.020   # per-hour process noise: Nash buckets (need extra spread)
R_FLOOR           = 1e-6    # lower bound on R so the Kalman gain never divides by ~0

# Priestley-Taylor PET constants
ALBEDO   = 0.20
ALPHA_PT = 1.26


# ---------------------------------------------------------------------------
# 1. LOAD -- forcing, observations, parameters, BMI config
# ---------------------------------------------------------------------------
def priestley_taylor_pet(srad_wm2, T_kelvin):
    """Hourly Priestley-Taylor PET (mm/h) from shortwave radiation + air temp."""
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    pet = ALPHA_PT * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def build_combined_forcing(cat_id, dir1, dir2, cat_out_dir):
    """Concatenate the two NWM-operational forcing CSVs into one file (reused if present)."""
    out_csv = cat_out_dir / f"{cat_id}_nwm_operational_combined.csv"
    if out_csv.exists():
        print(f"  Reusing combined forcing: {out_csv}")
        return out_csv
    df1 = pd.read_csv(Path(dir1) / f"{cat_id}.csv")
    df2 = pd.read_csv(Path(dir2) / f"{cat_id}.csv")
    combined = (pd.concat([df1, df2], ignore_index=True)
                .drop_duplicates(subset="time").sort_values("time"))
    combined.to_csv(out_csv, index=False)
    print(f"  Wrote combined forcing: {out_csv}")
    return out_csv


def load_forcing(path):
    """NWM operational CSV -> DataFrame with date (str), precip (mm/h), PET (mm/h)."""
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    df["total_precipitation"] = df["APCP_surface"] * 3600.0   # kg/m2/s -> mm/h
    df["potential_evaporation"] = priestley_taylor_pet(
        df["DSWRF_surface"].values, df["TMP_2maboveground"].values)
    return df


def load_observations(obs_csv):
    """qkrig obs CSV -> two dicts keyed by 'YYYY-MM-DD HH:MM:SS':
    observed flow (mm/h) and per-hour re-kriged variance ((mm/h)^2, used as R)."""
    df = pd.read_csv(obs_csv)
    for col in ("time", "qkrig_mm_hr", "qkrig_variance"):
        if col not in df.columns:
            raise ValueError(f"{obs_csv} is missing column '{col}' "
                             f"(found {list(df.columns)})")
    dates = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    obs_dict = dict(zip(dates, df["qkrig_mm_hr"].astype(float)))
    var_dict = dict(zip(dates, df["qkrig_variance"].astype(float)))
    return obs_dict, var_dict


def write_bmi_config(template_file, best_params, forcing_csv, out_path):
    """Stamp the calibrated parameters into the BMI config template; written once,
    shared by all members."""
    with open(template_file) as f:
        cfg = json.load(f)
    cfg["forcing_file"]          = str(forcing_csv)
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


# ---------------------------------------------------------------------------
# 2. BUILD -- the ensemble
# ---------------------------------------------------------------------------
def build_ensemble(n_members, cfg_path, forcing_df, rng):
    """Create n CFE members. Member 0 starts exactly at the config's initial
    state; members 1..n-1 get small random perturbations on the initial
    storages so the ensemble has spread from hour zero."""
    bmi_forcing = forcing_df.rename(columns={"date": "time"})
    models = []
    for i in range(n_members):
        m = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
        # All members share ONE in-memory forcing table (read-only). Without
        # this every member would re-read the CSV: slow and, at 600 members,
        # several GB of duplicated DataFrames.
        m.load_forcing_file = (lambda self=m: setattr(self, "forcing_data", bmi_forcing))
        m.initialize()
        if i > 0:
            # Direct attribute writes: BMI exposes only SOIL_CONCEPTUAL_STORAGE,
            # so GW/Nash must be set this way. Bound-clipping here may discard a
            # sliver of init noise; both choices deliberately mirror
            # calibrate_catchment_cfe_da_v2.py so results stay comparable.
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]
            m.soil_reservoir["storage_m"] = float(np.clip(
                m.soil_reservoir["storage_m"]
                * (1.0 + INIT_STATE_SIGMA * rng.standard_normal()), 1e-6, sm_max))
            m.gw_reservoir["storage_m"] = float(np.clip(
                m.gw_reservoir["storage_m"]
                * (1.0 + INIT_STATE_SIGMA * rng.standard_normal()), 1e-6, gw_max))
            # Nash buckets start at 0 -> multiplicative noise would do nothing;
            # use a tiny additive jitter instead.
            jitter = max(sm_max * 1e-6, 1e-9)
            m.nash_storage[0] = max(float(m.nash_storage[0])
                                    + jitter * rng.standard_normal(), 0.0)
            m.nash_storage[1] = max(float(m.nash_storage[1])
                                    + jitter * rng.standard_normal(), 0.0)
        models.append(m)
    return models


# ---------------------------------------------------------------------------
# 3./4. STEP, NOISE, and the EnKF ANALYSIS itself
# ---------------------------------------------------------------------------
def step_all_members(models, precip_mmh, pet_mmh, precip_bias, rng):
    """Advance every member one hour, each with its own perturbed forcing.
    Returns the per-member simulated runoff (mm/h)."""
    n = len(models)
    p_members = precip_mmh * precip_bias * rng.lognormal(-0.5 * PRECIP_SIGMA ** 2, PRECIP_SIGMA, n)
    e_members = np.maximum(pet_mmh * (1.0 + PET_SIGMA * rng.standard_normal(n)), 0.0)
    q_mmh = np.empty(n)
    for i, m in enumerate(models):
        m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                    float(p_members[i]) / 1000.0)            # mm/h -> m/h
        m.set_value("water_potential_evaporation_flux",
                    float(e_members[i]) / 1000.0 / 3600.0)   # mm/h -> m/s
        m.update()
        q_mmh[i] = m.get_value("land_surface_water__runoff_depth") * 1000.0  # m/h -> mm/h
    return q_mmh


def add_process_noise(models, rng):
    """Per-hour state perturbation on every member (the EnKF 'model error' / Q
    term). Without it the analysis pulls members together, spread collapses,
    Pyy -> 0 and the filter goes blind. Per-state sigmas follow each store's
    response time; Nash also gets a small additive floor so noise survives
    when a bucket is empty."""
    for m in models:
        sm_max = m.soil_reservoir["storage_max_m"]
        gw_max = m.gw_reservoir["storage_max_m"]
        nash_floor = max(sm_max * 1e-4, 1e-7)
        m.soil_reservoir["storage_m"] = float(np.clip(
            m.soil_reservoir["storage_m"]
            * (1.0 + SOIL_NOISE * rng.standard_normal()), 0.0, sm_max))
        m.gw_reservoir["storage_m"] = float(np.clip(
            m.gw_reservoir["storage_m"]
            * (1.0 + GW_NOISE * rng.standard_normal()), 0.0, gw_max))
        m.nash_storage[0] = max(
            float(m.nash_storage[0]) * (1.0 + NASH_NOISE * rng.standard_normal())
            + nash_floor * rng.standard_normal(), 0.0)
        m.nash_storage[1] = max(
            float(m.nash_storage[1]) * (1.0 + NASH_NOISE * rng.standard_normal())
            + nash_floor * rng.standard_normal(), 0.0)


def enkf_update(models, q_forecast_mmh, y_obs, krig_var, rng, diag):
    """One stochastic-EnKF analysis step (Burgers / van Leeuwen / Evensen 1998).

    - Pyy   = ensemble variance of simulated runoff (no heuristics)
    - Pxy   = cross-covariance of each state with runoff -> one Kalman gain
              per state; the data decides how the correction splits
    - R     = re-kriged per-hour kriging variance, used directly
    - the observation is perturbed once per member (textbook stochastic EnKF)
    - corrections pass through a mass-conserving cascade instead of being
      silently clipped at the physical bounds

    Returns True if an update was applied; accumulates diagnostics in `diag`.
    """
    n = len(models)
    q = np.asarray(q_forecast_mmh, dtype=float)
    if np.any(np.isnan(q)) or n < 2:
        return False
    R = max(float(krig_var), R_FLOOR)

    # Current state ensemble (meters of storage)
    sm = np.array([m.get_value("SOIL_CONCEPTUAL_STORAGE") for m in models], dtype=float)
    gw = np.array([m.gw_reservoir["storage_m"] for m in models], dtype=float)
    n0 = np.array([float(m.nash_storage[0]) for m in models], dtype=float)
    n1 = np.array([float(m.nash_storage[1]) for m in models], dtype=float)

    q_anom = q - q.mean()
    Pyy = float(q_anom @ q_anom) / (n - 1)
    denom = Pyy + R
    if denom < 1e-12:
        return False   # degenerate: no ensemble spread and no obs error

    # Kalman gain per state: K = Pxy / (Pyy + R)   [units: m per (mm/h)]
    K_sm = float((sm - sm.mean()) @ q_anom) / (n - 1) / denom
    K_gw = float((gw - gw.mean()) @ q_anom) / (n - 1) / denom
    K_n0 = float((n0 - n0.mean()) @ q_anom) / (n - 1) / denom
    K_n1 = float((n1 - n1.mean()) @ q_anom) / (n - 1) / denom

    # Each member is updated against its own perturbed copy of the observation
    obs_pert = y_obs + np.sqrt(R) * rng.standard_normal(n)
    innov = obs_pert - q   # innovation (mm/h), one per member

    for i, m in enumerate(models):
        sm_max = m.soil_reservoir["storage_max_m"]
        gw_max = m.gw_reservoir["storage_max_m"]

        sm_raw = sm[i] + K_sm * innov[i]
        gw_raw = gw[i] + K_gw * innov[i]
        n0_raw = n0[i] + K_n0 * innov[i]
        n1_raw = n1[i] + K_n1 * innov[i]

        # Mass-conserving cascade along CFE's flow path.
        # Positive overflow: soil -> Nash[0], GW -> Nash[1]
        so = max(sm_raw - sm_max, 0.0); sm_raw -= so; n0_raw += so
        go = max(gw_raw - gw_max, 0.0); gw_raw -= go; n1_raw += go
        # Negative underflow: Nash[1] -> Nash[0] -> soil -> GW -> (lost, logged)
        n1u = max(-n1_raw, 0.0); n1_raw += n1u; n0_raw -= n1u
        n0u = max(-n0_raw, 0.0); n0_raw += n0u; sm_raw -= n0u
        smu = max(-sm_raw, 0.0); sm_raw += smu; gw_raw -= smu
        gwu = max(-gw_raw, 0.0); gw_raw += gwu
        diag["mass_lost_mm"] += gwu * 1000.0

        if np.isnan(sm_raw) or np.isnan(gw_raw) or np.isnan(n0_raw) or np.isnan(n1_raw):
            continue

        # Soil goes through BMI; GW/Nash are not BMI-exposed, so these are
        # deliberate direct-attribute writes. The clips are belt-and-braces:
        # the cascade above already kept every value inside its bounds.
        m.set_value("SOIL_CONCEPTUAL_STORAGE", float(np.clip(sm_raw, 0.0, sm_max)))
        m.gw_reservoir["storage_m"] = float(np.clip(gw_raw, 0.0, gw_max))
        m.nash_storage[0] = float(max(n0_raw, 0.0))
        m.nash_storage[1] = float(max(n1_raw, 0.0))

    # Diagnostics for the end-of-run summary
    diag["n_updates"]    += 1
    diag["abs_innov_sum"] += abs(float(innov.mean()))
    diag["pyy_sum"]      += Pyy
    diag["K_sm_sum"]     += K_sm
    diag["K_gw_sum"]     += K_gw
    diag["K_n0_sum"]     += K_n0
    diag["K_n1_sum"]     += K_n1

    if Pyy > 1e-4 or diag["n_updates"] % 500 == 0:
        print(f"  [EnKF N={n}] update {diag['n_updates']} | obs={y_obs:.3f} "
              f"q_mean={q.mean():.3f} innov={float(innov.mean()):+.4f} "
              f"R={R:.2e} Pyy={Pyy:.2e} K_sm={K_sm:+.2e}")
    return True


# ---------------------------------------------------------------------------
# 5. SAVE -- scores and plot
# ---------------------------------------------------------------------------
def kge(obs, sim):
    """Kling-Gupta Efficiency on the overlap of obs and sim (NaNs skipped)."""
    ok = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[ok], sim[ok]
    if len(o) < 2:
        return float("nan")
    denom = float(np.sqrt(((o - o.mean()) ** 2).sum() * ((s - s.mean()) ** 2).sum()))
    if denom == 0 or o.std() == 0 or o.mean() == 0:
        return float("nan")
    r = float(((o - o.mean()) * (s - s.mean())).sum() / denom)
    alpha = float(s.std() / o.std())
    beta = float(s.mean() / o.mean())
    return 1.0 - float(np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def nse(obs, sim):
    """Nash-Sutcliffe Efficiency on the overlap of obs and sim (NaNs skipped)."""
    ok = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[ok], sim[ok]
    if len(o) < 2:
        return float("nan")
    denom = ((o - o.mean()) ** 2).sum()
    if denom == 0:
        return float("nan")
    return 1.0 - float(((o - s) ** 2).sum() / denom)


def plot_results(df, cat_id, n_members, kge_val, nse_val, out_png):
    """Two panels: full test year and Helene zoom, ensemble mean +/- 2 std."""
    dates = pd.to_datetime(df["date"])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10))
    for ax, lo, hi, title in (
        (ax1, dates.iloc[0], dates.iloc[-1],
         f"{cat_id} | EnKF N={n_members} | KGE={kge_val:.3f} NSE={nse_val:.3f}"),
        (ax2, pd.Timestamp("2024-09-20"), pd.Timestamp("2024-10-05"),
         "Helene zoom: Sep 20 - Oct 5, 2024"),
    ):
        w = (dates >= lo) & (dates <= hi)
        ax.fill_between(dates[w],
                        (df["sim_mm_h"] - 2 * df["sim_std_mm_h"])[w].clip(lower=0),
                        (df["sim_mm_h"] + 2 * df["sim_std_mm_h"])[w],
                        color="tomato", alpha=0.25, label="ensemble +/- 2 std")
        ax.plot(dates[w], df["sim_mm_h"][w], "tomato", lw=1.5, label="ensemble mean")
        ax.plot(dates[w], df["obs_mm_h"][w], "k", lw=1.0, label="qkrig obs")
        ax.set_ylabel("Discharge (mm/h)")
        ax.set_title(title)
        ax.legend(loc="upper left")
        twin = ax.twinx()
        twin.plot(dates[w], df["precip_mm_h"][w], "steelblue", lw=0.8, alpha=0.5)
        twin.set_ylim(twin.get_ylim()[::-1])   # rain hangs from the top
        twin.set_ylabel("Precip (mm/h)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# The whole story for one catchment
# ---------------------------------------------------------------------------
def run_catchment(cat_id, args):
    t_start = time.time()
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'=' * 70}\n{cat_id}: EnKF N={args.n_members}, seed={args.rng_seed}, "
          f"R = re-kriged variance (direct)\n{'=' * 70}")

    # ---- 1. LOAD -------------------------------------------------------
    params_file = Path(args.params_dir) / cat_id / f"{cat_id}_best_params.json"
    with open(params_file) as f:
        best_params = json.load(f)["best_parameters"]
    print(f"  Loaded calibrated parameters: {params_file}")

    forcing_csv = build_combined_forcing(
        cat_id, args.test_forcing_dir1, args.test_forcing_dir2, out_dir)
    forcing = load_forcing(forcing_csv)
    obs_dict, var_dict = load_observations(Path(args.obs_dir) / f"{cat_id}.csv")

    cfg_path = write_bmi_config(args.config_file, best_params, forcing_csv,
                                out_dir / f"{cat_id}_bmi_config_enkf.json")

    spin = forcing[(forcing["date"] >= SPINUP_START)
                   & (forcing["date"] <= SPINUP_END)].reset_index(drop=True)
    test = forcing[(forcing["date"] >= TEST_START)
                   & (forcing["date"] <= TEST_END)].reset_index(drop=True)
    print(f"  Spinup {SPINUP_START[:10]} -> {SPINUP_END[:10]} ({len(spin)} h) | "
          f"assimilation {TEST_START[:10]} -> {TEST_END[:10]} ({len(test)} h)")

    # ---- 2. BUILD ------------------------------------------------------
    rng = np.random.default_rng(args.rng_seed)
    print(f"  Building {args.n_members} ensemble members...")
    models = build_ensemble(args.n_members, str(cfg_path), forcing, rng)

    # Persistent per-member precip bias b_i -- drawn ONCE here at spinup start
    # (lognormal, mean 1) and held fixed for the whole run. This is the dominant
    # source of ensemble spread: each member sees a systematically wetter/drier
    # rainfall total all year, not just hourly noise that averages out.
    precip_bias = rng.lognormal(-0.5 * PRECIP_BIAS_SIGMA ** 2,
                                PRECIP_BIAS_SIGMA, args.n_members)

    # ---- 3. SPINUP (no DA) ----------------------------------------------
    print("  Spinup...")
    for k, row in enumerate(spin.itertuples()):
        step_all_members(models, float(row.total_precipitation),
                         float(row.potential_evaporation), precip_bias, rng)
        add_process_noise(models, rng)
        if (k + 1) % 1000 == 0:
            print(f"    spinup {k + 1}/{len(spin)} h ({time.time() - t_start:.0f}s)")

    # ---- 4. ASSIMILATE ---------------------------------------------------
    print("  Assimilating...")
    diag = {"n_updates": 0, "abs_innov_sum": 0.0, "pyy_sum": 0.0, "K_sm_sum": 0.0,
            "K_gw_sum": 0.0, "K_n0_sum": 0.0, "K_n1_sum": 0.0, "mass_lost_mm": 0.0}
    records = []
    member_state_rows = []
    member_q_rows = [] if args.save_members else None  # per-hour (n_members,) forecast runoff

    for k, row in enumerate(test.itertuples()):
        date = row.date
        precip = float(row.total_precipitation)
        pet = float(row.potential_evaporation)

        # Forecast step. q is recorded BEFORE the analysis below, so the
        # output series is the one-hour-ahead forecast -- comparing it with
        # obs(t) is standard sequential-filter verification, not circular.
        q = step_all_members(models, precip, pet, precip_bias, rng)
        if member_q_rows is not None:
            member_q_rows.append(q.astype(np.float32))

        y_obs = obs_dict.get(date, np.nan)
        krig_var = var_dict.get(date, np.nan)
        updated = False
        if not (args.no_da or np.isnan(y_obs) or np.isnan(krig_var)):
            updated = enkf_update(models, q, y_obs, krig_var, rng, diag)

        # Process noise on the analyzed states keeps spread alive for the
        # next cycle (and during obs gaps).
        add_process_noise(models, rng)

        # Save every member's carried-forward state through the Helene
        # window -- these are the restart points the forecast step will use.
        ts = pd.Timestamp(date)
        if HELENE_START <= ts <= HELENE_END:
            for i, m in enumerate(models):
                member_state_rows.append({
                    "time": date, "member": i,
                    "soil_m": float(m.soil_reservoir["storage_m"]),
                    "gw_m":   float(m.gw_reservoir["storage_m"]),
                    "nash0":  float(m.nash_storage[0]),
                    "nash1":  float(m.nash_storage[1]),
                })

        records.append({
            "date": date,
            "sim_mm_h": float(q.mean()),
            "sim_std_mm_h": float(q.std(ddof=1)),
            "obs_mm_h": float(y_obs) if not np.isnan(y_obs) else np.nan,
            "R_var": float(krig_var) if not np.isnan(krig_var) else np.nan,
            "precip_mm_h": precip,
            "da_updated": int(updated),
        })
        if (k + 1) % 500 == 0:
            print(f"    assim {k + 1}/{len(test)} h | updates so far: "
                  f"{diag['n_updates']} ({time.time() - t_start:.0f}s)")

    for m in models:
        m.finalize()

    # ---- 5. SAVE ---------------------------------------------------------
    df = pd.DataFrame(records)
    kge_val = kge(df["obs_mm_h"].values, df["sim_mm_h"].values)
    nse_val = nse(df["obs_mm_h"].values, df["sim_mm_h"].values)

    results_csv = out_dir / f"{cat_id}_enkf_results.csv"
    df.to_csv(results_csv, index=False)

    states_pq = out_dir / f"{cat_id}_member_states_helene.parquet"
    pd.DataFrame(member_state_rows).to_parquet(states_pq, index=False)

    if member_q_rows is not None:
        # Same pre-analysis (one-hour-ahead forecast) series as sim_mm_h, so
        # per-member skill computed from this file is directly comparable to
        # the ensemble-mean KGE/NSE above.
        member_q = pd.DataFrame(
            np.vstack(member_q_rows),
            columns=[f"member_{i:04d}" for i in range(args.n_members)])
        member_q.insert(0, "time", df["date"].values)
        member_q_pq = out_dir / f"{cat_id}_member_q.parquet"
        member_q.to_parquet(member_q_pq, index=False)
        print(f"  Saved: {member_q_pq.name} "
              f"({len(member_q)} h x {args.n_members} members)")

    n_up = max(diag["n_updates"], 1)
    summary = {
        "cat_id": cat_id,
        "n_members": args.n_members,
        "rng_seed": args.rng_seed,
        "r_definition": "R(t) = qkrig_variance(t) (re-kriged, direct)",
        "obs_dir": str(args.obs_dir),
        "params_file": str(params_file),
        "spinup": [SPINUP_START, SPINUP_END],
        "assimilation": [TEST_START, TEST_END],
        "kge": kge_val,
        "nse": nse_val,
        "n_updates": diag["n_updates"],
        "mean_abs_innovation_mm_h": diag["abs_innov_sum"] / n_up,
        "mean_pyy": diag["pyy_sum"] / n_up,
        "mean_K": {s: diag[f"K_{s}_sum"] / n_up for s in ("sm", "gw", "n0", "n1")},
        "mass_lost_at_gw_underflow_mm": diag["mass_lost_mm"],
        "runtime_s": round(time.time() - t_start, 1),
    }
    summary_json = out_dir / f"{cat_id}_enkf_summary.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    plot_png = out_dir / f"{cat_id}_enkf_plot.png"
    plot_results(df, cat_id, args.n_members, kge_val, nse_val, plot_png)

    print(f"  KGE={kge_val:.4f} NSE={nse_val:.4f} | {diag['n_updates']} updates | "
          f"mass lost {diag['mass_lost_mm']:.3f} mm | {summary['runtime_s']}s")
    print(f"  Saved: {results_csv.name}, {states_pq.name} "
          f"({len(member_state_rows)} state rows), {summary_json.name}, {plot_png.name}")
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="600-member EnKF over CFE for the 20% held-out experiment "
                    "(assimilation only, no forecast).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cat-id", required=True,
                        help="catchment id (cat-1016279) or 'all' for every CSV in --obs-dir")
    parser.add_argument("--n-members", type=int, default=600)
    parser.add_argument("--rng-seed", type=int, default=42)
    parser.add_argument("--obs-dir", default=DEFAULT_OBS_DIR)
    parser.add_argument("--params-dir", default=DEFAULT_PARAMS_DIR)
    parser.add_argument("--cfe-dir", default=DEFAULT_CFE_DIR)
    parser.add_argument("--config-file", default=DEFAULT_CONFIG_FILE)
    parser.add_argument("--test-forcing-dir1", default=DEFAULT_FORCING_DIR1)
    parser.add_argument("--test-forcing-dir2", default=DEFAULT_FORCING_DIR2)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--save-members", action="store_true",
                        help="archive every member's hourly forecast runoff to "
                             "<cat>_member_q.parquet (~25 MB/catchment at N=600)")
    parser.add_argument("--no-da", action="store_true",
                        help="open-loop mode: never assimilate (forcing noise and "
                             "process noise still applied). With --rng-seed 777042 "
                             "(= 42 + 777000) this reproduces run_cfe_enkf_forecast.py's "
                             "open-loop twin trajectory bit-exactly: same generator, "
                             "same draw order (build -> hourly step+noise), and the "
                             "skipped enkf_update never consumed draws in the twin either")
    args = parser.parse_args()

    global bmi_cfe
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    if args.cat_id == "all":
        cat_ids = sorted(p.stem for p in Path(args.obs_dir).glob("cat-*.csv"))
        print(f"Running all {len(cat_ids)} catchments from {args.obs_dir}")
    else:
        cat_ids = [args.cat_id]

    failed = []
    for cat_id in cat_ids:
        try:
            run_catchment(cat_id, args)
        except Exception:
            print(f"\nFAILED: {cat_id}")
            traceback.print_exc()
            failed.append(cat_id)

    if failed:
        print(f"\n{len(failed)}/{len(cat_ids)} catchments FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"\nAll {len(cat_ids)} catchment(s) completed.")


if __name__ == "__main__":
    main()
