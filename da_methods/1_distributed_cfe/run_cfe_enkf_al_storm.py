#!/usr/bin/env python3
"""
run_cfe_enkf_al_storm.py -- 600-member EnKF for CFE, Alabama storm 2026.

Identical engine to run_cfe_enkf_hrrr.py.  Only three things changed:
  1. Paths  -- Alabama calibration, krig obs, HRRR forcing, output.
  2. Time windows -- Oct 1-3 2026 spinup; Oct 4-8 2026 assimilation.
  3. Storm save window -- Oct 7-8 2026 (restart points for the 18-h forecast).

HRRR forcing structure on server (same convention as run_cfe_enkf_hrrr.py):
  <hrrr-dir>/camels_YYYYMMDD/cat-<id>.csv   (24 rows/day)
  e.g. /mnt/disk2/al_storm_2026/hrrr/camels_20261001/cat-1011916.csv
  Columns: time, APCP_1hr_acc_fcst, DSWRF, TMP

USAGE
    # smoke test -- 4 members, one catchment
    python3 run_cfe_enkf_al_storm.py --cat-id cat-1011916 --n-members 4 \
        --out-dir /tmp/enkf_al_smoke

    # all 1072 catchments
    python3 run_cfe_enkf_al_storm.py --cat-id all

OUTPUTS (per catchment, under <out-dir>/<cat-id>/)
    <cat>_enkf_results.csv                  hourly: sim mean/std, obs, R, precip, da_updated
    <cat>_member_states_al_storm.parquet    member states at every storm-window hour (Oct 7-8)
    <cat>_member_q.parquet                  per-member q (only with --save-members)
    <cat>_enkf_summary.json                 settings, KGE/NSE, EnKF diagnostics
    <cat>_enkf_plot.png                     full assimilation period + storm zoom
"""

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bmi_cfe = None

# ---------------------------------------------------------------------------
# Alabama storm 2026 constants
# ---------------------------------------------------------------------------
DEFAULT_OBS_DIR     = "/mnt/disk2/al_storm_2026/catchment_obs_r100"
DEFAULT_PARAMS_DIR  = "/mnt/disk2/svyas_clean_run/calibration_al_storm_2026"
DEFAULT_CFE_DIR     = "/mnt/disk2/suma_helen_poster/cfe_py"
DEFAULT_CONFIG_FILE = "/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
DEFAULT_HRRR_DIR    = "/mnt/disk2/al_storm_2026/hrrr"
DEFAULT_OUT_DIR     = "/mnt/disk2/svyas_clean_run/da_al_storm_2026"

SPINUP_START  = "2026-10-01 00:00:00"   # 0-hour spinup (WY2026 spin-up handles init)
SPINUP_END    = "2026-10-01 00:00:00"
TEST_START    = "2026-10-01 01:00:00"   # DA starts from first NWM timestep Oct 1
TEST_END      = "2026-10-08 17:00:00"   # real-time data cutoff
STORM_START   = pd.Timestamp("2026-10-07 00:00:00")   # save states for forecast restart
STORM_END     = pd.Timestamp("2026-10-08 17:00:00")

# Perturbation magnitudes (same as run_cfe_enkf_hrrr.py -- HRRR is tighter than NWM operational)
PRECIP_SIGMA     = 0.15
PET_SIGMA        = 0.10
INIT_STATE_SIGMA = 0.05
SOIL_NOISE       = 0.002
GW_NOISE         = 0.0015
NASH_NOISE       = 0.005
R_FLOOR          = 1e-6

ALBEDO   = 0.20
ALPHA_PT = 1.26


# ---------------------------------------------------------------------------
# 1. LOAD
# ---------------------------------------------------------------------------
def priestley_taylor_pet(srad_wm2, T_kelvin):
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    pet = ALPHA_PT * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def build_hrrr_forcing(cat_id, hrrr_dir, cat_out_dir):
    """Concatenate all daily HRRR folders for one catchment into a single CSV.
    Cached: re-read from disk if already built."""
    out_csv = cat_out_dir / f"{cat_id}_hrrr_combined.csv"
    if out_csv.exists():
        print(f"  Reusing combined HRRR forcing: {out_csv}")
        return out_csv

    day_dirs = sorted(Path(hrrr_dir).glob("camels_????????"))
    if not day_dirs:
        raise FileNotFoundError(f"No camels_YYYYMMDD folders found in {hrrr_dir}")

    dfs = []
    for d in day_dirs:
        f = d / f"{cat_id}.csv"
        if f.exists():
            dfs.append(pd.read_csv(f))

    if not dfs:
        raise FileNotFoundError(f"No HRRR CSVs found for {cat_id} in {hrrr_dir}")

    combined = (pd.concat(dfs, ignore_index=True)
                .drop_duplicates(subset="time")
                .sort_values("time")
                .reset_index(drop=True))
    combined.to_csv(out_csv, index=False)
    print(f"  Wrote combined HRRR forcing: {out_csv.name} "
          f"({len(combined)} h, {combined['time'].iloc[0]} -> {combined['time'].iloc[-1]})")
    return out_csv


def build_nwm_forcing(cat_id, nwm_dir, cat_out_dir):
    """Read flat NWM operational CSV (nwm_dir/cat-{cat_id}.csv), deduplicate by
    valid time (keep first occurrence = shortest-lead, analysis-like), cache result."""
    out_csv = cat_out_dir / f"{cat_id}_nwm_combined.csv"
    if out_csv.exists():
        print(f"  Reusing combined NWM forcing: {out_csv}")
        return out_csv

    src = Path(nwm_dir) / f"{cat_id}.csv"
    if not src.exists():
        raise FileNotFoundError(f"No NWM forcing CSV for {cat_id} in {nwm_dir}")

    df = (pd.read_csv(src)
          .drop_duplicates(subset="time", keep="last")   # keep f001 (most recent cycle = analysis-quality)
          .sort_values("time")
          .reset_index(drop=True))
    df.to_csv(out_csv, index=False)
    print(f"  Wrote combined NWM forcing: {out_csv.name} "
          f"({len(df)} h, {df['time'].iloc[0]} -> {df['time'].iloc[-1]})")
    return out_csv


def load_forcing(path):
    """Read forcing CSV → DataFrame with date (str), total_precipitation (mm/h), potential_evaporation (mm/h).

    Auto-detects HRRR vs NWM operational columns:
      HRRR: APCP_1hr_acc_fcst (mm/h already), DSWRF (W/m²), TMP (K)
      NWM:  precip_rate (kg/m²/s → ×3600 = mm/h), DSWRF_surface (W/m²), TMP_2maboveground (K)
    """
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")

    if "precip_rate" in df.columns:
        # NWM operational — precip_rate is kg/m²/s, ×3600 → mm/h
        precip = (df["precip_rate"].fillna(0.0) * 3600.0).clip(lower=0.0)
        srad   = df["DSWRF_surface"].fillna(0.0).values
        temp_k = df["TMP_2maboveground"].ffill().values
    else:
        # HRRR — APCP_1hr_acc_fcst is already mm/h
        precip = df["APCP_1hr_acc_fcst"].fillna(0.0).clip(lower=0.0)
        srad   = df["DSWRF"].fillna(0.0).values
        temp_k = df["TMP"].ffill().values

    df["total_precipitation"]  = precip
    df["potential_evaporation"] = priestley_taylor_pet(srad, temp_k)
    return df


def load_observations(obs_csv):
    df = pd.read_csv(obs_csv)
    for col in ("time", "qkrig_mm_hr", "qkrig_variance"):
        if col not in df.columns:
            raise ValueError(f"{obs_csv} missing column '{col}' (found {list(df.columns)})")
    dates    = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    obs_dict = dict(zip(dates, df["qkrig_mm_hr"].astype(float)))
    var_dict = dict(zip(dates, df["qkrig_variance"].astype(float)))
    return obs_dict, var_dict


def write_bmi_config(template_file, best_params, forcing_csv, out_path):
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
# 2. BUILD
# ---------------------------------------------------------------------------
def build_ensemble(n_members, cfg_path, forcing_df, rng):
    bmi_forcing = forcing_df.rename(columns={"date": "time"})
    models = []
    for i in range(n_members):
        m = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
        m.load_forcing_file = (lambda self=m: setattr(self, "forcing_data", bmi_forcing))
        m.initialize()
        if i > 0:
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]
            m.soil_reservoir["storage_m"] = float(np.clip(
                m.soil_reservoir["storage_m"]
                * (1.0 + INIT_STATE_SIGMA * rng.standard_normal()), 1e-6, sm_max))
            m.gw_reservoir["storage_m"] = float(np.clip(
                m.gw_reservoir["storage_m"]
                * (1.0 + INIT_STATE_SIGMA * rng.standard_normal()), 1e-6, gw_max))
            jitter = max(sm_max * 1e-6, 1e-9)
            m.nash_storage[0] = max(float(m.nash_storage[0])
                                    + jitter * rng.standard_normal(), 0.0)
            m.nash_storage[1] = max(float(m.nash_storage[1])
                                    + jitter * rng.standard_normal(), 0.0)
        models.append(m)
    return models


# ---------------------------------------------------------------------------
# 3./4. STEP, NOISE, EnKF ANALYSIS
# ---------------------------------------------------------------------------
def step_all_members(models, precip_mmh, pet_mmh, rng):
    n = len(models)
    p_members = precip_mmh * rng.lognormal(-0.5 * PRECIP_SIGMA ** 2, PRECIP_SIGMA, n)
    e_members = np.maximum(pet_mmh * (1.0 + PET_SIGMA * rng.standard_normal(n)), 0.0)
    q_mmh = np.empty(n)
    for i, m in enumerate(models):
        m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                    float(p_members[i]) / 1000.0)
        m.set_value("water_potential_evaporation_flux",
                    float(e_members[i]) / 1000.0 / 3600.0)
        m.update()
        q_mmh[i] = m.get_value("land_surface_water__runoff_depth") * 1000.0
    return q_mmh


def add_process_noise(models, rng):
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
    n = len(models)
    q = np.asarray(q_forecast_mmh, dtype=float)
    if np.any(np.isnan(q)) or n < 2:
        return False
    R = max(float(krig_var), R_FLOOR)

    sm = np.array([m.get_value("SOIL_CONCEPTUAL_STORAGE") for m in models], dtype=float)
    gw = np.array([m.gw_reservoir["storage_m"] for m in models], dtype=float)
    n0 = np.array([float(m.nash_storage[0]) for m in models], dtype=float)
    n1 = np.array([float(m.nash_storage[1]) for m in models], dtype=float)

    q_anom = q - q.mean()
    Pyy = float(q_anom @ q_anom) / (n - 1)
    denom = Pyy + R
    if denom < 1e-12:
        return False

    K_sm = float((sm - sm.mean()) @ q_anom) / (n - 1) / denom
    K_gw = float((gw - gw.mean()) @ q_anom) / (n - 1) / denom
    K_n0 = float((n0 - n0.mean()) @ q_anom) / (n - 1) / denom
    K_n1 = float((n1 - n1.mean()) @ q_anom) / (n - 1) / denom

    obs_pert = y_obs + np.sqrt(R) * rng.standard_normal(n)
    innov = obs_pert - q

    for i, m in enumerate(models):
        sm_max = m.soil_reservoir["storage_max_m"]
        gw_max = m.gw_reservoir["storage_max_m"]
        sm_raw = sm[i] + K_sm * innov[i]
        gw_raw = gw[i] + K_gw * innov[i]
        n0_raw = n0[i] + K_n0 * innov[i]
        n1_raw = n1[i] + K_n1 * innov[i]

        so = max(sm_raw - sm_max, 0.0); sm_raw -= so; n0_raw += so
        go = max(gw_raw - gw_max, 0.0); gw_raw -= go; n1_raw += go
        n1u = max(-n1_raw, 0.0); n1_raw += n1u; n0_raw -= n1u
        n0u = max(-n0_raw, 0.0); n0_raw += n0u; sm_raw -= n0u
        smu = max(-sm_raw, 0.0); sm_raw += smu; gw_raw -= smu
        gwu = max(-gw_raw, 0.0); gw_raw += gwu
        diag["mass_lost_mm"] += gwu * 1000.0

        if np.isnan(sm_raw) or np.isnan(gw_raw) or np.isnan(n0_raw) or np.isnan(n1_raw):
            continue
        m.set_value("SOIL_CONCEPTUAL_STORAGE", float(np.clip(sm_raw, 0.0, sm_max)))
        m.gw_reservoir["storage_m"] = float(np.clip(gw_raw, 0.0, gw_max))
        m.nash_storage[0] = float(max(n0_raw, 0.0))
        m.nash_storage[1] = float(max(n1_raw, 0.0))

    diag["n_updates"]     += 1
    diag["abs_innov_sum"] += abs(float(innov.mean()))
    diag["pyy_sum"]       += Pyy
    diag["K_sm_sum"]      += K_sm
    diag["K_gw_sum"]      += K_gw
    diag["K_n0_sum"]      += K_n0
    diag["K_n1_sum"]      += K_n1

    if Pyy > 1e-4 or diag["n_updates"] % 500 == 0:
        print(f"  [EnKF N={n}] update {diag['n_updates']} | obs={y_obs:.3f} "
              f"q_mean={q.mean():.3f} innov={float(innov.mean()):+.4f} "
              f"R={R:.2e} Pyy={Pyy:.2e} K_sm={K_sm:+.2e}")
    return True


# ---------------------------------------------------------------------------
# 5. SAVE
# ---------------------------------------------------------------------------
def kge(obs, sim):
    ok = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[ok], sim[ok]
    if len(o) < 2: return float("nan")
    denom = float(np.sqrt(((o - o.mean()) ** 2).sum() * ((s - s.mean()) ** 2).sum()))
    if denom == 0 or o.std() == 0 or o.mean() == 0: return float("nan")
    r = float(((o - o.mean()) * (s - s.mean())).sum() / denom)
    return 1.0 - float(np.sqrt((r - 1) ** 2 + (s.std() / o.std() - 1) ** 2
                               + (s.mean() / o.mean() - 1) ** 2))


def nse(obs, sim):
    ok = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[ok], sim[ok]
    if len(o) < 2: return float("nan")
    denom = ((o - o.mean()) ** 2).sum()
    return float("nan") if denom == 0 else 1.0 - float(((o - s) ** 2).sum() / denom)


def plot_results(df, cat_id, n_members, kge_val, nse_val, out_png):
    dates = pd.to_datetime(df["date"])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10))
    for ax, lo, hi, title in (
        (ax1, dates.iloc[0], dates.iloc[-1],
         f"{cat_id} | HRRR forcing | EnKF N={n_members} | "
         f"KGE={kge_val:.3f} NSE={nse_val:.3f}"),
        (ax2, pd.Timestamp("2026-10-05"), pd.Timestamp("2026-10-08 23:00"),
         "Alabama storm zoom: Oct 5-8, 2026"),
    ):
        w = (dates >= lo) & (dates <= hi)
        ax.fill_between(dates[w],
                        (df["sim_mm_h"] - 2 * df["sim_std_mm_h"])[w].clip(lower=0),
                        (df["sim_mm_h"] + 2 * df["sim_std_mm_h"])[w],
                        color="tomato", alpha=0.25, label="ensemble +/- 2 std")
        ax.plot(dates[w], df["sim_mm_h"][w], "tomato", lw=1.5, label="ensemble mean")
        ax.plot(dates[w], df["obs_mm_h"][w], "k",      lw=1.0, label="qkrig obs")
        ax.set_ylabel("Discharge (mm/h)")
        ax.set_title(title)
        ax.legend(loc="upper left")
        twin = ax.twinx()
        twin.plot(dates[w], df["precip_mm_h"][w], "steelblue", lw=0.8, alpha=0.5)
        twin.set_ylim(twin.get_ylim()[::-1])
        twin.set_ylabel("Precip (mm/h)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Main per-catchment runner
# ---------------------------------------------------------------------------
def run_catchment(cat_id, args):
    t_start = time.time()
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'=' * 70}\n{cat_id}: HRRR forcing | EnKF N={args.n_members}, "
          f"seed={args.rng_seed}\n{'=' * 70}")

    params_file = Path(args.params_dir) / cat_id / f"{cat_id}_best_params.json"
    if not params_file.exists():
        params_file = Path(args.params_dir) / f"{cat_id}_best_params.json"
    with open(params_file) as f:
        best_params = json.load(f)["best_parameters"]

    # delete stale HRRR cache when caller extends the time window to new days
    if getattr(args, "rebuild_hrrr_cache", False):
        stale = out_dir / f"{cat_id}_hrrr_combined.csv"
        if stale.exists():
            stale.unlink()

    forcing_csv = build_hrrr_forcing(cat_id, args.hrrr_dir, out_dir)
    forcing = load_forcing(forcing_csv)
    obs_dict, var_dict = load_observations(Path(args.obs_dir) / f"{cat_id}.csv")

    cfg_path = write_bmi_config(args.config_file, best_params, forcing_csv,
                                out_dir / f"{cat_id}_bmi_config_enkf.json")

    spin = forcing[(forcing["date"] >= args.spinup_start)
                   & (forcing["date"] <= args.spinup_end)].reset_index(drop=True)
    test = forcing[(forcing["date"] >= args.test_start)
                   & (forcing["date"] <= args.test_end)].reset_index(drop=True)
    print(f"  Spinup {args.spinup_start[:10]} -> {args.spinup_end[:10]} ({len(spin)} h) | "
          f"assim {args.test_start[:10]} -> {args.test_end[:10]} ({len(test)} h)")
    if len(spin) == 0 or len(test) == 0:
        raise SystemExit(
            f"  EMPTY WINDOW for {cat_id}: forcing spans "
            f"{forcing['date'].iloc[0]} .. {forcing['date'].iloc[-1]}")

    rng = np.random.default_rng(args.rng_seed)
    print(f"  Building {args.n_members} ensemble members...")
    models = build_ensemble(args.n_members, str(cfg_path), forcing, rng)

    print("  Spinup (no DA)...")
    for k, row in enumerate(spin.itertuples()):
        step_all_members(models, float(row.total_precipitation),
                         float(row.potential_evaporation), rng)
        add_process_noise(models, rng)
        if (k + 1) % 24 == 0:
            print(f"    spinup {k+1}/{len(spin)} h ({time.time()-t_start:.0f}s)")

    print("  Assimilating...")
    diag = {"n_updates": 0, "abs_innov_sum": 0.0, "pyy_sum": 0.0, "K_sm_sum": 0.0,
            "K_gw_sum": 0.0, "K_n0_sum": 0.0, "K_n1_sum": 0.0, "mass_lost_mm": 0.0}
    records = []
    member_state_rows = []
    member_q_rows = [] if args.save_members else None

    for k, row in enumerate(test.itertuples()):
        date   = row.date
        precip = float(row.total_precipitation)
        pet    = float(row.potential_evaporation)

        q = step_all_members(models, precip, pet, rng)
        if member_q_rows is not None:
            member_q_rows.append(q.astype(np.float32))

        y_obs    = obs_dict.get(date, np.nan)
        krig_var = var_dict.get(date, np.nan)
        updated  = False
        if not (np.isnan(y_obs) or np.isnan(krig_var)):
            updated = enkf_update(models, q, y_obs, krig_var, rng, diag)
        add_process_noise(models, rng)

        ts = pd.Timestamp(date)
        if STORM_START <= ts <= STORM_END:
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
            "sim_mm_h":     float(q.mean()),
            "sim_std_mm_h": float(q.std(ddof=1)),
            "obs_mm_h":     float(y_obs) if not np.isnan(y_obs) else np.nan,
            "R_var":        float(krig_var) if not np.isnan(krig_var) else np.nan,
            "precip_mm_h":  precip,
            "da_updated":   int(updated),
        })
        if (k + 1) % 24 == 0:
            print(f"    assim {k+1}/{len(test)} h | updates: "
                  f"{diag['n_updates']} ({time.time()-t_start:.0f}s)")

    for m in models:
        m.finalize()

    df = pd.DataFrame(records)
    kge_val = kge(df["obs_mm_h"].values, df["sim_mm_h"].values)
    nse_val = nse(df["obs_mm_h"].values, df["sim_mm_h"].values)

    df.to_csv(out_dir / f"{cat_id}_enkf_results.csv", index=False)
    pd.DataFrame(member_state_rows).to_parquet(
        out_dir / f"{cat_id}_member_states_al_storm.parquet", index=False)

    if member_q_rows is not None:
        member_q = pd.DataFrame(
            np.vstack(member_q_rows),
            columns=[f"member_{i:04d}" for i in range(args.n_members)])
        member_q.insert(0, "time", df["date"].values)
        pq_path = out_dir / f"{cat_id}_member_q.parquet"
        member_q.to_parquet(pq_path, index=False)
        print(f"  Saved: {pq_path.name} ({len(member_q)} h x {args.n_members} members)")

    n_up = max(diag["n_updates"], 1)
    summary = {
        "cat_id": cat_id, "forcing": "HRRR",
        "n_members": args.n_members, "rng_seed": args.rng_seed,
        "spinup": [args.spinup_start, args.spinup_end],
        "assimilation": [args.test_start, args.test_end],
        "storm_states_saved": [str(STORM_START), str(STORM_END)],
        "kge": kge_val, "nse": nse_val,
        "n_updates": diag["n_updates"],
        "mean_abs_innovation_mm_h": diag["abs_innov_sum"] / n_up,
        "mean_pyy": diag["pyy_sum"] / n_up,
        "mean_K": {s: diag[f"K_{s}_sum"] / n_up for s in ("sm", "gw", "n0", "n1")},
        "mass_lost_at_gw_underflow_mm": diag["mass_lost_mm"],
        "runtime_s": round(time.time() - t_start, 1),
    }
    with open(out_dir / f"{cat_id}_enkf_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    plot_results(df, cat_id, args.n_members, kge_val, nse_val,
                 out_dir / f"{cat_id}_enkf_plot.png")
    print(f"  KGE={kge_val:.4f} NSE={nse_val:.4f} | {diag['n_updates']} updates | "
          f"mass lost {diag['mass_lost_mm']:.3f} mm | {summary['runtime_s']}s")
    return summary


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="600-member HRRR-forced EnKF for CFE — Alabama storm 2026.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cat-id", required=True,
                        help="catchment id (e.g. cat-1011916) or 'all'")
    parser.add_argument("--n-members",    type=int, default=600)
    parser.add_argument("--rng-seed",     type=int, default=42)
    parser.add_argument("--obs-dir",      default=DEFAULT_OBS_DIR)
    parser.add_argument("--params-dir",   default=DEFAULT_PARAMS_DIR)
    parser.add_argument("--cfe-dir",      default=DEFAULT_CFE_DIR)
    parser.add_argument("--config-file",  default=DEFAULT_CONFIG_FILE)
    parser.add_argument("--hrrr-dir",     default=DEFAULT_HRRR_DIR,
                        help="root dir containing camels_YYYYMMDD subfolders")
    parser.add_argument("--spinup-start", default=SPINUP_START)
    parser.add_argument("--spinup-end",   default=SPINUP_END)
    parser.add_argument("--test-start",   default=TEST_START)
    parser.add_argument("--test-end",     default=TEST_END)
    parser.add_argument("--out-dir",      default=DEFAULT_OUT_DIR)
    parser.add_argument("--save-members", action="store_true")
    parser.add_argument("--rebuild-hrrr-cache", action="store_true",
                        help="delete _hrrr_combined.csv cache before build (use when adding new HRRR days)")
    args = parser.parse_args()

    # make STORM_END follow args.test_end so member states are saved up to the new cutoff
    global bmi_cfe, STORM_START, STORM_END
    STORM_END   = pd.Timestamp(args.test_end)
    STORM_START = STORM_END - pd.Timedelta("48h")
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    if args.cat_id == "all":
        cat_ids = sorted(p.stem for p in Path(args.obs_dir).glob("cat-*.csv"))
        have_p  = {q.parent.name for q in
                   Path(args.params_dir).glob("cat-*/cat-*_best_params.json")}
        have_h  = {d.name for d in Path(args.hrrr_dir).glob("camels_????????")}
        if not have_p:
            print(f"WARNING: no best_params found in {args.params_dir}")
        else:
            skip = [c for c in cat_ids if c not in have_p]
            cat_ids = [c for c in cat_ids if c in have_p]
            if skip:
                print(f"SKIPPING {len(skip)} catchment(s) with no best_params: "
                      f"{', '.join(skip[:5])}{'...' if len(skip) > 5 else ''}")
        print(f"Running {len(cat_ids)} catchments | HRRR days: {len(have_h)}")
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
        print(f"\n{len(failed)}/{len(cat_ids)} FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"\nAll {len(cat_ids)} catchment(s) completed.")


if __name__ == "__main__":
    main()
