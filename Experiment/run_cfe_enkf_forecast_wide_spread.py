#!/usr/bin/env python3
"""
run_cfe_enkf_forecast_wide_spread.py -- hourly-cycled 18-hour forecasts from the
600-member WIDE-SPREAD EnKF, Helene window. Companion to run_cfe_enkf_wide_spread.py.

WIDE-SPREAD VARIANT of run_cfe_enkf_forecast.py: identical forecast logic, but it
imports run_cfe_enkf_wide_spread (PRECIP_SIGMA=0.50, PET_SIGMA=0.20) instead of
run_cfe_enkf, and every member carries its PERSISTENT precip bias (PRECIP_BIAS_SIGMA
=0.40, drawn once) through spinup, assimilation AND the 18-hour forecast -- so the
forecasts are consistent with the wide-spread analysis (and the routed gauge run).
The DA ensemble and the open-loop twin each draw their own bias from their own RNG.

WHAT THIS SCRIPT DOES, IN ONE PARAGRAPH
    It reproduces the exact same 600-member assimilation as run_cfe_enkf_wide_spread.py
    (same seed, same RNG stream -> bit-identical analysis states), but stops
    treating Sep 20 - Sep 29, 2024 as plain assimilation: at every hour in that
    window it (1) assimilates the new qkrig observation, (2) saves all member
    states, (3) freezes the ensemble, runs an 18-hour open-loop forecast for
    every member, (4) rewinds to the frozen state and moves to the next hour.
    The result is 240 hourly forecast issues x 18 leads x N members per
    catchment -- the operational "issue a forecast every hour" picture.

THE TIMELINE
    SPINUP        2023-02-01 .. 2023-09-30      no DA (spread develops)
    ASSIMILATE    2023-10-01 .. 2024-09-19 23h  hourly EnKF (as run_cfe_enkf)
    ISSUE SWEEP   2024-09-20 00h .. 2024-09-29 23h, for every hour t0:
                      step + EnKF update at t0 + process noise   (analysis)
                      save analysis states of every member at t0
                      snapshot full model state (4 stores + GIUH queue)
                      forecast t0+1 .. t0+18, no DA               (open loop)
                      restore snapshot, continue to t0+1

FORECAST PHYSICS CHOICES (flags below)
    - forcing during forecast = NWM operational (perfect-forcing hindcast) by
      default, OR -- with --forecast-forcing-dir -- the HRRR forecast forcing
      (a real forecast). Either way it is perturbed per member exactly like the
      assimilation step: the member's PERSISTENT precip bias x fresh hourly
      lognormal noise (PERTURB_FORECAST_FORCING). Assimilation ALWAYS stays NWM.
    - no process noise during forecast (ADD_FORECAST_PROCESS_NOISE = False):
      spread comes from analysis-state diversity + forcing noise only.
    - the forecast RNG is a separate, per-issue-time generator, so forecasting
      consumes NOTHING from the assimilation RNG stream: the analysis
      trajectory is bit-identical to run_cfe_enkf_wide_spread.py with the same
      seed (the Helene member states here must match the wide-spread run's
      <cat>_member_states_helene.parquet -- a built-in correctness check).

WHY SNAPSHOT/RESTORE IS SAFE
    bmi_cfe.update() takes its forcing from set_value() calls (no internal
    forcing index), so the complete evolving state of a member is exactly:
        soil_reservoir["storage_m"], gw_reservoir["storage_m"],
        nash_storage[0..1], runoff_queue_m_per_timestep (GIUH queue).
    The vol* accumulators keep counting forecast water but are pure
    diagnostics; this script does not report mass balance.

USAGE
    # one catchment, full 600 members (~6-8 min)
    ~/miniforge3/envs/t-route/bin/python3 run_cfe_enkf_forecast_wide_spread.py --cat-id cat-1016279

    # all 21 catchments, one after the other
    ~/miniforge3/envs/t-route/bin/python3 run_cfe_enkf_forecast_wide_spread.py --cat-id all

    # quick smoke test
    ~/miniforge3/envs/t-route/bin/python3 run_cfe_enkf_forecast_wide_spread.py \
        --cat-id cat-1016279 --n-members 4 --out-dir /tmp/enkf_fc_smoke

OPEN-LOOP TWIN
    A second, never-assimilated 600-member ensemble (same build, same physics,
    own RNG stream) runs in parallel through spinup + the whole test period and
    issues its own 18-hour forecast at every issue time. The DA-vs-open-loop
    gap at each lead measures exactly what assimilated initial conditions buy.
    PAIRED DESIGN: at each issue time both forecast sets reuse the SAME
    per-issue forcing-perturbation draws, so the only difference between the
    DA forecast and the open-loop forecast is the initial state.

OUTPUTS (per catchment, under <out-dir>/<cat-id>/)
    <cat>_forecast_leadtime_da.parquet        one row per (issue_time, lead_hour):
                                      issue_time, lead_hour, valid_time,
                                      member_0000..NNNN  (mm/h, float32)
                                      -- 240 x 18 rows at full sweep
    <cat>_forecast_leadtime_openloop.parquet  same shape, open-loop twin
    <cat>_analysis_states_sweep.parquet  DA member states at every issue hour
                                      (time, member, soil_m, gw_m, nash0, nash1)
    <cat>_forecast_summary.json       settings + per-lead KGE/NSE of the
                                      ensemble-mean forecast vs qkrig obs,
                                      for both da and openloop
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
matplotlib.use("Agg")  # headless -- must be set before bmi_cfe imports pyplot

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cfe_enkf_wide_spread as enkf  # reuse: loaders, ensemble, EnKF step, metrics

# ---------------------------------------------------------------------------
# Sweep constants (CLI-overridable)
# ---------------------------------------------------------------------------
ISSUE_START = "2024-09-20 00:00:00"   # first forecast issue hour
ISSUE_END   = "2024-09-29 23:00:00"   # last forecast issue hour (10 days)
LEAD_HOURS  = 18                      # forecast horizon per issue

PERTURB_FORECAST_FORCING  = True      # per-member met noise during forecast
ADD_FORECAST_PROCESS_NOISE = False    # no state noise during forecast

DEFAULT_OUT_DIR = "/mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_wide_spread"
DEFAULT_DA_TRAJ_DIR = "/mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_da"
DEFAULT_OL_TRAJ_DIR = "/mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_ol"


# ---------------------------------------------------------------------------
# HRRR forecast forcing (optional, --forecast-forcing-dir) -- the 18h leads only
# ---------------------------------------------------------------------------
def load_hrrr_forecast_forcing(hrrr_dir, cat_id, lo, hi):
    """Concatenate the daily HRRR files camels_<YYYYMMDD>/<cat>.csv spanning
    [lo, hi] into one valid-time-indexed forcing series.

    HRRR differs from the NWM loader (enkf.load_forcing) in two ways:
      - precip column is APCP_1hr_acc_fcst and is ALREADY mm/h (1-hour
        accumulation) -- NO *3600 (NWM's APCP_surface is kg/m2/s and needs it);
      - met columns are DSWRF / TMP (vs NWM's DSWRF_surface / TMP_2maboveground).
    PET is the same Priestley-Taylor formula as NWM (reused from enkf)."""
    days = pd.date_range(lo.normalize(), hi.normalize(), freq="D")
    frames = []
    for d in days:
        f = Path(hrrr_dir) / f"camels_{d:%Y%m%d}" / f"{cat_id}.csv"
        if f.exists():
            frames.append(pd.read_csv(f))
    if not frames:
        raise FileNotFoundError(
            f"No HRRR forcing for {cat_id} in {hrrr_dir} over {lo:%Y-%m-%d}..{hi:%Y-%m-%d}")
    df = pd.concat(frames, ignore_index=True)
    df["date"] = pd.to_datetime(df["time"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    df = df.drop_duplicates(subset="date").sort_values("date").reset_index(drop=True)
    df["total_precipitation"] = df["APCP_1hr_acc_fcst"]             # already mm/h
    df["potential_evaporation"] = enkf.priestley_taylor_pet(
        df["DSWRF"].values, df["TMP"].values)
    return df


# ---------------------------------------------------------------------------
# Snapshot / restore -- the full evolving state of a member (see header)
# ---------------------------------------------------------------------------
def snapshot_states(models):
    return [(float(m.soil_reservoir["storage_m"]),
             float(m.gw_reservoir["storage_m"]),
             np.array(m.nash_storage, copy=True),
             np.array(m.runoff_queue_m_per_timestep, copy=True))
            for m in models]


def restore_states(models, snap):
    for m, (sm, gw, nash, queue) in zip(models, snap):
        m.soil_reservoir["storage_m"] = sm
        m.gw_reservoir["storage_m"] = gw
        m.nash_storage = nash.copy()
        m.runoff_queue_m_per_timestep = queue.copy()


def forecast_step(models, precip_mmh, pet_mmh, precip_bias, rng):
    """One open-loop forecast hour: per-member perturbed (or deterministic)
    forcing, NO DA, NO process noise. Returns per-member runoff (mm/h).
    WIDE-SPREAD: perturbed precip carries each member's PERSISTENT bias
    (precip_bias) x fresh hourly lognormal noise -- the same forcing model the
    member saw during assimilation, so the forecast stays consistent."""
    n = len(models)
    if PERTURB_FORECAST_FORCING:
        p = precip_mmh * precip_bias * rng.lognormal(-0.5 * enkf.PRECIP_SIGMA ** 2,
                                                     enkf.PRECIP_SIGMA, n)
        e = np.maximum(pet_mmh * (1.0 + enkf.PET_SIGMA * rng.standard_normal(n)), 0.0)
    else:
        p = np.full(n, precip_mmh)
        e = np.full(n, max(pet_mmh, 0.0))
    q_mmh = np.empty(n)
    for i, m in enumerate(models):
        m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                    float(p[i]) / 1000.0)             # mm/h -> m/h
        m.set_value("water_potential_evaporation_flux",
                    float(e[i]) / 1000.0 / 3600.0)    # mm/h -> m/s
        m.update()
        q_mmh[i] = m.get_value("land_surface_water__runoff_depth") * 1000.0
    if ADD_FORECAST_PROCESS_NOISE:
        enkf.add_process_noise(models, rng)
    return q_mmh


# ---------------------------------------------------------------------------
# The whole story for one catchment
# ---------------------------------------------------------------------------
def run_catchment(cat_id, args):
    t_start = time.time()
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'=' * 70}\n{cat_id}: EnKF N={args.n_members} + hourly {args.lead_hours}h "
          f"forecasts {args.issue_start[:10]} .. {args.issue_end[:10]}\n{'=' * 70}")

    # ---- 1. LOAD (identical to run_cfe_enkf_wide_spread) ------------------
    params_file = Path(args.params_dir) / cat_id / f"{cat_id}_best_params.json"
    with open(params_file) as f:
        best_params = json.load(f)["best_parameters"]
    forcing_csv = enkf.build_combined_forcing(
        cat_id, args.test_forcing_dir1, args.test_forcing_dir2, out_dir)
    forcing = enkf.load_forcing(forcing_csv)
    obs_dict, var_dict = enkf.load_observations(Path(args.obs_dir) / f"{cat_id}.csv")
    cfg_path = enkf.write_bmi_config(args.config_file, best_params, forcing_csv,
                                     out_dir / f"{cat_id}_bmi_config_enkf.json")

    spin = forcing[(forcing["date"] >= enkf.SPINUP_START)
                   & (forcing["date"] <= enkf.SPINUP_END)].reset_index(drop=True)
    # assimilation-only stretch: TEST_START up to the hour BEFORE the sweep
    assim = forcing[(forcing["date"] >= enkf.TEST_START)
                    & (forcing["date"] < args.issue_start)].reset_index(drop=True)
    sweep = forcing[(forcing["date"] >= args.issue_start)
                    & (forcing["date"] <= args.issue_end)].reset_index(drop=True)
    # forecast leads run past ISSUE_END; index the whole forcing table by date
    forcing_by_date = forcing.set_index("date")
    # Forecast forcing: HRRR if --forecast-forcing-dir is given (assimilation
    # always stays on the NWM table above); else reuse NWM (perfect-forcing).
    if args.forecast_forcing_dir:
        fc_lo = pd.Timestamp(args.issue_start)
        fc_hi = pd.Timestamp(args.issue_end) + pd.Timedelta(hours=args.lead_hours)
        fc_forcing_by_date = load_hrrr_forecast_forcing(
            args.forecast_forcing_dir, cat_id, fc_lo, fc_hi).set_index("date")
        print(f"  Forecast forcing: HRRR ({len(fc_forcing_by_date)} h) from {args.forecast_forcing_dir}")
    else:
        fc_forcing_by_date = forcing_by_date
        print("  Forecast forcing: NWM (perfect-forcing hindcast)")
    print(f"  Spinup {len(spin)} h | assimilation {len(assim)} h | "
          f"issue sweep {len(sweep)} h x {args.lead_hours} leads")

    # ---- 2. BUILD + 3. SPINUP + plain ASSIMILATION ------------------------
    # Same seed, same draw order as run_cfe_enkf_wide_spread.py -> identical trajectories.
    rng = np.random.default_rng(args.rng_seed)
    print(f"  Building {args.n_members} ensemble members...")
    models = enkf.build_ensemble(args.n_members, str(cfg_path), forcing, rng)
    # WIDE-SPREAD: persistent per-member precip bias, drawn ONCE right after the
    # build (matching run_cfe_enkf_wide_spread's RNG order) and held all run.
    precip_bias = rng.lognormal(-0.5 * enkf.PRECIP_BIAS_SIGMA ** 2,
                                enkf.PRECIP_BIAS_SIGMA, args.n_members)

    # Open-loop twin: identical construction, NEVER assimilates. It gets its
    # own RNG stream so the DA ensemble's draw sequence stays bit-identical to
    # run_cfe_enkf_wide_spread.py with the same seed.
    ol_rng = np.random.default_rng(args.rng_seed + 777_000)
    print(f"  Building {args.n_members} open-loop twin members...")
    ol_models = enkf.build_ensemble(args.n_members, str(cfg_path), forcing, ol_rng)
    ol_precip_bias = ol_rng.lognormal(-0.5 * enkf.PRECIP_BIAS_SIGMA ** 2,
                                      enkf.PRECIP_BIAS_SIGMA, args.n_members)

    print("  Spinup (both ensembles)...")
    for k, row in enumerate(spin.itertuples()):
        p, e = float(row.total_precipitation), float(row.potential_evaporation)
        enkf.step_all_members(models, p, e, precip_bias, rng)
        enkf.add_process_noise(models, rng)
        enkf.step_all_members(ol_models, p, e, ol_precip_bias, ol_rng)
        enkf.add_process_noise(ol_models, ol_rng)
        if (k + 1) % 2000 == 0:
            print(f"    spinup {k + 1}/{len(spin)} h ({time.time() - t_start:.0f}s)")

    diag = {"n_updates": 0, "abs_innov_sum": 0.0, "pyy_sum": 0.0, "K_sm_sum": 0.0,
            "K_gw_sum": 0.0, "K_n0_sum": 0.0, "K_n1_sum": 0.0, "mass_lost_mm": 0.0}

    def assimilate_one_hour(row):
        """One DA cycle, exactly as in run_cfe_enkf: step -> update -> noise."""
        q = enkf.step_all_members(models, float(row.total_precipitation),
                                  float(row.potential_evaporation), precip_bias, rng)
        y_obs = obs_dict.get(row.date, np.nan)
        krig_var = var_dict.get(row.date, np.nan)
        if not (np.isnan(y_obs) or np.isnan(krig_var)):
            enkf.enkf_update(models, q, y_obs, krig_var, rng, diag)
        enkf.add_process_noise(models, rng)
        return q

    def step_openloop_one_hour(row):
        """Open-loop twin: step + process noise, never any DA update. Returns the
        per-member runoff (mm/h) -- the open-loop member_q."""
        q = enkf.step_all_members(ol_models, float(row.total_precipitation),
                                  float(row.potential_evaporation), ol_precip_bias, ol_rng)
        enkf.add_process_noise(ol_models, ol_rng)
        return q

    # member_q buffers (--save-members): per-hour DA and open-loop runoff over
    # assimilation + sweep. Same draw order as run_cfe_enkf_wide_spread.py, so this is
    # bit-identical to that script's member_q (DA) and its --no-da twin (open loop) --
    # the trajectory warm-up routing gets the same antecedent inflow either way.
    save_mq = args.save_members
    da_mq_rows = [] if save_mq else None
    ol_mq_rows = [] if save_mq else None
    mq_times   = [] if save_mq else None

    print("  Assimilating up to the issue window (open loop runs alongside)...")
    for k, row in enumerate(assim.itertuples()):
        q_da = assimilate_one_hour(row)
        q_ol = step_openloop_one_hour(row)
        if save_mq:
            mq_times.append(row.date)
            da_mq_rows.append(q_da.astype(np.float32))
            ol_mq_rows.append(q_ol.astype(np.float32))
        if (k + 1) % 2000 == 0:
            print(f"    assim {k + 1}/{len(assim)} h | updates: {diag['n_updates']} "
                  f"({time.time() - t_start:.0f}s)")

    # ---- 4. ISSUE SWEEP ----------------------------------------------------
    print("  Issue sweep (assimilate -> save states -> forecast both -> rewind)...")
    member_names = [f"member_{i:04d}" for i in range(args.n_members)]
    fc_meta, fc_rows_da, fc_rows_ol = [], [], []   # forecast output buffers
    state_rows = []                    # analysis states at every issue hour

    def forecast_18h(model_set, t0_ts, fc_bias, fc_rng, sink):
        """Freeze model_set, run the lead window, rewind. Appends to sink.
        fc_bias = the ensemble's persistent per-member precip bias."""
        snap = snapshot_states(model_set)
        for lead in range(1, args.lead_hours + 1):
            tv = (t0_ts + pd.Timedelta(hours=lead)).strftime("%Y-%m-%d %H:%M:%S")
            frow = fc_forcing_by_date.loc[tv]
            qf = forecast_step(model_set, float(frow["total_precipitation"]),
                               float(frow["potential_evaporation"]), fc_bias, fc_rng)
            sink.append(qf.astype(np.float32))
        restore_states(model_set, snap)

    for k, row in enumerate(sweep.itertuples()):
        t0 = row.date
        q_da = assimilate_one_hour(row)       # analysis state at t0
        q_ol = step_openloop_one_hour(row)    # twin advances without DA
        if save_mq:
            mq_times.append(t0)
            da_mq_rows.append(q_da.astype(np.float32))
            ol_mq_rows.append(q_ol.astype(np.float32))

        # save every member's analysis state (the forecast initial conditions)
        for i, m in enumerate(models):
            state_rows.append({
                "time": t0, "member": i,
                "soil_m": float(m.soil_reservoir["storage_m"]),
                "gw_m":   float(m.gw_reservoir["storage_m"]),
                "nash0":  float(m.nash_storage[0]),
                "nash1":  float(m.nash_storage[1]),
            })

        # Per-issue forecast RNG, independent of the assimilation stream.
        # PAIRED DESIGN: the open-loop forecast re-seeds the SAME generator,
        # so both forecasts see identical forcing perturbations -- the only
        # difference between them is the initial state (DA vs never-DA).
        t0_ts = pd.Timestamp(t0)
        for lead in range(1, args.lead_hours + 1):
            fc_meta.append((t0, lead,
                            (t0_ts + pd.Timedelta(hours=lead)).strftime("%Y-%m-%d %H:%M:%S")))
        forecast_18h(models, t0_ts, precip_bias,
                     np.random.default_rng(args.rng_seed * 1_000_003 + k), fc_rows_da)
        forecast_18h(ol_models, t0_ts, ol_precip_bias,
                     np.random.default_rng(args.rng_seed * 1_000_003 + k), fc_rows_ol)

        if (k + 1) % 24 == 0:
            print(f"    issued {k + 1}/{len(sweep)} forecast pairs "
                  f"({time.time() - t_start:.0f}s)")

    for m in models:
        m.finalize()
    for m in ol_models:
        m.finalize()

    # ---- 5. SAVE -----------------------------------------------------------
    meta = pd.DataFrame(fc_meta, columns=["issue_time", "lead_hour", "valid_time"])

    def save_forecast(rows, tag):
        fc = pd.concat([meta.copy(),
                        pd.DataFrame(np.vstack(rows), columns=member_names)], axis=1)
        pq = out_dir / f"{cat_id}_forecast_leadtime_{tag}.parquet"
        fc.to_parquet(pq, index=False)
        # per-lead skill of the ensemble-mean forecast vs qkrig obs (metadata
        # only -- the parquet keeps every member)
        ens = fc[member_names].mean(axis=1)
        obs = fc["valid_time"].map(obs_dict)
        skill = {}
        for lead, g in fc.groupby("lead_hour").groups.items():
            o = obs.loc[g].values.astype(float)
            s = ens.loc[g].values
            skill[int(lead)] = {"kge": enkf.kge(o, s), "nse": enkf.nse(o, s),
                                "n": int(np.isfinite(o).sum())}
        return pq, skill

    fc_pq, skill_by_lead = save_forecast(fc_rows_da, "da")
    ol_pq, skill_by_lead_ol = save_forecast(fc_rows_ol, "openloop")

    states_pq = out_dir / f"{cat_id}_analysis_states_sweep.parquet"
    pd.DataFrame(state_rows).to_parquet(states_pq, index=False)

    if save_mq:
        # DA member_q -> --da-traj-dir, open-loop member_q -> --ol-traj-dir, each as
        # <traj-dir>/<cat>/<cat>_member_q.parquet (exactly what route_forecast.py reads).
        def save_member_q(rows, traj_dir):
            d = Path(traj_dir) / cat_id
            d.mkdir(parents=True, exist_ok=True)
            mq = pd.DataFrame(np.vstack(rows), columns=member_names)
            mq.insert(0, "time", mq_times)
            pq = d / f"{cat_id}_member_q.parquet"
            mq.to_parquet(pq, index=False)
            return pq
        save_member_q(da_mq_rows, args.da_traj_dir)
        save_member_q(ol_mq_rows, args.ol_traj_dir)
        print(f"  Saved member_q ({len(mq_times)} h): "
              f"{args.da_traj_dir}/{cat_id}/, {args.ol_traj_dir}/{cat_id}/")

    summary = {
        "cat_id": cat_id,
        "n_members": args.n_members,
        "rng_seed": args.rng_seed,
        "issue_window": [args.issue_start, args.issue_end],
        "n_issues": len(sweep),
        "lead_hours": args.lead_hours,
        "perturb_forecast_forcing": PERTURB_FORECAST_FORCING,
        "forecast_process_noise": ADD_FORECAST_PROCESS_NOISE,
        "r_definition": "R(t) = qkrig_variance(t) (re-kriged, direct)",
        "forecast_forcing": ("HRRR:" + str(args.forecast_forcing_dir)
                             if args.forecast_forcing_dir else "NWM (perfect-forcing)"),
        "obs_dir": str(args.obs_dir),
        "params_file": str(params_file),
        "n_assim_updates": diag["n_updates"],
        "skill_by_lead_ens_mean_da": skill_by_lead,
        "skill_by_lead_ens_mean_openloop": skill_by_lead_ol,
        "runtime_s": round(time.time() - t_start, 1),
    }
    with open(out_dir / f"{cat_id}_forecast_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    l1, lN = skill_by_lead.get(1, {}), skill_by_lead.get(args.lead_hours, {})
    o1, oN = skill_by_lead_ol.get(1, {}), skill_by_lead_ol.get(args.lead_hours, {})
    print(f"  Done: {len(sweep)} issues x {args.lead_hours} leads | "
          f"DA  lead-1 KGE={l1.get('kge', float('nan')):.3f} "
          f"lead-{args.lead_hours}={lN.get('kge', float('nan')):.3f} | "
          f"OL  lead-1={o1.get('kge', float('nan')):.3f} "
          f"lead-{args.lead_hours}={oN.get('kge', float('nan')):.3f} | "
          f"{summary['runtime_s']}s")
    print(f"  Saved: {fc_pq.name}, {ol_pq.name}, {states_pq.name} "
          f"({len(state_rows)} state rows)")
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Hourly-cycled 18h forecasts from the 600-member EnKF "
                    "(Sep 20-29, 2024).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cat-id", required=True,
                        help="catchment id (cat-1016279) or 'all'")
    parser.add_argument("--n-members", type=int, default=600)
    parser.add_argument("--rng-seed", type=int, default=42)
    parser.add_argument("--lead-hours", type=int, default=LEAD_HOURS)
    parser.add_argument("--issue-start", default=ISSUE_START)
    parser.add_argument("--issue-end", default=ISSUE_END)
    parser.add_argument("--obs-dir", default=enkf.DEFAULT_OBS_DIR)
    parser.add_argument("--params-dir", default=enkf.DEFAULT_PARAMS_DIR)
    parser.add_argument("--cfe-dir", default=enkf.DEFAULT_CFE_DIR)
    parser.add_argument("--config-file", default=enkf.DEFAULT_CONFIG_FILE)
    parser.add_argument("--test-forcing-dir1", default=enkf.DEFAULT_FORCING_DIR1)
    parser.add_argument("--test-forcing-dir2", default=enkf.DEFAULT_FORCING_DIR2)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--forecast-forcing-dir", default=None,
                        help="HRRR forecast-forcing dir (camels_<YYYYMMDD>/<cat>.csv). "
                             "If set, the 18h forecast leads use HRRR (precip already "
                             "mm/h); assimilation stays on NWM. If unset, forecasts use "
                             "the NWM forcing (perfect-forcing hindcast).")
    parser.add_argument("--save-members", action="store_true",
                        help="also dump per-hour DA + open-loop member_q (runoff) "
                             "trajectories for the trajectory warm-up routing -- removes "
                             "the need for a separate run_cfe_enkf_wide_spread.py run")
    parser.add_argument("--da-traj-dir", default=DEFAULT_DA_TRAJ_DIR,
                        help="where to write the DA member_q (with --save-members)")
    parser.add_argument("--ol-traj-dir", default=DEFAULT_OL_TRAJ_DIR,
                        help="where to write the open-loop member_q (with --save-members)")
    args = parser.parse_args()

    # bmi_cfe is imported here and handed to the run_cfe_enkf module, whose
    # build_ensemble/step functions reference it as a module global.
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    enkf.bmi_cfe = _bmi_cfe

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
        print(f"\nFailed catchments: {failed}")
        sys.exit(1)
    print(f"\nAll {len(cat_ids)} catchment(s) completed.")


if __name__ == "__main__":
    main()
