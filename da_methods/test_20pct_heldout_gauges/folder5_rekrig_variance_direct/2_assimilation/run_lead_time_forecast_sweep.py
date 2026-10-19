"""
Forecast lead-time evaluation.

For each issue-time t0 in a sampling schedule across the test period, fork the
ensemble at t0 and free-run 18 hours of forecast. Two parallel trajectories are
maintained throughout the test period:

  da:        Production setup. DA on, forcing perturbed, process noise on.
  openloop:  No DA from t=0. Forcing perturbed, process noise on.

At each issue time t0, the state of each scenario's ensemble is copied into a
forecast ensemble that free-runs 18 hours with:
  - DA off            (no observations consumed during the forecast window)
  - process noise off (no obs to collapse toward, so no anti-collapse needed)
  - forcing perturbed (lognormal precip, Gaussian PET — proxy for forecast
                       precip uncertainty in lieu of actual HEFS forecasts)

Issue-time schedule:
  - Base cadence: every --base-step-h hours across the full test period
                  (default 6h → ~1600 issue times for the year)
  - Densified to hourly across --dense-start..--dense-end if provided
                  (default: 2024-09-24 → 2024-09-28, the Helene window)

Outputs (stacked across all issue times):
  <out>/<cat>/<cat>_lead_time_forecasts_da.csv
  <out>/<cat>/<cat>_lead_time_forecasts_openloop.csv
  Columns: issue_time, lead_hour, valid_time, member_00 .. member_19
           (units: mm/h)

Post-processing pipeline (separate scripts, run after this):
  route_lead_time_forecasts.py — pushes each (t0, lead_hour) forecast through
                                  T-route to the gauge
  plot_lead_time_decay.py      — error-vs-lead-time curve, DA vs open-loop

Imports EnKFAssimilator from calibrate_catchment_cfe_da_v2.py so the DA math
matches production exactly.
"""

import argparse
import os
import sys
import json
import importlib.util
import numpy as np
import pandas as pd
from pathlib import Path

ALBEDO   = 0.20
ALPHA_PT = 1.26

SPINUP_START = "2023-02-01 00:00:00"
SPINUP_END   = "2023-09-30 23:00:00"
TEST_START   = "2023-10-01 00:00:00"
TEST_END     = "2024-10-31 23:00:00"

# Forecast lead time (hours after each issue time)
FORECAST_LEAD_HOURS = 18

# Default densification window — Hurricane Helene
DEFAULT_DENSE_START = "2024-09-24 00:00:00"
DEFAULT_DENSE_END   = "2024-09-28 23:00:00"


def priestley_taylor_pet(srad_wm2, T_kelvin, alpha=ALPHA_PT):
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    pet = alpha * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def load_test_forcing(test_forcing_file):
    df = pd.read_csv(test_forcing_file)
    df['date'] = pd.to_datetime(df['time']).dt.strftime('%Y-%m-%d %H:%M:%S')
    df['total_precipitation'] = df['APCP_surface'] * 3600.0
    df['potential_evaporation'] = priestley_taylor_pet(
        df['DSWRF_surface'].values, df['TMP_2maboveground'].values
    )
    return df


def import_enkf_class(script_path):
    """Load EnKFAssimilator from the production script without triggering its main()."""
    spec = importlib.util.spec_from_file_location("prod_v2", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.EnKFAssimilator


# ---------------- State snapshot / restore -----------------------------------
#
# The EnKF only ever touches these 4 states, so they're the only states that
# diverge between members. Copy these from src → dst to fork the ensemble.

def snapshot_states(models):
    return [
        {
            "soil_m":  float(m.soil_reservoir["storage_m"]),
            "gw_m":    float(m.gw_reservoir["storage_m"]),
            "nash0_m": float(m.nash_storage[0]),
            "nash1_m": float(m.nash_storage[1]),
        }
        for m in models
    ]


def restore_states(models, snapshot):
    for m, s in zip(models, snapshot):
        m.soil_reservoir["storage_m"] = s["soil_m"]
        m.gw_reservoir["storage_m"]   = s["gw_m"]
        m.nash_storage[0]             = s["nash0_m"]
        m.nash_storage[1]             = s["nash1_m"]


# ---------------- Issue-time schedule ----------------------------------------

def build_issue_time_schedule(dates_list, base_step_h, dense_start, dense_end):
    """Return a sorted list of issue-time strings sampled from dates_list.

    Every base_step_h hours across the full test period, plus every hour inside
    [dense_start, dense_end] if those are provided.
    """
    dates_dt = pd.to_datetime(dates_list)
    selected = set()

    # Base cadence (every base_step_h hours from the first date)
    base_mask = (np.arange(len(dates_dt)) % base_step_h == 0)
    for d in dates_dt[base_mask]:
        selected.add(d.strftime('%Y-%m-%d %H:%M:%S'))

    # Densified hourly window
    if dense_start is not None and dense_end is not None:
        ds = pd.Timestamp(dense_start)
        de = pd.Timestamp(dense_end)
        dense_mask = (dates_dt >= ds) & (dates_dt <= de)
        for d in dates_dt[dense_mask]:
            selected.add(d.strftime('%Y-%m-%d %H:%M:%S'))

    return sorted(selected)


# ---------------- CFE helpers ------------------------------------------------

def build_models(N, bmi_cfe, tmp_cfg, test_forcing_file, enkf,
                 apply_init_perturbation=True, member_zero_clean=True):
    """Build N CFE BMI instances with optional init-state perturbation."""
    def custom_load_forcing(self_cfe):
        df = load_test_forcing(test_forcing_file)
        self_cfe.forcing_data = df.rename(columns={"date": "time"})

    models = []
    for i in range(N):
        m = bmi_cfe.BMI_CFE(cfg_file=tmp_cfg)
        m.load_forcing_file = custom_load_forcing.__get__(m)
        m.initialize()
        if apply_init_perturbation and not (member_zero_clean and i == 0):
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]
            sm0 = m.soil_reservoir["storage_m"]
            gw0 = m.gw_reservoir["storage_m"]
            n00 = float(m.nash_storage[0])
            n10 = float(m.nash_storage[1])
            r = enkf.init_state_perturb_frac
            m.soil_reservoir["storage_m"] = float(np.clip(
                sm0 * (1 + r * enkf.rng.standard_normal()), 1e-6, sm_max))
            m.gw_reservoir["storage_m"]   = float(np.clip(
                gw0 * (1 + r * enkf.rng.standard_normal()), 1e-6, gw_max))
            jitter = max(sm_max * 1e-6, 1e-9)
            m.nash_storage[0] = max(n00 + jitter * enkf.rng.standard_normal(), 0.0)
            m.nash_storage[1] = max(n10 + jitter * enkf.rng.standard_normal(), 0.0)
        models.append(m)
    return models


def step_ensemble(models, p_arr, e_arr):
    """Advance each member one hour with member-specific (P, PET). Returns Q (mm/h)."""
    q = np.empty(len(models), dtype=float)
    for i, m in enumerate(models):
        m.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux',
                    float(p_arr[i]) / 1000.0)
        m.set_value('water_potential_evaporation_flux',
                    float(e_arr[i]) / 1000.0 / 3600.0)
        m.update()
        q[i] = m.get_value('land_surface_water__runoff_depth') * 1000.0  # m/h → mm/h
    return q


def do_sanity_check(prod_models, fcst_models, dates_list, forcing_by_date, h, N):
    """One-shot verification that the 4-state snapshot is sufficient to reproduce
    next-hour forecast Q.

    At main-loop hour h (after prod's hour-h step has completed including DA +
    process noise), snapshot prod, restore into fcst, then step BOTH ensembles
    one hour using identical unperturbed forcing. If snapshot/restore captures
    everything that matters, prod's and fcst's hour-(h+1) Q should match member-
    by-member to machine precision.

    Note: this mutates prod's state (h+1 step with non-perturbed forcing and no
    DA), so the caller should exit the script right after.
    """
    if h + 1 >= len(dates_list):
        print("[sanity-check] not enough horizon for h+1 step — skipping.")
        return False
    snap = snapshot_states(prod_models)
    restore_states(fcst_models, snap)

    next_date = dates_list[h + 1]
    p_next, e_next = forcing_by_date[next_date]
    p_uniform = np.full(N, p_next, dtype=float)
    e_uniform = np.full(N, e_next, dtype=float)

    q_prod = step_ensemble(prod_models, p_uniform, e_uniform)
    q_fcst = step_ensemble(fcst_models, p_uniform, e_uniform)

    diff = q_prod - q_fcst
    max_abs = float(np.max(np.abs(diff)))
    print("[sanity-check] forecast Q at h+1 with identical unperturbed forcing:")
    print(f"  prod first 5 members: {q_prod[:5]}")
    print(f"  fcst first 5 members: {q_fcst[:5]}")
    print(f"  max |prod - fcst|:    {max_abs:.3e} mm/h")
    if max_abs < 1e-9:
        print("  OK: snapshot/restore preserves state to machine precision.")
        return True
    print("  WARN: nonzero divergence — likely a state not in the 4-state snapshot.")
    return False


def run_forecast(fcst_models, enkf_fcst, dates_list, forcing_by_date,
                 t0_idx, n_lead):
    """Run an n_lead-hour free-forecast starting from the current state of fcst_models.

    Forcing is perturbed per member (lognormal precip + Gaussian PET, same as
    production). No DA. No process noise.

    Returns: (lead_hours, valid_times, q_matrix shape (n_lead, N)).
    """
    N = len(fcst_models)
    q_matrix = np.full((n_lead, N), np.nan, dtype=float)
    lead_hours = []
    valid_times = []
    for lead in range(1, n_lead + 1):
        idx = t0_idx + lead
        if idx >= len(dates_list):
            break
        valid_date = dates_list[idx]
        p, e = forcing_by_date[valid_date]
        p_arr, e_arr = enkf_fcst.perturb_forcing(p, e)
        q = step_ensemble(fcst_models, p_arr, e_arr)
        q_matrix[lead - 1, :] = q
        lead_hours.append(lead)
        valid_times.append(valid_date)
    return lead_hours, valid_times, q_matrix


# ---------------- Main run ---------------------------------------------------

def run(args, bmi_cfe, EnKFAssimilator):
    cat_id = args.cat_id
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Combined test forcing (same idiom as the per-member script)
    f1 = os.path.join(args.test_forcing_dir1, f'{cat_id}.csv')
    f2 = os.path.join(args.test_forcing_dir2, f'{cat_id}.csv')
    df1 = pd.read_csv(f1)
    df2 = pd.read_csv(f2)
    combined = pd.concat([df1, df2], ignore_index=True)
    combined = combined.drop_duplicates(subset='time').sort_values('time')
    test_forcing_file = str(out_dir / f'{cat_id}_nwm_operational_combined.csv')
    combined.to_csv(test_forcing_file, index=False)

    # Pre-staged calibrated params
    best_params_file = out_dir / f'{cat_id}_best_params.json'
    if not best_params_file.exists():
        print(f"Missing {best_params_file}. Pre-stage from Run 3.")
        return
    with open(best_params_file) as f:
        best = json.load(f)['best_parameters']

    # Build temp CFE config with calibrated params
    with open(args.config_file) as f:
        cfg = json.load(f)
    cfg['forcing_file']          = test_forcing_file
    cfg['soil_params']['bb']     = best['bb']
    cfg['soil_params']['smcmax'] = best['smcmax']
    cfg['soil_params']['satdk']  = best['satdk']
    cfg['slop']                  = best['slop']
    cfg['max_gw_storage']        = best['max_gw_storage']
    cfg['expon']                 = best['expon']
    cfg['Cgw']                   = best['Cgw']
    cfg['K_lf']                  = best['K_lf']
    cfg['K_nash']                = best['K_nash']
    cfg['partition_scheme']      = "Schaake" if best['scheme'] <= 0.5 else "Xinanjiang"
    tmp_cfg = str(out_dir / f'{cat_id}_bmi_config_temp_leadtime.json')
    with open(tmp_cfg, 'w') as f:
        json.dump(cfg, f)

    # Three EnKF instances:
    #   enkf_da     — runs the DA-on trajectory; consumes obs each hour
    #   enkf_ol     — runs the open-loop trajectory; never calls update_states
    #   enkf_fcst   — used only inside the 18-hour free-runs (perturb_forcing only)
    # All three share defaults with production. Distinct RNG seeds keep them
    # independent so the openloop and forecast spreads are not coupled to DA.
    obs_file = os.path.join(args.obs_dir, f'{cat_id}.csv')
    seed_base = args.rng_seed if args.rng_seed is not None else 0
    enkf_da = EnKFAssimilator(
        n_members=args.enkf_members, obs_error_std=args.enkf_obs_error_std,
        obs_file=obs_file, use_vrugt_r=(not args.no_vrugt_r),
        vrugt_alpha=args.vrugt_alpha, vrugt_scale=args.vrugt_scale,
        rng_seed=seed_base if args.rng_seed is not None else None,
    )
    enkf_ol = EnKFAssimilator(
        n_members=args.enkf_members, obs_error_std=args.enkf_obs_error_std,
        obs_file=obs_file, use_vrugt_r=(not args.no_vrugt_r),
        vrugt_alpha=args.vrugt_alpha, vrugt_scale=args.vrugt_scale,
        rng_seed=(seed_base + 1) if args.rng_seed is not None else None,
    )
    enkf_fcst = EnKFAssimilator(
        n_members=args.enkf_members, obs_error_std=args.enkf_obs_error_std,
        obs_file=obs_file, use_vrugt_r=(not args.no_vrugt_r),
        vrugt_alpha=args.vrugt_alpha, vrugt_scale=args.vrugt_scale,
        rng_seed=(seed_base + 2) if args.rng_seed is not None else None,
    )
    # Optional R override — replace every per-hour obs variance with a constant.
    # Matches DualEarth/new_EnKF.py (R=0.07). Set via --hardcoded-r on the CLI.
    # Effect: at storm peak ≈ existing Vrugt R, at low flow ≈ 15-17× larger →
    # near-zero gain at low flow where kriging obs is noisy.
    if args.hardcoded_r is not None:
        for enkf in (enkf_da, enkf_ol, enkf_fcst):
            for date in enkf.obs_var_dict:
                enkf.obs_var_dict[date] = args.hardcoded_r
        print(f"[lead-time] R hardcoded to {args.hardcoded_r} mm^2/h^2 "
              f"(overrides Vrugt + kriging variance formula)")

    N = enkf_da.n_members
    print(f"[lead-time] {cat_id} | N={N} | DA on/off both run | "
          f"forecast lead = {FORECAST_LEAD_HOURS}h")

    # Build three ensembles (init-perturbed except the forecast ensemble, whose
    # state will be overwritten at each issue time before each free-run)
    prod_models   = build_models(N, bmi_cfe, tmp_cfg, test_forcing_file, enkf_da)
    openloop_models = build_models(N, bmi_cfe, tmp_cfg, test_forcing_file, enkf_ol)
    fcst_models   = build_models(N, bmi_cfe, tmp_cfg, test_forcing_file,
                                  enkf_fcst, apply_init_perturbation=False)

    # ----- Spinup with perturbed forcing on both trajectories -----
    df = load_test_forcing(test_forcing_file)
    sp_mask = (df['date'] >= SPINUP_START) & (df['date'] <= SPINUP_END)
    df_sp = df[sp_mask]
    print(f"[lead-time] spinup: {len(df_sp)} hours")
    for p, e in zip(df_sp['total_precipitation'], df_sp['potential_evaporation']):
        p_da, e_da = enkf_da.perturb_forcing(p, e)
        step_ensemble(prod_models, p_da, e_da)
        enkf_da.add_process_noise(prod_models)
        p_ol, e_ol = enkf_ol.perturb_forcing(p, e)
        step_ensemble(openloop_models, p_ol, e_ol)
        enkf_ol.add_process_noise(openloop_models)

    # ----- Test period: step both, fork forecasts at issue times -----
    t_mask = (df['date'] >= TEST_START) & (df['date'] <= TEST_END)
    df_test = df[t_mask].reset_index(drop=True)
    dates_list = list(df_test['date'].values)
    forcing_by_date = dict(zip(
        df_test['date'].values,
        zip(df_test['total_precipitation'].values,
            df_test['potential_evaporation'].values)))

    issue_times = build_issue_time_schedule(
        dates_list, args.base_step_h,
        args.dense_start, args.dense_end)
    issue_set = set(issue_times)
    print(f"[lead-time] test period: {len(dates_list)} hours | "
          f"{len(issue_times)} issue times "
          f"(base={args.base_step_h}h, dense={args.dense_start or 'none'}"
          f"..{args.dense_end or 'none'})")

    # Pre-build the date→index map so forecast windows are fast to look up
    date_to_idx = {d: i for i, d in enumerate(dates_list)}

    da_rows = []  # (issue_time, lead_hour, valid_time, *member_values)
    ol_rows = []
    sanity_done = False

    for h, current_date in enumerate(dates_list):
        p, e = forcing_by_date[current_date]

        # --- DA trajectory step ---
        p_da, e_da = enkf_da.perturb_forcing(p, e)
        q_da = step_ensemble(prod_models, p_da, e_da)
        enkf_da.update_states(prod_models, current_date, q_da)
        enkf_da.add_process_noise(prod_models)

        # --- Open-loop trajectory step (no DA) ---
        p_ol, e_ol = enkf_ol.perturb_forcing(p, e)
        q_ol = step_ensemble(openloop_models, p_ol, e_ol)
        enkf_ol.add_process_noise(openloop_models)

        # --- One-shot sanity check at the first scheduled issue time ---
        if args.sanity_check and not sanity_done and current_date in issue_set:
            ok = do_sanity_check(
                prod_models, fcst_models, dates_list, forcing_by_date, h, N)
            sanity_done = True
            print(f"[sanity-check] exiting (re-run without --sanity-check for "
                  f"the full sweep). Result: {'OK' if ok else 'FAIL'}")
            for m in prod_models + openloop_models + fcst_models:
                m.finalize()
            sys.exit(0 if ok else 1)

        # --- Fork forecasts at scheduled issue times ---
        if current_date in issue_set:
            for scenario_label, src_models, rows_acc in [
                ('da', prod_models, da_rows),
                ('openloop', openloop_models, ol_rows),
            ]:
                snap = snapshot_states(src_models)
                restore_states(fcst_models, snap)
                leads, valids, qm = run_forecast(
                    fcst_models, enkf_fcst,
                    dates_list, forcing_by_date,
                    t0_idx=h, n_lead=FORECAST_LEAD_HOURS)
                for k, (lead, vt) in enumerate(zip(leads, valids)):
                    rows_acc.append((current_date, lead, vt, *qm[k, :]))

    for m in prod_models + openloop_models + fcst_models:
        m.finalize()

    # ----- Save forecast CSVs -----
    cols = ['issue_time', 'lead_hour', 'valid_time'] + [f'member_{i:02d}' for i in range(N)]
    da_path = out_dir / f'{cat_id}_lead_time_forecasts_da.csv'
    ol_path = out_dir / f'{cat_id}_lead_time_forecasts_openloop.csv'
    pd.DataFrame(da_rows, columns=cols).to_csv(da_path, index=False)
    pd.DataFrame(ol_rows, columns=cols).to_csv(ol_path, index=False)
    print(f"[lead-time] saved {da_path}")
    print(f"[lead-time] saved {ol_path}")

    # Issue-time schedule (small, useful for downstream scripts)
    sched_path = out_dir / f'{cat_id}_lead_time_issue_times.csv'
    pd.DataFrame({'issue_time': issue_times}).to_csv(sched_path, index=False)
    print(f"[lead-time] saved {sched_path}")

    print(f"[lead-time] DA assimilator: updates={enkf_da.n_updates} | "
          f"avg Pyy={enkf_da.avg_pyy:.6f} | "
          f"mass lost={enkf_da.total_overflow_lost_mm:.3f} mm")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',             required=True)
    parser.add_argument('--forcing-dir',        required=True)
    parser.add_argument('--obs-dir',            required=True)
    parser.add_argument('--cfe-dir',            required=True)
    parser.add_argument('--config-file',        required=True)
    parser.add_argument('--param-bounds',       required=True,
                        help='Kept for symmetry; not read in this script')
    parser.add_argument('--out-dir',            required=True)
    parser.add_argument('--test-forcing-dir1',  required=True)
    parser.add_argument('--test-forcing-dir2',  required=True)
    parser.add_argument('--enkf-members',       type=int, default=20)
    parser.add_argument('--enkf-obs-error-std', type=float, default=0.05)
    parser.add_argument('--no-vrugt-r',         action='store_true')
    parser.add_argument('--vrugt-alpha',        type=float, default=0.10)
    parser.add_argument('--vrugt-scale',        type=float, default=0.001)
    parser.add_argument('--rng-seed',           type=int, default=None)
    parser.add_argument('--base-step-h',        type=int, default=6,
                        help='Issue-time cadence in hours across the full test '
                             'period (default 6h)')
    parser.add_argument('--dense-start',        type=str, default=DEFAULT_DENSE_START,
                        help='Start of hourly-densification window '
                             '(default 2024-09-24 00:00:00 — Helene)')
    parser.add_argument('--dense-end',          type=str, default=DEFAULT_DENSE_END,
                        help='End of hourly-densification window '
                             '(default 2024-09-28 23:00:00 — Helene)')
    parser.add_argument('--sanity-check',       action='store_true',
                        help='At the first issue time, fork fcst_models and step '
                             'both prod and fcst one hour with identical unperturbed '
                             'forcing. Confirms the 4-state snapshot is sufficient '
                             'to reproduce next-hour Q to machine precision, then '
                             'exits. Use before kicking off the full sweep.')
    parser.add_argument('--hardcoded-r',        type=float, default=None,
                        help='If set, override the Vrugt + kriging-variance R formula '
                             'and use this constant value (mm^2/h^2) at every hour. '
                             'Matches the DualEarth/new_EnKF.py reference (R=0.07). '
                             'Suppresses DA at low flow where kriging is noisy while '
                             'preserving DA strength at storm peaks. Applied to all '
                             'three EnKF instances (da, openloop, fcst) for consistency.')
    parser.add_argument('--prod-script',        default=None,
                        help='Path to calibrate_catchment_cfe_da_v2.py for '
                             'importing EnKFAssimilator. Defaults to next-to-this-file.')
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    prod_script = args.prod_script or os.path.join(here, 'calibrate_catchment_cfe_da_v2.py')
    if not os.path.exists(prod_script):
        raise FileNotFoundError(
            f"Could not find production script at {prod_script}. "
            f"Pass --prod-script to override.")
    EnKFAssimilator = import_enkf_class(prod_script)

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe

    run(args, _bmi_cfe, EnKFAssimilator)


if __name__ == '__main__':
    main()
