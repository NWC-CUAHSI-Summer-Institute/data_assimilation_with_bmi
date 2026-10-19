"""
Run the production DA pipeline for ONE catchment and save per-member streamflow.

Identical to calibrate_catchment_cfe_da_v2.py in algorithm (init perturbation +
forcing perturbation + process noise + true EnKF with Vrugt R), but writes one
column per ensemble member instead of only the ensemble mean. This is for
the per-member factor-decomposition plot.

To keep the comparison apples-to-apples with the sensitivity runs (which skip
spinup), this script also SKIPS the spinup loop. Each member starts directly
from CFE's default initial state, plus the t=0 init perturbation.

Output:
  <out-dir>/<cat-id>/<cat-id>_production_per_member.csv
  Columns: date, member_00, member_01, ..., member_19

Imports EnKFAssimilator from calibrate_catchment_cfe_da_v2 to guarantee the
DA math matches production exactly.

Usage:
  python3 run_production_per_member.py \
    --cat-id cat-1016300 \
    --forcing-dir /mnt/.../nwm_retro_catchment_forcings \
    --obs-dir /mnt/.../catchment_ts_03463300_with_variance \
    --cfe-dir /mnt/.../cfe_py \
    --config-file /mnt/.../cat_03463300_bmi_config_cfe.json \
    --param-bounds /mnt/.../CFE_parameter_bounds.json \
    --out-dir /mnt/.../v2_production_per_member \
    --test-forcing-dir1 ... \
    --test-forcing-dir2 ... \
    --enkf-members 20
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

TEST_START = "2023-10-01 00:00:00"
TEST_END   = "2024-10-31 23:00:00"


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


def run(args, bmi_cfe, EnKFAssimilator):
    cat_id = args.cat_id
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Combined test forcing
    f1 = os.path.join(args.test_forcing_dir1, f'{cat_id}.csv')
    f2 = os.path.join(args.test_forcing_dir2, f'{cat_id}.csv')
    df1 = pd.read_csv(f1)
    df2 = pd.read_csv(f2)
    combined = pd.concat([df1, df2], ignore_index=True)
    combined = combined.drop_duplicates(subset='time').sort_values('time')
    test_forcing_file = str(out_dir / f'{cat_id}_nwm_operational_combined.csv')
    combined.to_csv(test_forcing_file, index=False)

    # Load best_params
    best_params_file = out_dir / f'{cat_id}_best_params.json'
    if not best_params_file.exists():
        print(f"Missing {best_params_file}. Pre-stage from Run 3.")
        return
    with open(best_params_file) as f:
        saved = json.load(f)
    best = saved['best_parameters']

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
    tmp_cfg = str(out_dir / f'{cat_id}_bmi_config_temp_per_member.json')
    with open(tmp_cfg, 'w') as f:
        json.dump(cfg, f)

    # Build the EnKF assimilator (same defaults as production)
    obs_file = os.path.join(args.obs_dir, f'{cat_id}.csv')
    enkf = EnKFAssimilator(
        n_members=args.enkf_members,
        obs_error_std=args.enkf_obs_error_std,
        obs_file=obs_file,
        use_vrugt_r=(not args.no_vrugt_r),
        vrugt_alpha=args.vrugt_alpha,
        vrugt_scale=args.vrugt_scale,
        rng_seed=args.rng_seed,
    )
    N = enkf.n_members
    print(f"[per-member] {cat_id} | N={N} | use_vrugt_r={enkf.use_vrugt_r}")
    if getattr(args, 'hardcoded_r', None) is not None:
        enkf.obs_var_dict = {k: args.hardcoded_r for k in enkf.obs_var_dict}
        print(f"[per-member] obs_var_dict overridden: R = {args.hardcoded_r} (fixed)")

    # Custom forcing loader for CFE BMI
    def custom_load_forcing(self_cfe):
        df = load_test_forcing(test_forcing_file)
        self_cfe.forcing_data = df.rename(columns={"date": "time"})

    # Build N CFE instances; perturb initial states (member 0 left clean)
    models = []
    init_states_per_member = []  # capture initial-state snapshot after perturbation
    for i in range(N):
        m = bmi_cfe.BMI_CFE(cfg_file=tmp_cfg)
        m.load_forcing_file = custom_load_forcing.__get__(m)
        m.initialize()
        if i > 0:
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
        # Snapshot the (possibly-perturbed) initial state for this member
        init_states_per_member.append({
            "soil_m":   float(m.soil_reservoir["storage_m"]),
            "gw_m":     float(m.gw_reservoir["storage_m"]),
            "nash0_m":  float(m.nash_storage[0]),
            "nash1_m":  float(m.nash_storage[1]),
            "soil_max_m": float(m.soil_reservoir["storage_max_m"]),
            "gw_max_m":   float(m.gw_reservoir["storage_max_m"]),
        })
        models.append(m)

    # SKIP SPINUP — go directly to the test period (matches sensitivity-script setup)
    df = load_test_forcing(test_forcing_file)
    t_mask = (df['date'] >= TEST_START) & (df['date'] <= TEST_END)
    df_test = df[t_mask]
    n_hours = len(df_test)
    q_matrix = np.full((n_hours, N), np.nan, dtype=float)
    precip_matrix = np.full((n_hours, N), np.nan, dtype=float)  # mm/h per member
    pet_matrix    = np.full((n_hours, N), np.nan, dtype=float)  # mm/h per member
    dates_out = df_test['date'].values

    for h, (p, e, current_date) in enumerate(zip(
            df_test['total_precipitation'],
            df_test['potential_evaporation'],
            df_test['date'])):

        # Perturb forcing per member
        p_arr, e_arr = enkf.perturb_forcing(p, e)
        # Capture the per-member perturbed forcing for later plotting/diagnostics
        precip_matrix[h, :] = p_arr
        pet_matrix[h, :]    = e_arr

        # Advance each member one hour
        ensemble_q = np.empty(N, dtype=float)
        for i, m in enumerate(models):
            m.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux',
                        float(p_arr[i]) / 1000)
            m.set_value('water_potential_evaporation_flux',
                        float(e_arr[i]) / 1000 / 3600)
            m.update()
            ensemble_q[i] = m.get_value('land_surface_water__runoff_depth') * 1000  # mm/h

        # Record pre-analysis forecast per member (this is the "actual production" value)
        q_matrix[h, :] = ensemble_q

        # DA step (uses production EnKFAssimilator)
        enkf.update_states(models, current_date, ensemble_q)

        # Process noise after DA
        enkf.add_process_noise(models)

    for m in models:
        m.finalize()

    # ----- Save per-member CSVs -----
    def _save_matrix(matrix, suffix):
        out = {'date': dates_out}
        for i in range(N):
            out[f'member_{i:02d}'] = matrix[:, i]
        df_out = pd.DataFrame(out)
        out_path = out_dir / f'{cat_id}_production_per_member_{suffix}.csv'
        df_out.to_csv(out_path, index=False)
        return out_path

    # Streamflow outputs
    q_path = out_dir / f'{cat_id}_production_per_member.csv'
    out = {'date': dates_out}
    for i in range(N):
        out[f'member_{i:02d}'] = q_matrix[:, i]
    pd.DataFrame(out).to_csv(q_path, index=False)
    print(f"[per-member] saved {q_path}")

    # Perturbed forcings (so each member's actual inputs are recoverable)
    precip_path = _save_matrix(precip_matrix, "precip")
    pet_path    = _save_matrix(pet_matrix,    "pet")
    print(f"[per-member] saved {precip_path}")
    print(f"[per-member] saved {pet_path}")

    # Initial states (small, JSON is fine)
    init_states_path = out_dir / f'{cat_id}_production_per_member_initial_states.json'
    with open(init_states_path, 'w') as f:
        json.dump({
            "catchment_id": cat_id,
            "n_members": N,
            "members": {f"member_{i:02d}": init_states_per_member[i] for i in range(N)},
        }, f, indent=2)
    print(f"[per-member] saved {init_states_path}")

    print(f"[per-member] ensemble mean spread over test period: "
          f"{float(q_matrix.std(axis=1).mean()):.5f} mm/h")


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
    parser.add_argument('--hardcoded-r',        type=float, default=None,
                        help='Fix R to this constant value for all timesteps (overrides Vrugt formula)')
    parser.add_argument('--vrugt-alpha',        type=float, default=0.10)
    parser.add_argument('--vrugt-scale',        type=float, default=0.001)
    parser.add_argument('--rng-seed',           type=int, default=None)
    parser.add_argument('--prod-script',        default=None,
                        help='Path to calibrate_catchment_cfe_da_v2.py for importing '
                             'EnKFAssimilator. Defaults to the script located next to this file.')
    args = parser.parse_args()

    # Locate the production script (default: same dir as this script)
    if args.prod_script is None:
        here = os.path.dirname(os.path.abspath(__file__))
        prod_script = os.path.join(here, 'calibrate_catchment_cfe_da_v2.py')
    else:
        prod_script = args.prod_script
    if not os.path.exists(prod_script):
        raise FileNotFoundError(f"Could not find production script at {prod_script}. "
                                f"Pass --prod-script to override.")
    EnKFAssimilator = import_enkf_class(prod_script)

    # Dynamic CFE import (same idiom as production script)
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe

    run(args, _bmi_cfe, EnKFAssimilator)


if __name__ == '__main__':
    main()
