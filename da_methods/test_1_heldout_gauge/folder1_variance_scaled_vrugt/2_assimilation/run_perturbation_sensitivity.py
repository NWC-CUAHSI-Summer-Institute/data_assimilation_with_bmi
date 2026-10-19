"""
Perturbation-source sensitivity analysis for the multi-catchment CFE DA experiment.

Ensemble-spread attribution: run 20 CFE members through the test
period with ONLY ONE perturbation source active at a time, no DA, no spinup.

Three sub-experiments (selected via --source):
    init     : perturb initial states once at t=0 (multiplicative ~5%)
    forcing  : perturb hourly precip (lognormal sigma=0.15) and PET (Gaussian sigma=0.10)
    process  : apply per-hour process noise on states (soil 0.2%, GW 0.15%, Nash 0.5%)

Each sub-experiment runs N=20 CFE members. DA is always OFF. Spinup is skipped — each
member starts from CFE's default initial state (modified by the perturbation if --source
init). This is the cleanest setup for attributing test-period spread to its source.

Output:
  <out-dir>/<cat-id>/<cat-id>_sensitivity_<source>.csv
  Columns: date, member_00, member_01, ..., member_19

Each member column is that member's hourly Q_sim (mm/h) over the test period.

Reads the same calibrated best_params.json as the production script. Does not need
the kriging obs file (no DA), but takes the same --obs-dir argument so it can also
log Q_obs alongside for plotting reference.

Usage:
  python3 run_perturbation_sensitivity.py \
    --cat-id cat-1016300 \
    --source forcing \
    --forcing-dir  <NWM retro dir> \
    --obs-dir      <kriging dir> \
    --cfe-dir      <cfe_py> \
    --config-file  <CFE BMI config JSON> \
    --param-bounds <CFE parameter bounds JSON> \
    --out-dir      <output dir> \
    --test-forcing-dir1 <NWM operational 2023-Feb2024> \
    --test-forcing-dir2 <NWM operational Feb2024-Sep2025>
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
N_MEMBERS = 20

# We only run the test window for this sensitivity analysis.
TEST_START = "2023-10-01 00:00:00"
TEST_END   = "2024-10-31 23:00:00"

# Perturbation magnitudes — same as production defaults
INIT_STATE_FRAC   = 0.05
PRECIP_PERTURB    = 0.15   # lognormal sigma
PET_PERTURB       = 0.10   # Gaussian sigma
SOIL_PROC_NOISE   = 0.002
GW_PROC_NOISE     = 0.0015
NASH_PROC_NOISE   = 0.005

# Set at runtime
CAT_ID            = None
FORCING_FILE      = None
TEST_FORCING_FILE = None
OBS_FILE          = None
CFE_CONFIG_FILE   = None
OUT_DIR           = None
bmi_cfe           = None


# ---------- Forcing helpers (same as production) ----------
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
    df['date'] = pd.to_datetime(df['time']).dt.strftime('%Y-%m-%d %H:%M:%S')
    df['total_precipitation'] = df['APCP_surface'] * 3600.0  # kg/m^2/s -> mm/h
    df['potential_evaporation'] = priestley_taylor_pet(
        df['DSWRF_surface'].values,
        df['TMP_2maboveground'].values,
    )
    return df


# ---------- Perturbation helpers ----------
def perturb_precip_lognormal(P, rng, sigma=PRECIP_PERTURB):
    """Lognormal multiplier with mean 1; preserves zero precip."""
    mu = -0.5 * sigma ** 2
    return float(P * rng.lognormal(mu, sigma))


def perturb_pet_gaussian(E, rng, sigma=PET_PERTURB):
    """Gaussian multiplicative noise on PET, clipped at 0."""
    return max(float(E * (1.0 + sigma * rng.standard_normal())), 0.0)


def apply_init_state_perturb(m, rng, frac=INIT_STATE_FRAC):
    """One-time multiplicative perturbation on the four states (called per member at t=0)."""
    sm_max = m.soil_reservoir["storage_max_m"]
    gw_max = m.gw_reservoir["storage_max_m"]
    sm0 = m.soil_reservoir["storage_m"]
    gw0 = m.gw_reservoir["storage_m"]
    n00 = float(m.nash_storage[0])
    n10 = float(m.nash_storage[1])
    m.soil_reservoir["storage_m"] = float(np.clip(
        sm0 * (1.0 + frac * rng.standard_normal()), 1e-6, sm_max))
    m.gw_reservoir["storage_m"] = float(np.clip(
        gw0 * (1.0 + frac * rng.standard_normal()), 1e-6, gw_max))
    # Nash buckets often start at 0; small additive jitter so they're not identical.
    jitter = max(sm_max * 1e-6, 1e-9)
    m.nash_storage[0] = max(n00 + jitter * rng.standard_normal(), 0.0)
    m.nash_storage[1] = max(n10 + jitter * rng.standard_normal(), 0.0)


def apply_process_noise(m, rng):
    """Per-timestep process noise on all 4 states for one member."""
    sm_max = m.soil_reservoir["storage_max_m"]
    gw_max = m.gw_reservoir["storage_max_m"]
    sm0 = m.soil_reservoir["storage_m"]
    gw0 = m.gw_reservoir["storage_m"]
    n00 = float(m.nash_storage[0])
    n10 = float(m.nash_storage[1])
    nash_floor = max(sm_max * 1e-4, 1e-7)
    m.soil_reservoir["storage_m"] = float(np.clip(
        sm0 * (1.0 + SOIL_PROC_NOISE * rng.standard_normal()), 0.0, sm_max))
    m.gw_reservoir["storage_m"] = float(np.clip(
        gw0 * (1.0 + GW_PROC_NOISE * rng.standard_normal()), 0.0, gw_max))
    m.nash_storage[0] = max(
        n00 * (1.0 + NASH_PROC_NOISE * rng.standard_normal())
        + nash_floor * rng.standard_normal(), 0.0)
    m.nash_storage[1] = max(
        n10 * (1.0 + NASH_PROC_NOISE * rng.standard_normal())
        + nash_floor * rng.standard_normal(), 0.0)


# ---------- Main per-source run ----------
def run_sensitivity(source, best_param_dict):
    """Run N members through the test period with only `source` perturbation active."""

    def custom_load_forcing(self_cfe):
        df = load_test_forcing()
        self_cfe.forcing_data = df.rename(columns={"date": "time"})

    # Build a per-catchment temporary CFE config with the calibrated params
    with open(CFE_CONFIG_FILE) as f:
        cfg = json.load(f)
    cfg['forcing_file']          = TEST_FORCING_FILE
    cfg['soil_params']['bb']     = best_param_dict['bb']
    cfg['soil_params']['smcmax'] = best_param_dict['smcmax']
    cfg['soil_params']['satdk']  = best_param_dict['satdk']
    cfg['slop']                  = best_param_dict['slop']
    cfg['max_gw_storage']        = best_param_dict['max_gw_storage']
    cfg['expon']                 = best_param_dict['expon']
    cfg['Cgw']                   = best_param_dict['Cgw']
    cfg['K_lf']                  = best_param_dict['K_lf']
    cfg['K_nash']                = best_param_dict['K_nash']
    cfg['partition_scheme']      = "Schaake" if best_param_dict['scheme'] <= 0.5 else "Xinanjiang"

    tmp_cfg = str(OUT_DIR / f'{CAT_ID}_bmi_config_temp_sensitivity_{source}.json')
    with open(tmp_cfg, 'w') as f:
        json.dump(cfg, f)

    # Deterministic seed per (catchment, source) so each run is reproducible
    seed = hash((CAT_ID, source)) & 0x7fffffff
    rng = np.random.default_rng(seed)

    print(f"[sensitivity] {CAT_ID} | source={source} | N={N_MEMBERS} | seed={seed}")

    # Build N members. Initial-state perturbation only happens when source == 'init'.
    models = []
    for i in range(N_MEMBERS):
        m = bmi_cfe.BMI_CFE(cfg_file=tmp_cfg)
        m.load_forcing_file = custom_load_forcing.__get__(m)
        m.initialize()
        if source == 'init' and i > 0:
            apply_init_state_perturb(m, rng)
        models.append(m)

    # Skip spinup. Go directly to the test period.
    df = load_test_forcing()
    t_mask = (df['date'] >= TEST_START) & (df['date'] <= TEST_END)
    df_test = df[t_mask]

    # Pre-allocate per-member Q arrays
    n_hours = len(df_test)
    q_matrix = np.full((n_hours, N_MEMBERS), np.nan, dtype=float)
    dates_out = df_test['date'].values

    for h, (p, e) in enumerate(zip(df_test['total_precipitation'],
                                    df_test['potential_evaporation'])):
        for i, m in enumerate(models):
            if source == 'forcing':
                p_i = perturb_precip_lognormal(p, rng)
                e_i = perturb_pet_gaussian(e, rng)
            else:
                p_i = float(p)
                e_i = float(e)
            m.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux', p_i / 1000)
            m.set_value('water_potential_evaporation_flux',                            e_i / 1000 / 3600)
            m.update()
            q_matrix[h, i] = m.get_value('land_surface_water__runoff_depth') * 1000  # m/h -> mm/h

        # Process noise applied after each hour, but only when source == 'process'
        if source == 'process':
            for m in models:
                apply_process_noise(m, rng)

    for m in models:
        m.finalize()

    # Write per-member CSV: date + 20 columns
    out = {'date': dates_out}
    for i in range(N_MEMBERS):
        out[f'member_{i:02d}'] = q_matrix[:, i]
    df_out = pd.DataFrame(out)
    out_path = OUT_DIR / f'{CAT_ID}_sensitivity_{source}.csv'
    df_out.to_csv(out_path, index=False)
    print(f"[sensitivity] saved {out_path}")
    print(f"[sensitivity] mean ensemble spread (std across members) over test period: "
          f"{float(q_matrix.std(axis=1).mean()):.5f} mm/h")


def main():
    global CAT_ID, FORCING_FILE, TEST_FORCING_FILE, OBS_FILE, CFE_CONFIG_FILE
    global OUT_DIR, bmi_cfe

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',             required=True)
    parser.add_argument('--source',             required=True,
                        choices=['init', 'forcing', 'process'],
                        help="Which perturbation source to isolate")
    parser.add_argument('--forcing-dir',        required=True)
    parser.add_argument('--obs-dir',            required=True,
                        help='Kept for symmetry with production script; not used since DA is off')
    parser.add_argument('--cfe-dir',            required=True)
    parser.add_argument('--config-file',        required=True)
    parser.add_argument('--param-bounds',       required=True,
                        help='Kept for symmetry; not read in this script')
    parser.add_argument('--out-dir',            required=True)
    parser.add_argument('--test-forcing-dir1',  default=None)
    parser.add_argument('--test-forcing-dir2',  default=None)
    args = parser.parse_args()

    CAT_ID            = args.cat_id
    FORCING_FILE      = os.path.join(args.forcing_dir, f'{CAT_ID}.csv')
    OBS_FILE          = os.path.join(args.obs_dir,     f'{CAT_ID}.csv')
    CFE_CONFIG_FILE   = args.config_file
    OUT_DIR           = Path(args.out_dir) / CAT_ID
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Build combined test forcing CSV (same idiom as production)
    if args.test_forcing_dir1 and args.test_forcing_dir2:
        f1 = os.path.join(args.test_forcing_dir1, f'{CAT_ID}.csv')
        f2 = os.path.join(args.test_forcing_dir2, f'{CAT_ID}.csv')
        if os.path.exists(f1) and os.path.exists(f2):
            df1 = pd.read_csv(f1)
            df2 = pd.read_csv(f2)
            combined = pd.concat([df1, df2], ignore_index=True)
            combined = combined.drop_duplicates(subset='time').sort_values('time')
            TEST_FORCING_FILE = str(OUT_DIR / f'{CAT_ID}_nwm_operational_combined.csv')
            combined.to_csv(TEST_FORCING_FILE, index=False)
        else:
            print(f"Warning: test forcing files not found for {CAT_ID}, exiting")
            return

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    best_params_file = OUT_DIR / f'{CAT_ID}_best_params.json'
    if not best_params_file.exists():
        print(f"No best params file at {best_params_file}. Pre-stage from Run 3.")
        return
    with open(best_params_file) as f:
        saved = json.load(f)
    best_param_dict = saved['best_parameters']

    run_sensitivity(args.source, best_param_dict)


if __name__ == '__main__':
    main()
