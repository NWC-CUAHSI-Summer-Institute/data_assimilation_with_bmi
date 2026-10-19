"""
Per-catchment CFE test run with optional EnKF data assimilation.

This script runs the test period (Oct 2023 – Oct 2024, includes Hurricane
Helene) on top of pre-calibrated CFE parameters, with optional state
assimilation against per-catchment kriging streamflow observations.

  The CFE calibration step (DDS against the 2020–2022 retro period) is NOT
  part of this script. It is performed once, separately, via the calibrate-cfe
  workflow at NWC-CUAHSI-Summer-Institute/calibrate-cfe (held-out gauges branch).
  This script loads the resulting <cat-id>_best_params.json from --out-dir and
  uses those parameters fixed through the entire test run.

================================================================================
ENSEMBLE KALMAN FILTER (true stochastic EnKF — Burgers/van Leeuwen/Evensen 1998)
================================================================================
Each hour during the test period:
  - N CFE members run in parallel with perturbed precip + PET
  - forecast variance Pyy = ensemble variance of Q (no heuristic)
  - Kalman gain VECTOR: one K per state, from cross-covariance Pxy(state, Q)
    (no hardcoded state-split weights — the data decides each hour)
  - observation perturbed N times: obs_i = obs + sqrt(R) × N(0,1)
  - each member updated with its own obs_i and its own innovation
  - per-state per-timestep process noise prevents ensemble collapse
  - output time series is the ensemble mean

Measured improvement on USGS gauge 03463300 (21 sub-catchments, Helene year):
  mean Test KGE 0.692 (no DA) → 0.817 (with DA).

================================================================================
STATES UPDATED (4 states, applied per ensemble member with mass cascade)
================================================================================
  1. soil_reservoir["storage_m"]      via BMI: SOIL_CONCEPTUAL_STORAGE
  2. gw_reservoir["storage_m"]        direct attribute (not exposed via BMI)
  3. nash_storage[0]                  direct attribute (upstream Nash bucket)
  4. nash_storage[1]                  direct attribute (feeds streamflow)

Mass-conserving cascade — when a state hits its bound, the excess/deficit is
redirected along CFE's physical flow direction instead of being silently clipped:
    Positive corrections (water added):
       soil full → excess pushed to Nash[0]
       GW   full → excess pushed to Nash[1]
    Negative corrections (water removed):
       Nash[1] < 0 → deficit absorbed from Nash[0]
       Nash[0] < 0 → deficit absorbed from soil
       soil    < 0 → deficit absorbed from GW
       GW      < 0 → accept loss (logged in mm via total_overflow_lost_mm)

================================================================================
INPUTS
================================================================================
Forcing (training):  per-catchment NWM retro CSV
                     columns: time, APCP_surface, DSWRF_surface, TMP_2maboveground, ...
Forcing (testing):   per-catchment NWM operational CSVs (two dirs concatenated)
                     APCP_surface in kg/m²/s → converted to mm/h inside
Obs (kriging):       per-catchment CSV. Three column layouts auto-handled:
                       a) datetime, qkrig                (format A)
                       b) datetime, qkrig, variance       (format A with-variance)
                       c) time, qkrig_mm_hr               (format B standard)
                       d) time, qkrig_mm_hr, qkrig_variance  (format B with-variance)
                     If a per-hour variance column is present it is used as R;
                     otherwise R = --enkf-obs-error-std^2 is used as fallback.

Time splits:
  Spinup (cal):   2019-01-01 00:00 → 2019-12-31 23:00
  Calibration:   2020-01-01 00:00 → 2022-12-31 23:00
  Spinup (test): 2023-02-01 00:00 → 2023-09-30 23:00
  Test (Helene): 2023-10-01 00:00 → 2024-10-31 23:00

================================================================================
KEY CLI FLAGS
================================================================================
  --enkf-enabled           Turn DA on
  --enkf-members 20        N ensemble members (>=2 required)
  --enkf-obs-error-std 0.05  Fallback obs std (mm/h) when no variance column
                              and --no-vrugt-r is set

OBSERVATION ERROR VARIANCE R
  Production default: Vrugt et al. 2005 (SODA paper) heteroscedastic R
  scaled by the kriging variance:
      R(t) = (alpha * y_obs(t))^2 + scale * sigma^2_krig(t)
  Flow-magnitude term keeps R small at low flow (so DA fires) and larger
  at peaks (where obs is also more uncertain). Kriging variance brings
  per-hour catchment-specific info in, scaled so it does not dominate.
  Defaults: alpha=0.10 (10% relative error per Vrugt), scale=0.001.

  Pass --no-vrugt-r to revert to the raw kriging variance behavior (or
  --enkf-obs-error-std^2 fallback when no variance column). Not
  recommended for this basin — raw kriging variance is 8-10x the obs
  magnitude, which collapses the Kalman gain.

Test-loop perturbation defaults (see EnKFAssimilator.__init__):
    precip:         lognormal multiplier with sigma=0.15, mean=1
                    (literature standard; guarantees non-negative)
    PET:            ×N(1, 0.10), clipped at 0
    initial state:  ×N(1, 0.05) on soil/GW; small additive jitter on Nash
    observation:    additive noise with std sqrt(R) (Burgers/Evensen)
    process noise per state (per-hour multiplicative; prevents ensemble
    collapse — required by textbook EnKF):
      soil   : 0.002  (0.2%)    slow-evolving, small noise sufficient
      GW     : 0.0015 (0.15%)   slowest state
      Nash   : 0.005  (0.5%)    + additive floor; needs more because nash
                                spread does not develop naturally
    Magnitudes match NWC-CUAHSI/data_assimilation_with_bmi reference impl.

================================================================================
USAGE
================================================================================
Production run (true EnKF, N=20, Vrugt R on by default):
  python3 calibrate_catchment_cfe_da_v2.py \\
    --cat-id cat-1016300 \\
    --forcing-dir  /mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings \\
    --obs-dir      /mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance \\
    --cfe-dir      /mnt/disk2/suma_helen_poster/cfe_py \\
    --config-file  /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \\
    --param-bounds /mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json \\
    --out-dir      /mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt \\
    --test-forcing-dir1 /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings \\
    --test-forcing-dir2 /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings \\
    --enkf-enabled --enkf-members 20

  Use the with-variance obs dir so the per-hour kriging variance enters R(t)
  via the Vrugt + kriging-scaling formula. Pass --no-vrugt-r to disable the
  Vrugt R and fall back to raw kriging variance (not recommended).

NOTE: pre-stage the calibrated <cat-id>_best_params.json into the out-dir
before running:
  mkdir -p <out-dir>/<cat-id>
  cp <run3-dir>/<cat-id>/<cat-id>_best_params.json <out-dir>/<cat-id>/

To run all 21 catchments, launch in batches of ~5 (each catchment uses N CFE
instances simultaneously; 21 × 20 = 420 processes if run all at once).
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ALBEDO   = 0.20
ALPHA_PT = 1.26

TIME_SPLIT = {
    "spinup-for-calibration": {
        "start": "2019-01-01 00:00:00",
        "end":   "2019-12-31 23:00:00",
    },
    "calibration": {
        "start": "2020-01-01 00:00:00",
        "end":   "2022-12-31 23:00:00",
    },
    "spinup-for-testing": {
        "start": "2023-02-01 00:00:00",
        "end":   "2023-09-30 23:00:00",
    },
    "testing": {
        "start": "2023-10-01 00:00:00",
        "end":   "2024-10-31 23:00:00",
    },
}

# Set at runtime
CAT_ID            = None
FORCING_FILE      = None
TEST_FORCING_FILE = None
OBS_FILE          = None
CFE_CONFIG_FILE   = None
PARAM_BOUNDS_FILE = None
OUT_DIR           = None
bmi_cfe           = None
ENKF_ENABLED      = False
ENKF_CONFIG       = {}


class EnKFAssimilator:
    """
    True stochastic Ensemble Kalman Filter (Burgers / van Leeuwen / Evensen 1998)
    for CFE state updates from per-catchment kriging observations.

    `update_states(models, ...)` runs the N-member ensemble update:
    Kalman gain per state derived from cross-covariance Pxy / (Pyy + R),
    observations perturbed per member, each member updated independently.

    Updates 4 CFE states per member (soil, GW, Nash[0], Nash[1]) with a
    mass-conserving overflow/underflow cascade so corrections are redirected
    along CFE's flow path instead of being silently clipped.
    """
    
    def __init__(self, n_members=20, obs_error_std=0.05, obs_file=None,
                 precip_perturb_frac=0.15, pet_perturb_frac=0.10,
                 init_state_perturb_frac=0.05,
                 soil_process_noise_frac=0.002,
                 gw_process_noise_frac=0.0015,
                 nash_process_noise_frac=0.005,
                 use_vrugt_r=True, vrugt_alpha=0.10, vrugt_scale=0.001,
                 rng_seed=None):
        """
        Initialize EnKF assimilator.

        Args:
            n_members (int): Number of ensemble members (>=2 for true EnKF)
            obs_error_std (float): Default obs std dev when CSV lacks 'variance'
            obs_file (str): Path to kriging observations CSV
            precip_perturb_frac (float): Lognormal sigma for precip noise (mean=1)
            pet_perturb_frac (float): Multiplicative Gaussian noise std for PET
            init_state_perturb_frac (float): Multiplicative noise std on init states
            soil_process_noise_frac (float): Per-hour multiplicative process noise std
                for soil reservoir. Default 0.002 (0.2%) — soil is slow-evolving so a
                small noise is sufficient and avoids cumulative-random-walk drift.
                Matches NWC-CUAHSI reference implementation.
            gw_process_noise_frac (float): Per-hour multiplicative process noise std
                for GW reservoir. Default 0.0015 (0.15%) — slowest state, smallest noise.
            nash_process_noise_frac (float): Per-hour multiplicative process noise std
                for Nash[0] and Nash[1]. Default 0.005 (0.5%) — needed larger than
                soil/GW because Nash spread does not develop naturally (buckets start
                at 0 in CFE init). Combined with a small additive floor so noise
                survives when the bucket is empty.
            use_vrugt_r (bool): If True (default, production setting), compute
                observation error variance R as a Vrugt 2005 (SODA) heteroscedastic
                function of flow magnitude, scaled by the kriging variance:
                R(t) = (alpha * y_obs)^2 + scale * sigma^2_krig. Fixes the case
                where raw kriging variance is too large to allow DA to fire, and
                gives the most uniform performance across catchments. Set to False
                to revert to the raw kriging variance / fallback std behavior.
            vrugt_alpha (float): Relative-error fraction in the Vrugt R. Default 0.10
                (10% of flow). Standard hydrology DA value (Vrugt et al. 2005).
            vrugt_scale (float): Scaling on the kriging-variance term. Default 0.001.
                Picked so the kriging contribution is the same order of magnitude as
                the flow-magnitude term at typical flows for this basin.
            rng_seed (int or None): Seed for the ensemble random generator
        """
        self.n_members = n_members
        self.obs_error_std = obs_error_std
        self.obs_error_var = obs_error_std ** 2
        self.precip_perturb_frac     = precip_perturb_frac
        self.pet_perturb_frac        = pet_perturb_frac
        self.init_state_perturb_frac = init_state_perturb_frac
        self.soil_process_noise_frac = soil_process_noise_frac
        self.gw_process_noise_frac   = gw_process_noise_frac
        self.nash_process_noise_frac = nash_process_noise_frac
        self.use_vrugt_r             = use_vrugt_r
        self.vrugt_alpha             = vrugt_alpha
        self.vrugt_scale             = vrugt_scale
        self.rng = np.random.default_rng(rng_seed)
        
        # Load kriging observations
        # Handles three obs-CSV layouts:
        #   Format A:                     columns = 'datetime', 'qkrig', optional 'variance'
        #   Format B (no variance):       columns = 'time', 'qkrig_mm_hr'
        #   Format B with-variance:       columns = 'time', 'qkrig_mm_hr', 'qkrig_variance'
        self.obs_dict = {}
        self.obs_var_dict = {}
        if obs_file and os.path.exists(obs_file):
            obs_df = pd.read_csv(obs_file)
            if 'time' in obs_df.columns and 'qkrig_mm_hr' in obs_df.columns:
                obs_df = obs_df.rename(columns={'time': 'datetime', 'qkrig_mm_hr': 'qkrig'})
            if 'qkrig_variance' in obs_df.columns:
                obs_df = obs_df.rename(columns={'qkrig_variance': 'variance'})
            obs_df['date'] = pd.to_datetime(obs_df['datetime']).dt.strftime('%Y-%m-%d %H:%M:%S')
            has_variance = 'variance' in obs_df.columns
            for _, row in obs_df.iterrows():
                y_obs    = row['qkrig']
                krig_var = row['variance'] if has_variance else self.obs_error_var
                self.obs_dict[row['date']] = y_obs
                # Vrugt et al. 2005 (SODA) heteroscedastic R, scaled by kriging
                # Vrugt 2005 heteroscedastic R scaled by kriging variance:
                #   R(t) = (alpha * y_obs(t))^2 + scale * sigma^2_krig(t)
                # Flow-magnitude term makes R small at low flow (so DA fires),
                # larger at peaks (where obs is also more uncertain). Kriging
                # variance brings per-hour catchment-specific info in, scaled
                # down so it doesn't dominate.
                if self.use_vrugt_r and not pd.isna(y_obs):
                    self.obs_var_dict[row['date']] = (
                        (self.vrugt_alpha * y_obs) ** 2
                        + self.vrugt_scale * krig_var
                    )
                else:
                    self.obs_var_dict[row['date']] = krig_var
        
        self.n_updates = 0
        self.total_increment = 0.0
        self.total_overflow_lost_mm = 0.0  # water lost at GW underflow (mass-conservation breach)

        # True-EnKF diagnostics (test path)
        self.avg_pyy = 0.0
        self.avg_K_sm = 0.0
        self.avg_K_gw = 0.0
        self.avg_K_n0 = 0.0
        self.avg_K_n1 = 0.0

    # ------------------------------------------------------------------
    # Ensemble helpers (used by the test loop)
    # ------------------------------------------------------------------
    def perturb_forcing(self, precip_mm_h, pet_mm_h):
        """Return arrays of length n_members with perturbed (P, PET).

        Precip uses **lognormal** noise (literature standard for precip in
        hydrology DA — Clark 2008, Renard 2010). Guarantees non-negative
        precip and matches the heavy-tailed empirical distribution of
        precip errors. Lognormal is parameterized so the multiplier has
        E[noise] = 1 (mu = -sigma²/2).

        PET uses multiplicative Gaussian (more symmetric distribution).
        """
        N = self.n_members
        sigma_p = self.precip_perturb_frac
        mu_p    = -0.5 * sigma_p ** 2  # makes E[exp(mu + sigma·Z)] = 1
        p = precip_mm_h * self.rng.lognormal(mu_p, sigma_p, N)
        # Lognormal × non-negative is non-negative; no clipping needed.

        e = pet_mm_h * (1.0 + self.pet_perturb_frac * self.rng.standard_normal(N))
        return p, np.maximum(e, 0.0)

    def perturb_initial_state(self, value, upper=None, floor=1e-6):
        """Multiplicatively perturb a single initial-state scalar across N members."""
        v = value * (1.0 + self.init_state_perturb_frac * self.rng.standard_normal(self.n_members))
        v = np.maximum(v, floor)
        if upper is not None:
            v = np.minimum(v, upper)
        return v

    def add_process_noise(self, models):
        """Inject per-timestep state perturbations on every member.

        Process noise (the Q matrix in classical Kalman / "model error" term in
        EnKF) is required by textbook EnKF to prevent ensemble collapse over
        long runs. Without it, the analysis pulls all members toward similar
        states, spread shrinks, Pyy → 0, and the filter goes blind.

        Each state gets its own σ (calibrated by physical response time):
          soil  : 0.2%/hr  — slow-evolving; small noise prevents random-walk drift
          GW    : 0.15%/hr — slowest state
          Nash  : 0.5%/hr  + small additive floor — needs more because Nash
                  spread does not develop naturally (buckets start at 0).
        Magnitudes match NWC-CUAHSI/data_assimilation_with_bmi reference impl.
        """
        sigma_soil = self.soil_process_noise_frac
        sigma_gw   = self.gw_process_noise_frac
        sigma_nash = self.nash_process_noise_frac
        if sigma_soil <= 0 and sigma_gw <= 0 and sigma_nash <= 0:
            return
        for m in models:
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]
            sm0 = m.soil_reservoir["storage_m"]
            gw0 = m.gw_reservoir["storage_m"]
            n00 = float(m.nash_storage[0])
            n10 = float(m.nash_storage[1])
            # Additive floor for nash so noise survives at zero.
            nash_floor = max(sm_max * 1e-4, 1e-7)
            if sigma_soil > 0:
                m.soil_reservoir["storage_m"] = float(np.clip(
                    sm0 * (1.0 + sigma_soil * self.rng.standard_normal()),
                    0.0, sm_max))
            if sigma_gw > 0:
                m.gw_reservoir["storage_m"] = float(np.clip(
                    gw0 * (1.0 + sigma_gw * self.rng.standard_normal()),
                    0.0, gw_max))
            if sigma_nash > 0:
                m.nash_storage[0] = max(
                    n00 * (1.0 + sigma_nash * self.rng.standard_normal())
                    + nash_floor * self.rng.standard_normal(), 0.0)
                m.nash_storage[1] = max(
                    n10 * (1.0 + sigma_nash * self.rng.standard_normal())
                    + nash_floor * self.rng.standard_normal(), 0.0)

    # ------------------------------------------------------------------
    # True stochastic EnKF (Burgers / van Leeuwen / Evensen 1998)
    # Used by run_testing_period — replaces 30/15/20/35 + 0.3-heuristic with
    # ensemble-derived Kalman gains.
    # ------------------------------------------------------------------
    def update_states(self, models, current_date, forecast_runoffs_mm_h):
        """
        Perform a true stochastic EnKF update across N ensemble members.

        Args:
            models: list of N CFE BMI model instances (one per ensemble member)
            current_date (str): 'YYYY-MM-DD HH:MM:SS'
            forecast_runoffs_mm_h: array-like of N forecast runoffs (mm/h)

        Returns:
            dict with update statistics (mean innovation, Pyy, K per state).
        """
        stats = {"updated": False, "innovation_mean": 0.0, "pyy": 0.0,
                 "K_sm": 0.0, "K_gw": 0.0, "K_n0": 0.0, "K_n1": 0.0}

        if current_date not in self.obs_dict:
            return stats
        obs = self.obs_dict[current_date]
        obs_var = self.obs_var_dict[current_date]
        if pd.isna(obs):
            return stats

        q = np.asarray(forecast_runoffs_mm_h, dtype=float)
        if np.any(np.isnan(q)):
            return stats

        N = self.n_members
        if N < 2:
            return stats  # need ≥2 members for a meaningful covariance

        # ---- Collect ensemble of 4 states (meters) ----
        sm = np.array([m.get_value('SOIL_CONCEPTUAL_STORAGE') for m in models], dtype=float)
        gw = np.array([m.gw_reservoir["storage_m"] for m in models], dtype=float)
        n0 = np.array([float(m.nash_storage[0]) for m in models], dtype=float)
        n1 = np.array([float(m.nash_storage[1]) for m in models], dtype=float)

        # ---- Ensemble means and anomalies ----
        sm_a = sm - sm.mean()
        gw_a = gw - gw.mean()
        n0_a = n0 - n0.mean()
        n1_a = n1 - n1.mean()
        q_mean = q.mean()
        q_a    = q - q_mean

        # Pyy = ensemble variance of Q (mm/h)^2 ; degenerate if all members agree
        Pyy = float((q_a * q_a).sum() / (N - 1))
        denom = Pyy + obs_var
        if denom < 1e-12:
            return stats

        # Pxy = state/Q cross-covariance (units: m × mm/h)
        Pxy_sm = float((sm_a * q_a).sum() / (N - 1))
        Pxy_gw = float((gw_a * q_a).sum() / (N - 1))
        Pxy_n0 = float((n0_a * q_a).sum() / (N - 1))
        Pxy_n1 = float((n1_a * q_a).sum() / (N - 1))

        K_sm = Pxy_sm / denom   # units: m / (mm/h)
        K_gw = Pxy_gw / denom
        K_n0 = Pxy_n0 / denom
        K_n1 = Pxy_n1 / denom

        # Perturb the observation N times (Burgers/Evensen)
        obs_pert = obs + np.sqrt(max(obs_var, 0.0)) * self.rng.standard_normal(N)

        # Per-member innovation (mm/h); state increments (m)
        innov = obs_pert - q
        sm_delta = K_sm * innov
        gw_delta = K_gw * innov
        n0_delta = K_n0 * innov
        n1_delta = K_n1 * innov

        # ---- Apply per member with overflow-aware cascade ----
        for i, m in enumerate(models):
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]

            sm_raw = sm[i] + sm_delta[i]
            gw_raw = gw[i] + gw_delta[i]
            n0_raw = n0[i] + n0_delta[i]
            n1_raw = n1[i] + n1_delta[i]

            # Positive overflow cascade: soil→Nash[0], GW→Nash[1]
            so = max(sm_raw - sm_max, 0.0); sm_raw -= so; n0_raw += so
            go = max(gw_raw - gw_max, 0.0); gw_raw -= go; n1_raw += go

            # Negative underflow cascade: Nash[1]→Nash[0]→soil→GW→loss
            n1u = max(-n1_raw, 0.0); n1_raw += n1u; n0_raw -= n1u
            n0u = max(-n0_raw, 0.0); n0_raw += n0u; sm_raw -= n0u
            smu = max(-sm_raw, 0.0); sm_raw += smu; gw_raw -= smu
            gwu = max(-gw_raw, 0.0); gw_raw += gwu
            self.total_overflow_lost_mm += gwu * 1000.0

            if (np.isnan(sm_raw) or np.isnan(gw_raw)
                    or np.isnan(n0_raw) or np.isnan(n1_raw)):
                continue

            sm_new = float(np.clip(sm_raw, 0.0, sm_max))
            gw_new = float(np.clip(gw_raw, 0.0, gw_max))
            n0_new = float(np.clip(n0_raw, 0.0, None))
            n1_new = float(np.clip(n1_raw, 0.0, None))

            m.set_value('SOIL_CONCEPTUAL_STORAGE', sm_new)
            m.gw_reservoir["storage_m"] = gw_new
            m.nash_storage[0] = n0_new
            m.nash_storage[1] = n1_new

        # ---- Diagnostics ----
        self.n_updates += 1
        innov_mean = float(innov.mean())
        self.total_increment += abs(innov_mean)  # repurpose for ensemble: mean innov magnitude
        # Running averages of Pyy and Kalman gains for end-of-run reporting
        a = self.n_updates
        self.avg_pyy  = ((a - 1) * self.avg_pyy  + Pyy)  / a
        self.avg_K_sm = ((a - 1) * self.avg_K_sm + K_sm) / a
        self.avg_K_gw = ((a - 1) * self.avg_K_gw + K_gw) / a
        self.avg_K_n0 = ((a - 1) * self.avg_K_n0 + K_n0) / a
        self.avg_K_n1 = ((a - 1) * self.avg_K_n1 + K_n1) / a

        stats.update(updated=True, innovation_mean=innov_mean, pyy=Pyy,
                     K_sm=K_sm, K_gw=K_gw, K_n0=K_n0, K_n1=K_n1)

        # Log sparingly — every hour is noisy. Print at storms (large Pyy).
        if Pyy > 1e-4 or self.n_updates % 200 == 0:
            print(f"  [EnKF N={N}] {current_date} | obs={obs:.3f} | q_mean={q_mean:.3f} | "
                  f"innov_mean={innov_mean:+.4f} | Pyy={Pyy:.6f} | "
                  f"K_sm={K_sm:+.3e} K_gw={K_gw:+.3e} K_n0={K_n0:+.3e} K_n1={K_n1:+.3e}")

        return stats



def priestley_taylor_pet(srad_wm2, T_kelvin, alpha=ALPHA_PT):
    """PT-PET from shortwave radiation + temperature. Returns mm/h."""
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    pet = alpha * (delta / (delta + gamma)) * Rn_mj / lam
    return np.maximum(pet, 0.0)


def load_forcing(path=None):
    """Load NWM forcing CSV (retro or operational) and add PET column. Returns DataFrame."""
    path = path or FORCING_FILE
    df = pd.read_csv(path)
    df = df.rename(columns={"time": "date", "APCP_surface": "total_precipitation"})
    df['date'] = pd.to_datetime(df['date']).dt.strftime('%Y-%m-%d %H:%M:%S')
    df['potential_evaporation'] = priestley_taylor_pet(
        df['DSWRF_surface'].values,
        df['TMP_2maboveground'].values
    )
    return df


def load_test_forcing():
    """Load NWM operational forcing CSV. APCP_surface is in kg/m²/s → convert to mm/h."""
    df = pd.read_csv(TEST_FORCING_FILE)
    df['date'] = pd.to_datetime(df['time']).dt.strftime('%Y-%m-%d %H:%M:%S')
    df['total_precipitation'] = df['APCP_surface'] * 3600.0  # kg/m²/s → mm/h
    df['potential_evaporation'] = priestley_taylor_pet(
        df['DSWRF_surface'].values,
        df['TMP_2maboveground'].values
    )
    return df


def _kge_score(obs, sim):
    """Kling-Gupta Efficiency on the overlap of obs and sim (NaNs skipped)."""
    m = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[m], sim[m]
    if len(o) < 2:
        return float('nan')
    denom = float(np.sqrt(((o - o.mean()) ** 2).sum() * ((s - s.mean()) ** 2).sum()))
    if denom == 0:
        return float('nan')
    r = float(((o - o.mean()) * (s - s.mean())).sum() / denom)
    alpha = float(s.std() / o.std()) if o.std() != 0 else float('nan')
    beta = float(s.mean() / o.mean()) if o.mean() != 0 else float('nan')
    return 1.0 - float(np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def _nse_score(obs, sim):
    """Nash-Sutcliffe Efficiency on the overlap of obs and sim (NaNs skipped)."""
    m = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[m], sim[m]
    if len(o) < 2:
        return float('nan')
    denom = ((o - o.mean()) ** 2).sum()
    if denom == 0:
        return float('nan')
    return 1.0 - float(((o - s) ** 2).sum() / denom)


def run_testing_period(best_param_dict):
    """Run CFE over test period using NWM operational forcing with optional EnKF-DA."""

    def custom_load_forcing(self_cfe):
        df = load_test_forcing()
        self_cfe.forcing_data = df.rename(columns={"date": "time"})

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

    tmp_cfg = str(OUT_DIR / f'{CAT_ID}_bmi_config_temp_test.json')
    with open(tmp_cfg, 'w') as f:
        json.dump(cfg, f)

    # Initialize EnKF assimilator first so its RNG/perturbation knobs are available
    enkf_test = None
    N_members = 1
    if ENKF_ENABLED:
        enkf_test = EnKFAssimilator(
            n_members=ENKF_CONFIG.get('n_members', 20),
            obs_error_std=ENKF_CONFIG.get('obs_error_std', 0.05),
            obs_file=OBS_FILE,
            use_vrugt_r=ENKF_CONFIG.get('use_vrugt_r', True),
            vrugt_alpha=ENKF_CONFIG.get('vrugt_alpha', 0.10),
            vrugt_scale=ENKF_CONFIG.get('vrugt_scale', 0.001),
        )
        N_members = enkf_test.n_members
        if N_members < 2:
            raise ValueError(f"--enkf-members must be >=2 for true EnKF (got {N_members})")
        print(f"[EnKF] Initializing {N_members} ensemble members for {CAT_ID}")

    # Build N CFE model instances. With EnKF, each member starts from a slightly
    # perturbed initial state so the ensemble has spread from t=0.
    models = []
    for i in range(N_members):
        m = bmi_cfe.BMI_CFE(cfg_file=tmp_cfg)
        m.load_forcing_file = custom_load_forcing.__get__(m)
        m.initialize()
        if ENKF_ENABLED and i > 0:
            sm_max = m.soil_reservoir["storage_max_m"]
            gw_max = m.gw_reservoir["storage_max_m"]
            sm0 = m.soil_reservoir["storage_m"]
            gw0 = m.gw_reservoir["storage_m"]
            n00 = float(m.nash_storage[0])
            n10 = float(m.nash_storage[1])
            r = enkf_test.init_state_perturb_frac
            m.soil_reservoir["storage_m"] = float(np.clip(
                sm0 * (1 + r * enkf_test.rng.standard_normal()), 1e-6, sm_max))
            m.gw_reservoir["storage_m"]   = float(np.clip(
                gw0 * (1 + r * enkf_test.rng.standard_normal()), 1e-6, gw_max))
            # Nash buckets often start at 0 → use small additive jitter instead of multiplicative
            jitter = max(sm_max * 1e-6, 1e-9)
            m.nash_storage[0] = max(n00 + jitter * enkf_test.rng.standard_normal(), 0.0)
            m.nash_storage[1] = max(n10 + jitter * enkf_test.rng.standard_normal(), 0.0)
        models.append(m)

    df = load_test_forcing()

    # ----- Spinup (with perturbed forcing per member if ensemble) -----
    sp_mask = (df['date'] >= TIME_SPLIT['spinup-for-testing']['start']) & \
              (df['date'] <= TIME_SPLIT['spinup-for-testing']['end'])
    df_sp = df[sp_mask]
    for p, e in zip(df_sp['total_precipitation'], df_sp['potential_evaporation']):
        if ENKF_ENABLED:
            p_arr, e_arr = enkf_test.perturb_forcing(p, e)
        for i, m in enumerate(models):
            p_i = float(p_arr[i]) if ENKF_ENABLED else p
            e_i = float(e_arr[i]) if ENKF_ENABLED else e
            m.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux', p_i / 1000)
            m.set_value('water_potential_evaporation_flux',                            e_i / 1000 / 3600)
            m.update()
        # Process noise during spinup keeps spread alive before DA starts
        if ENKF_ENABLED:
            enkf_test.add_process_noise(models)

    # ----- Test period -----
    t_mask = (df['date'] >= TIME_SPLIT['testing']['start']) & \
             (df['date'] <= TIME_SPLIT['testing']['end'])
    df_test = df[t_mask]
    outputs = models[0].get_output_var_names()
    out_lists = {o: [] for o in outputs}

    for p, e, current_date in zip(df_test['total_precipitation'],
                                   df_test['potential_evaporation'],
                                   df_test['date']):
        if ENKF_ENABLED:
            p_arr, e_arr = enkf_test.perturb_forcing(p, e)

        # Run each member forward one hour, collect ensemble Q_sim
        ensemble_q_mm_h = np.empty(N_members, dtype=float)
        for i, m in enumerate(models):
            p_i = float(p_arr[i]) if ENKF_ENABLED else p
            e_i = float(e_arr[i]) if ENKF_ENABLED else e
            m.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux', p_i / 1000)
            m.set_value('water_potential_evaporation_flux',                            e_i / 1000 / 3600)
            m.update()
            ensemble_q_mm_h[i] = m.get_value('land_surface_water__runoff_depth') * 1000  # m/h → mm/h

        # True ensemble EnKF state update
        if ENKF_ENABLED and enkf_test is not None:
            enkf_test.update_states(models, current_date, ensemble_q_mm_h)

        # Process noise on analyzed states — keeps ensemble spread alive across
        # the long test run so the filter doesn't collapse (Pyy → 0).
        if ENKF_ENABLED and enkf_test is not None:
            enkf_test.add_process_noise(models)

        # Recorded Q_sim at time t is the FORECAST (pre-analysis): it was computed
        # by model.update() above, before this hour's DA touched the states. DA
        # affects the recorded series only from t+1 onward, via the analyzed states
        # carried into the next iteration. This is standard sequential-filter
        # forecast verification — comparing Q_sim(t) to obs(t) is not circular.
        for o in outputs:
            vals = [m.get_value(o) for m in models]
            out_lists[o].append(float(np.mean(vals)))

    for m in models:
        m.finalize()

    sim_test = np.array(out_lists['land_surface_water__runoff_depth']) * 1000  # m/h → mm/h
    test_dates = pd.to_datetime(df_test['date'].values)

    # Load kriging obs for test period
    obs_raw = pd.read_csv(OBS_FILE)
    if 'time' in obs_raw.columns and 'qkrig_mm_hr' in obs_raw.columns:
        obs_raw = obs_raw.rename(columns={'time': 'datetime', 'qkrig_mm_hr': 'qkrig'})
    obs_raw['date'] = pd.to_datetime(obs_raw['datetime']).dt.strftime('%Y-%m-%d %H:%M:%S')
    obs_raw = obs_raw.rename(columns={"qkrig": "obs_mm_h"})
    test_dates_df = pd.DataFrame({'date': test_dates.strftime('%Y-%m-%d %H:%M:%S')})
    merged = test_dates_df.merge(obs_raw[['date', 'obs_mm_h']], on='date', how='left')
    obs_test = merged['obs_mm_h'].values

    kge_test = _kge_score(obs_test, sim_test)
    nse_test = _nse_score(obs_test, sim_test)
    
    da_suffix = " | EnKF-DA" if ENKF_ENABLED else ""
    print(f"Test KGE: {kge_test:.4f} | Test NSE: {nse_test:.4f}{da_suffix}")

    df_out = pd.DataFrame({
        'date':        test_dates.strftime('%Y-%m-%d %H:%M:%S'),
        'sim_mm_h':    sim_test,
        'obs_mm_h':    obs_test,
        'precip_mm_h': df_test['total_precipitation'].values,
    })
    df_out.to_csv(OUT_DIR / f'{CAT_ID}_test_results.csv', index=False)

    # Plot test period
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10))

    ax1.plot(test_dates, sim_test, 'tomato', lw=1.5, label='simulated')
    ax1.plot(test_dates, obs_test, 'k', lw=1, label='observed (kriging)')
    da_label = " (EnKF-DA)" if ENKF_ENABLED else ""
    ax1.set_ylabel('Discharge (mm/h)')
    ax1.set_title(f'{CAT_ID} | Test (Oct 2023–Oct 2024){da_label} | KGE={kge_test:.4f} | NSE={nse_test:.4f}')
    ax1.legend()
    ax1_twin = ax1.twinx()
    ax1_twin.plot(test_dates, df_test['total_precipitation'].values, 'steelblue', lw=0.8, alpha=0.5)
    ax1_twin.set_ylim([80, 0])
    ax1_twin.set_ylabel('Precip (mm/h)')

    helene = (test_dates >= pd.Timestamp('2024-09-20')) & (test_dates <= pd.Timestamp('2024-10-05'))
    ax2.plot(test_dates[helene], sim_test[helene], 'tomato', lw=2, label='simulated')
    ax2.plot(test_dates[helene], obs_test[helene], 'k', lw=1.5, label='observed')
    ax2.set_ylabel('Discharge (mm/h)')
    ax2.set_title('Helene Zoom: Sep 20 – Oct 5, 2024')
    ax2.legend()
    ax2_twin = ax2.twinx()
    ax2_twin.bar(test_dates[helene], df_test['total_precipitation'].values[helene],
                 color='steelblue', alpha=0.4, width=0.04)
    ax2_twin.set_ylim([50, 0])
    ax2_twin.set_ylabel('Precip (mm/h)')

    plt.tight_layout()
    plt.savefig(OUT_DIR / f'{CAT_ID}_test_plot.png', dpi=150, bbox_inches='tight')
    plt.close()

    if ENKF_ENABLED and enkf_test is not None:
        print(f"[EnKF Stats - Test] N={enkf_test.n_members} | updates: {enkf_test.n_updates} | "
              f"mean |innov|: {enkf_test.total_increment/max(enkf_test.n_updates,1):.4f} mm/h | "
              f"avg Pyy: {enkf_test.avg_pyy:.6f} | "
              f"avg K_sm={enkf_test.avg_K_sm:+.3e} K_gw={enkf_test.avg_K_gw:+.3e} "
              f"K_n0={enkf_test.avg_K_n0:+.3e} K_n1={enkf_test.avg_K_n1:+.3e} | "
              f"Mass lost at GW underflow: {enkf_test.total_overflow_lost_mm:.4f} mm")

    return kge_test, nse_test


def main():
    global CAT_ID, FORCING_FILE, TEST_FORCING_FILE, OBS_FILE
    global CFE_CONFIG_FILE, PARAM_BOUNDS_FILE, OUT_DIR, bmi_cfe
    global ENKF_ENABLED, ENKF_CONFIG

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',            required=True,  help='Catchment ID, e.g. cat-1016300')
    parser.add_argument('--forcing-dir',        required=True,  help='Dir with per-catchment NWM retro CSVs')
    parser.add_argument('--obs-dir',            required=True,  help='Dir with per-catchment kriging obs CSVs')
    parser.add_argument('--cfe-dir',            required=True,  help='Path to cfe_py directory')
    parser.add_argument('--config-file',        required=True,  help='Base CFE BMI config JSON')
    parser.add_argument('--param-bounds',       required=True,  help='Parameter bounds JSON')
    parser.add_argument('--out-dir',            required=True,  help='Output directory')
    parser.add_argument('--test-forcing-dir1',  default=None,   help='NWM operational forcings dir 1 (2023-Feb2024)')
    parser.add_argument('--test-forcing-dir2',  default=None,   help='NWM operational forcings dir 2 (Feb2024-2025)')
    parser.add_argument('--enkf-enabled',       action='store_true', help='Enable EnKF-based Data Assimilation')
    parser.add_argument('--enkf-members',       type=int, default=20,   help='Number of ensemble members (default: 20)')
    parser.add_argument('--enkf-obs-error-std', type=float, default=0.05, help='Observation error std dev in mm/h (default: 0.05)')
    # Vrugt R is ON by default in the production config. Pass --no-vrugt-r to disable.
    parser.add_argument('--no-vrugt-r',         action='store_true',
                        help='Disable Vrugt 2005 heteroscedastic R (revert to raw kriging variance)')
    parser.add_argument('--vrugt-alpha',        type=float, default=0.10,
                        help='Relative-error fraction in Vrugt R (default: 0.10)')
    parser.add_argument('--vrugt-scale',        type=float, default=0.001,
                        help='Scaling on kriging-variance term in Vrugt R (default: 0.001)')
    args = parser.parse_args()

    CAT_ID            = args.cat_id
    FORCING_FILE      = os.path.join(args.forcing_dir, f'{CAT_ID}.csv')
    OBS_FILE          = os.path.join(args.obs_dir,     f'{CAT_ID}.csv')
    CFE_CONFIG_FILE   = args.config_file
    PARAM_BOUNDS_FILE = args.param_bounds
    OUT_DIR           = Path(args.out_dir) / CAT_ID
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Set EnKF configuration
    ENKF_ENABLED = args.enkf_enabled
    ENKF_CONFIG = {
        'n_members': args.enkf_members,
        'obs_error_std': args.enkf_obs_error_std,
        'use_vrugt_r': not args.no_vrugt_r,   # default ON; --no-vrugt-r disables
        'vrugt_alpha': args.vrugt_alpha,
        'vrugt_scale': args.vrugt_scale,
    }

    if ENKF_ENABLED:
        print(f"\n{'='*80}")
        print(f"EnKF Data Assimilation ENABLED")
        print(f"  Ensemble members: {ENKF_CONFIG['n_members']}")
        print(f"  Observation error std dev: {ENKF_CONFIG['obs_error_std']} mm/h")
        if ENKF_CONFIG['use_vrugt_r']:
            print(f"  Vrugt R: alpha={ENKF_CONFIG['vrugt_alpha']}, scale={ENKF_CONFIG['vrugt_scale']}")
        else:
            print(f"  Vrugt R: disabled (using raw kriging variance / fallback)")
        print(f"{'='*80}\n")

    # Build combined test forcing file if both dirs provided
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
            print(f"Warning: test forcing files not found for {CAT_ID}, skipping test period")

    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    best_params_file = OUT_DIR / f'{CAT_ID}_best_params.json'

    if not TEST_FORCING_FILE or not os.path.exists(TEST_FORCING_FILE):
        print("No test forcing available; nothing to do. Pre-stage the operational "
              "forcing dirs via --test-forcing-dir1/--test-forcing-dir2.")
        return

    if not best_params_file.exists():
        print(f"No best params file at {best_params_file}. Pre-stage the Run 3 "
              f"best_params.json into the out-dir before running.")
        return

    with open(best_params_file) as f:
        saved = json.load(f)
    best_param_dict = saved['best_parameters']
    print(f"\nLoaded calibrated parameters from {best_params_file}")
    print(f"Running test period (Oct 2023 – Oct 2024)...")
    run_testing_period(best_param_dict)


if __name__ == '__main__':
    main()
