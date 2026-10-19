# DA Experiment Index

Four experiments comparing observation error variance (R) formulations for CFE + EnKF
data assimilation at USGS gauge 03463300 (South Toe River Near Celo, NC).

Each folder is self-contained with its own `1_calibrate/`, `2_assimilation/`, `3_routing/`,
`4_evaluation/` pipeline. Calibrated parameters (best_params.json) are shared across all
folders from the external `calibrate-cfe` repo — R formula is irrelevant at calibration time.

---

## Results Summary

| Folder | R Formula | Full KGE | Full NSE | Helene KGE | Helene NSE | Helene peak % USGS |
|---|---|---|---|---|---|---|
| [1 — Variance scaled Vrugt](folder1_variance_scaled_vrugt/) | `(0.10·y)² + 0.001·σ²_krig` | +0.277 | **+0.555** | +0.213 | **+0.459** | ~34% |
| [2 — Fixed R=0.07](folder2_fixed_r_007/) | `0.07` (constant) | **+0.503** | +0.163 | **+0.439** | -0.009 | **~65%** |
| [3 — Dynamic Vrugt seeded](folder3_dynamic_vrugt_seeded/) | `(0.10·y)² + 0.001·σ²_krig` (seed=42, spliced obs) | +0.261 | +0.485 | +0.189 | +0.372 | ~37% |
| [4 — Dynamic variance direct](folder4_dynamic_variance_direct/) | `σ²_krig` directly (seed=42) | +0.200 | +0.428 | +0.132 | +0.303 | ~35% |

USGS Helene peak: 1885.7 m³/s

**Finding:** Fixed R=0.07 (Folder 2) gives the best KGE and captures 65% of the Helene
peak — constant small R keeps Kalman gain high throughout the flood. The Vrugt formula
inflates R at high flows, dampening DA updates exactly when they matter most, so Folder 1
captures only 34% of the peak despite having the best NSE. Raw σ²_krig (Folder 4) is
worst across all metrics.

---

## Completeness Matrix

| Step | Folder 1 | Folder 2 | Folder 3 | Folder 4 |
|---|---|---|---|---|
| 1. Calibrate | ✅ shared | ✅ shared | ✅ shared | ✅ shared |
| 2a. Forcing arm (30 members) | ✅ | ✅ (same as F1) | ✅ | ✅ |
| 2b. Hydro-state arm (20 members) | ✅ | ✅ (same as F1) | ✅ | ✅ |
| 2c. 18hr forecast cycles | ✅ | ✅ | ✅ | ✅ |
| 2d. 600-member crossed ensemble | ✅ | ✅ | ✅ | ✅ |
| 3. Route analysis trajectory | ✅ | ✅ | ✅ | ✅ |
| 3. Route 18hr forecast cycles | ✅ | ✅ | ✅ | ✅ |
| 3. Route crossed ensemble | ✅ | ✅ | ✅ | ✅ |
| 4a. Error decay | ✅ | ✅ | ✅ | ✅ |
| 4b. Ensemble vs obs | ✅ | ✅ | ✅ | ✅ |
| 4c. Reconstructed timeseries | ✅ | ✅ | ✅ | ✅ |

---

## Server Data Paths

| Folder | Analysis results | Forecast cycles | Ensemble | Routed analysis |
|---|---|---|---|---|
| 1 | `da_results/v2_true_enkf_vrugt/` | `v2_lead_time_forecast/` | `v2_crossed_ensemble_vrugt/` | `vrugt_dynamic_routed/` |
| 2 | `da_results/v2_fixed_r007_analysis/` | `v2_lead_time_forecast_hardcoded_r/` | `v2_crossed_ensemble/` | `fixed_r007_routed/` |
| 3 | `1400_sites_helene/da_results_dynamic_vrugt_seeded/` | `da_results/da_forecast_dynamic_vrugt_seeded/` | `da_results/da_crossed_dynamic_vrugt_seeded/` | `dynamic_vrugt_seeded_routed/` |
| 4 | `1400_sites_helene/da_results_dynamic_novrugt_seeded/` | `da_results/da_forecast_dynamic_novrugt_seeded/` | `da_results/da_crossed_dynamic_novrugt_seeded/` | `dynamic_novrugt_seeded_routed/` |

All paths are under `/mnt/disk2/suma_helen_poster/` unless prefixed with `1400_sites_helene/`.
