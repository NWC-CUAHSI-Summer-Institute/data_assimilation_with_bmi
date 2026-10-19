# Folder 2 — Fixed R = 0.07

## R Formula
```
R = 0.07  (constant, all timesteps)
```
Constant observation error variance — no flow-magnitude or kriging-variance dependence.

## Pipeline

### 1. Calibrate
Shared best_params.json from external `calibrate-cfe` repo (same as all folders).

### 2. Assimilation
| Sub-step | Script | Server output |
|---|---|---|
| 2a. Forcing arm (30 members, R=0.07) | `run_perturbation_da_on.py --hardcoded-r 0.07` | `v2_perturbation_da_on/` |
| 2b. Hydro-state arm (20 members, R=0.07) | `run_perturbation_da_on.py --hardcoded-r 0.07` | `v2_perturbation_da_on/` |
| 18hr forecast cycles | `run_lead_time_forecast_sweep.py --hardcoded-r 0.07` | `v2_lead_time_forecast_hardcoded_r/` |
| 600-member crossed ensemble | `batch_run_crossed_all_cats.sh` (--hardcoded-r 0.07) | `v2_crossed_ensemble/` |

> Note: 2a/2b arm runs default to R=0.07 — same figures as Folder 1 arms.

### 3. Route
| Script | Input | Output | Status |
|---|---|---|---|
| `run_route.py` | `v2_fixed_r007_analysis/` | `fixed_r007_routed/routed_Q_test.csv` | ✅ |
| `route_lead_time_forecasts.py` | `v2_lead_time_forecast_hardcoded_r/` | `v2_lead_time_forecast_hardcoded_r_routed/` | ✅ |
| `run_route_crossed_ensemble.py` | `v2_crossed_ensemble/` | `v2_crossed_ensemble_routed/` | ✅ |

### 4. Evaluate
| Sub-step | Status |
|---|---|
| 4a. Error decay | ✅ `v2_lead_time_forecast_hardcoded_r_routed/error_fixed_target_mean.png` |
| 4b. 600-member ensemble vs obs | ✅ `v2_crossed_ensemble_routed/` |
| 4c. Ensemble mean + spread | ✅ `v2_lead_time_forecast_hardcoded_r_routed/lead_time_reconstructed_timeseries*.png` |

## Key Results (deterministic analysis routing)
- Full period KGE: **+0.503**  NSE: **+0.163**
- Helene KGE: **+0.439**  NSE: **-0.009**
- Helene peak: 1227.1 m³/s (65% of USGS 1885.7 m³/s)

## Data Source
`/mnt/disk2/suma_helen_poster/da_results/v2_fixed_r007_analysis/` (analysis trajectory)
`/mnt/disk2/suma_helen_poster/da_results/fixed_r007_routed/` (routed output)
