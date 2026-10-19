# Folder 1 — Variance Scaled Using Vrugt Formula

## R Formula
```
R(t) = (0.10 × y_obs(t))² + 0.001 × σ²_krig(t)
```
Observation error variance scales with flow magnitude (Vrugt 2005 heteroscedastic).
σ²_krig from: `/mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance/`

## Pipeline

### 1. Calibrate
Shared best_params.json from external `calibrate-cfe` repo (R formula not involved).

### 2. Assimilation
| Sub-step | Script | Server output |
|---|---|---|
| 2a. Forcing arm (30 members) | `run_perturbation_da_on.py` | `v2_perturbation_da_on/` |
| 2b. Hydro-state arm (20 members) | `run_perturbation_da_on.py` | `v2_perturbation_da_on/` |
| 18hr forecast cycles | `run_lead_time_forecast_sweep.py` | `v2_lead_time_forecast/` |
| 600-member crossed ensemble | `run_crossed_ensemble.py` (no --hardcoded-r) | `v2_crossed_ensemble_vrugt/` |

### 3. Route
| Script | Input | Output |
|---|---|---|
| `run_route.py` | `v2_true_enkf_vrugt/` | `vrugt_dynamic_routed/routed_Q_test.csv` |
| `route_lead_time_forecasts.py` | `v2_lead_time_forecast/` | `v2_lead_time_forecast_routed/` |
| `run_route_crossed_ensemble.py` | `v2_crossed_ensemble_vrugt/` | `v2_crossed_ensemble_vrugt_routed/` |

### 4. Evaluate
| Sub-step | Script | Status |
|---|---|---|
| 4a. Error decay | `plot_forecast_error_fixed_target.py` | ✅ |
| 4b. 600-member ensemble vs obs | `plot_production_ensemble_forecast.py` | ⏳ pending Vrugt ensemble run |
| 4c. Ensemble mean + spread (18hr cycles) | `plot_lead_time_decay_gauge.py` | ✅ |

## Key Results (deterministic analysis routing)
- Full period KGE: **+0.277**
- Helene KGE: **+0.213**
- Helene peak: 752 m³/s (40% of USGS 1885.7 m³/s)
