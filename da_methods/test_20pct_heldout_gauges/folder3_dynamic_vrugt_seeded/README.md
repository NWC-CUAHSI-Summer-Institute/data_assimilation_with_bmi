# Folder 3 — Dynamic Vrugt Seeded

## R Formula
```
R(t) = (0.10 × y_obs(t))² + 0.001 × σ²_krig(t)
```
Same Vrugt 2005 heteroscedastic formula as Folder 1, but using dynamic kriging
variance and a fixed RNG seed. See the controlled-comparison note below.

## Design Differences vs Folder 1

| | Folder 1 | Folder 3 |
|---|---|---|
| R formula | Vrugt | Vrugt (same) |
| RNG seed | `hash(cat_id)` | `42` (fixed) |
| σ²_krig source | `catchment_ts_03463300_with_variance` (static ~4.17) | `catchment_ts_03463300_spliced_dyn_helene` (dynamic, ~1.6 at Helene peak) |

> **Controlled comparison caveat:** F3 uses a different RNG seed and a different
> kriging variance dataset than F1. The Qkrig flow observations are identical across
> both datasets — only σ²_krig differs. Since the Vrugt formula weights σ²_krig by
> 0.001, the variance difference has negligible effect on R. The remaining gap in KGE
> between F1 (+0.277) and F3 (+0.261) is explained by the seed difference, not by any
> change in the R formula. **F3 vs F4 is the clean controlled comparison** for this
> obs dataset and seed. See the main [README](../README.md) for full details.

## Key Results

| Metric | Value |
|---|---|
| Full period KGE | **+0.261** |
| Full period NSE | **+0.485** |
| Helene KGE | **+0.189** |
| Helene NSE | **+0.372** |
| Helene peak (routed) | **693.8 m³/s (37% of USGS 1885.7)** |
| Ensemble p95 at Helene peak | 759 m³/s |
| Ensemble p50 at Helene peak | 484 m³/s |

## Pipeline Status

| Step | Status | Server path |
|---|---|---|
| 1. Calibrate | ✅ shared | `1400_sites_helene/da_results_dynamic_vrugt_seeded/{cat}/{cat}_best_params.json` |
| 2a/2b. Perturbation arms | ✅ | `suma_helen_poster/da_results/da_arms_dynamic_vrugt_seeded/` |
| 2c. 18hr forecast cycles | ✅ | `suma_helen_poster/da_results/da_forecast_dynamic_vrugt_seeded/` |
| 2d. 600-member ensemble | ✅ | `suma_helen_poster/da_results/da_crossed_dynamic_vrugt_seeded/` |
| 3. Route analysis | ✅ | `suma_helen_poster/da_results/dynamic_vrugt_seeded_routed/` |
| 3. Route forecast cycles | ✅ | `suma_helen_poster/da_results/da_forecast_dynamic_vrugt_seeded_routed/` |
| 3. Route ensemble | ✅ | `suma_helen_poster/da_results/da_crossed_dynamic_vrugt_seeded_routed/` |
| 4a. Error decay | ✅ | `f3_lead_time_decay_gauge_pooled.png` / `_by_regime.png` |
| 4b. Ensemble vs obs | ✅ | `f3_helene_ensemble_vs_usgs.png` / `f3_helene_ensemble_twopanel.png` |
| 4c. Reconstructed timeseries | ✅ | `f3_reconstructed_timeseries.png` / `_helene.png` |

All server paths are under `/mnt/disk2/` unless prefixed with `1400_sites_helene/`.

## Running This Experiment

```bash
# On the server — all 3 stages
bash 2_assimilation/batch_run_all_f3.sh

# Or stage by stage
bash 2_assimilation/batch_run_all_f3.sh arms      # 2a/2b
bash 2_assimilation/batch_run_all_f3.sh forecast  # 2c
bash 2_assimilation/batch_run_all_f3.sh ensemble  # 2d

# 4c timeseries plots
bash 4_evaluation/4c_timeseries/run_4c_f3.sh
```
