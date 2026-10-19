# Folder 4 — Dynamic Variance Direct

## R Formula
```
R(t) = σ²_krig(t)
```
Raw kriging variance used directly as observation error — no flow-magnitude scaling
and no Vrugt formula. RNG seed = 42.

## Design Differences vs Other Folders

| | Folder 1 | Folder 2 | Folder 3 | Folder 4 |
|---|---|---|---|---|
| R formula | Vrugt | Fixed 0.07 | Vrugt | **Direct σ²** |
| RNG seed | hash(cat_id) | hash(cat_id) | 42 | **42** |
| σ²_krig source | static ~4.17 | static ~4.17 | dynamic | **dynamic** |

> **Controlled comparison:** F4 shares the same RNG seed (42), same obs dataset
> (`spliced_dyn_helene`), and same 21 catchments as F3. The only difference between
> F3 and F4 is the R formula — making **F3 vs F4 the clean controlled comparison**
> for this experimental setup.
>
> F4 cannot be cleanly compared to F1 or F2 because both the seed and the R formula
> differ. See the main [README](../README.md) for full details.

## Why F4 Performs Worst

With `R = σ²_krig`, the Kalman gain is:
```
K = P / (P + σ²_krig)
```
During Helene, `σ²_krig` is dynamic (~1.6–2.5 at peak in spliced obs, vs ~4.17
static). Even with lower variance, `σ²_krig` is still much larger than the fixed
R=0.07 used by F2. This keeps K lower than F2 throughout the flood, weakening DA
updates at the most critical timesteps. The absence of the flow-magnitude term means
there is no systematic relationship between R and observed flow — the DA strength is
driven entirely by the spatial kriging uncertainty, which is not well-correlated with
the hydrological need for correction.

## Key Results

| Metric | Value |
|---|---|
| Full period KGE | **+0.200** |
| Full period NSE | **+0.428** |
| Helene KGE | **+0.132** |
| Helene NSE | **+0.303** |
| Helene peak (routed) | **667.9 m³/s (35% of USGS 1885.7)** |
| Ensemble p95 at Helene peak | 744 m³/s |
| Ensemble p50 at Helene peak | 552 m³/s |

## Pipeline Status

| Step | Status | Server path |
|---|---|---|
| 1. Calibrate | ✅ shared | `1400_sites_helene/da_results_dynamic_novrugt_seeded/{cat}/{cat}_best_params.json` |
| 2a/2b. Perturbation arms | ✅ | `suma_helen_poster/da_results/da_arms_dynamic_novrugt_seeded/` |
| 2c. 18hr forecast cycles | ✅ | `suma_helen_poster/da_results/da_forecast_dynamic_novrugt_seeded/` |
| 2d. 600-member ensemble | ✅ | `suma_helen_poster/da_results/da_crossed_dynamic_novrugt_seeded/` |
| 3. Route analysis | ✅ | `suma_helen_poster/da_results/dynamic_novrugt_seeded_routed/` |
| 3. Route forecast cycles | ✅ | `suma_helen_poster/da_results/da_forecast_dynamic_novrugt_seeded_routed/` |
| 3. Route ensemble | ✅ | `suma_helen_poster/da_results/da_crossed_dynamic_novrugt_seeded_routed/` |
| 4a. Error decay | ✅ | `f4_lead_time_decay_gauge_pooled.png` / `_by_regime.png` |
| 4b. Ensemble vs obs | ✅ | `f4_helene_ensemble_vs_usgs.png` / `f4_helene_ensemble_twopanel.png` |
| 4c. Reconstructed timeseries | ✅ | `f4_reconstructed_timeseries.png` / `_helene.png` |

All server paths are under `/mnt/disk2/` unless prefixed with `1400_sites_helene/`.

## Running This Experiment

```bash
# On the server — all 3 stages
bash 2_assimilation/batch_run_all_f4.sh

# Or stage by stage
bash 2_assimilation/batch_run_all_f4.sh arms      # 2a/2b
bash 2_assimilation/batch_run_all_f4.sh forecast  # 2c
bash 2_assimilation/batch_run_all_f4.sh ensemble  # 2d

# 4c timeseries plots
bash 4_evaluation/4c_timeseries/run_4c_f4.sh
```
