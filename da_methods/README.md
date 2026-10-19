# CFE + EnKF Data Assimilation — Hurricane Helene Forecasting

**Gauge:** USGS 03463300 — South Toe River Near Celo, NC  
**Watershed:** 21 NWM catchments · 113.18 km² · outlet reach wb-1016283  
**Event:** Hurricane Helene, September 24–29, 2024 (USGS peak: 1885.7 m³/s)  
**Test period:** October 2023 – October 2024 (13 months)  
**Model:** CFE (Conceptual Functional Equivalent) hydrological model via BMI interface  
**DA method:** Ensemble Kalman Filter (EnKF) assimilating Qkrig streamflow observations  

---

## What Makes This Setup Different

### Distributed CFE — not lumped

CFE is run as **21 independent instances — one per catchment** in the South Toe River
watershed. Each instance has its own calibrated parameters (`best_params.json`), its own
NWM per-catchment forcings, and its own runoff output.

This is different from the lumped CFE configuration used in earlier work (e.g. summer
institute projects), where a single CFE instance represents the entire watershed with
one parameter set and one forcing input. The distributed setup captures spatial
variability in soil, land cover, and catchment response that a lumped model cannot.

### Custom T-route wrapper — not ngen-troute

Two T-route use cases exist in this repo, both using a **custom wrapper** that calls
T-route's Muskingum-Cunge routing function (`compute_network_structured`) directly:

| Approach | What initializes T-route | When used |
|---|---|---|
| **Qkrig → T-route** | Kriged streamflow pseudo-observations | Evaluate kriging directly as routing input |
| **Distributed CFE → T-route** | Per-catchment CFE runoff (mm/h) | Evaluate calibrated + DA-corrected CFE |

Both differ from **ngen-troute** — the standard T-route integration inside the full
NextGen framework — which expects NWM-formatted YAML configs and full hydrofabric
network inputs. The custom wrapper was built because the CFE-BMI DA pipeline outputs
per-catchment CSVs that are incompatible with the ngen-troute CLI interface.

### DA builds on top

The EnKF data assimilation layer (this repo's primary contribution) sits on top of
distributed CFE and the custom T-route wrapper:

```
Layer 1 — Distributed CFE          per-catchment runoff (21 instances)
Layer 2 — Custom T-route wrapper   routes CFE or Qkrig output to outlet
Layer 3 — EnKF DA                  corrects CFE states using Qkrig observations
```

---

## Pipeline Design

### 1. Calibrate CFE with Qkrig

CFE parameters are calibrated per catchment using Qkrig kriged streamflow as
pseudo-observations. Calibration is shared across all four R-formula experiments —
the R formula is **not involved at calibration time**.

Script: `1_calibrate/calibrate_catchment_cfe_da_v2.py`  
Best params saved as: `{cat-id}_best_params.json`

| Flag | Default | Purpose |
|---|---|---|
| `--enkf-enabled` | off | Turn DA on during calibration run |
| `--enkf-members 20` | 20 | Ensemble size |
| `--no-vrugt-r` | off | Use raw kriging variance instead of Vrugt formula |
| `--vrugt-alpha 0.10` | 0.10 | α in Vrugt R formula |
| `--vrugt-scale 0.001` | 0.001 | Kriging variance scaling factor |

---

### 2. Assimilation Script for CFE

Script: `2_assimilation/run_perturbation_da_on.py`

#### 2a — Update Met Forcings (Forcing Arm, 30 members)

Stochastic perturbation of meteorological forcing to quantify uncertainty from
precipitation and PET inputs, with DA-corrected initial states:

- Precip: lognormal multiplier σ=0.15
- PET: Gaussian multiplier σ=0.10, clipped at 0

Produces `{cat-id}_da_forcing_arm.csv` — 30-member spaghetti over Helene window.  
**Figure 2a:** ensemble spread from met forcing uncertainty shows how sensitive
18-hr forecasts are to precipitation perturbations alone.

#### 2b — Update Hydro States using Kriging Error Variance in Vrugt R (Hydro-state Arm, 20 members)

Stochastic perturbation of initial hydrologic states (soil moisture, groundwater
storage) with deterministic forcing. DA uses Qkrig observations with Vrugt R:

```
R(t) = (α · y_obs(t))² + scale · σ²_krig(t)
       α = 0.10,  scale = 0.001
```

Produces `{cat-id}_da_hydro_arm.csv` — 20-member spaghetti over Helene window.  
**Figure 2b:** ensemble spread from initial state uncertainty shows how sensitive
18-hr forecasts are to hydrologic state perturbations alone.

#### 2c — 18-hour Forecast Cycle (save and restart)

At each initialization time t0:
1. EnKF analyses the ensemble state using Qkrig observation
2. Analyzed state snapshotted to `_da_snapshots.parquet` — enables restart
3. DA switched off — ensemble free-runs for 18 hours
4. At t0+1, DA resumes; a new 18-hour fork begins from the next snapshot

Script: `2_assimilation/run_lead_time_forecast_sweep.py`

#### 2d — Crossed Ensemble Design (600 members)

The forcing arm (2a) and hydro-state arm (2b) are crossed to form a 600-member ensemble
that captures both sources of uncertainty simultaneously:

```
member_k  →  forcing draw  k // 20   (selects 1 of 30 met perturbations)
          →  hydro draw    k % 20    (selects 1 of 20 initial-state draws)
```

Script: `2_assimilation/run_crossed_ensemble.py`

---

### 3. run_route.py with DA

Per-catchment runoff (mm/h) from each ensemble member is routed through the
channel network at every hour in the 18-hr forecast cycle for each initialization
timestep via T-route Muskingum-Cunge (`compute_network_structured`) to the
terminal reach wb-1016283 at USGS gauge 03463300.

Scripts: `3_routing/run_route.py`, `run_route_crossed_ensemble.py`,
`route_lead_time_forecasts.py`

GPKG: `/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg`

---

### 4. Evaluate at the Gauge

#### 4a — Forecast Error Decay: Convergence to Open Loop

How forecast skill degrades from lead hour 1 → 18. DA-initialized forecasts are
compared against the open loop (no DA) at each lead hour. Shows whether the
benefit of DA-corrected initial states persists or converges to open loop skill
by the end of the 18-hour window.

Scripts: `4_evaluation/4a_error_decay/`

#### 4b — Compare 600-member Ensemble vs Observations during Helene

Full 600-member crossed ensemble spread plotted against USGS obs at gauge 03463300
during Hurricane Helene. Shows whether the ensemble brackets the observed peak and
how the DA-corrected ensemble performs vs open loop.

Scripts: `4_evaluation/4b_ensemble_vs_obs/`

#### 4c — Ensemble Mean with Spread from All 18-hour Forecast Initializations

Reconstructed timeseries: for each valid_time during Helene, all forecasts that
land on that moment (across all 151 issue times × all lead hours × all members)
are pooled. Median + 5th/95th percentile envelope plotted against USGS obs.
Shows overall forecast skill as one continuous picture across the storm.

Scripts: `4_evaluation/4c_timeseries/`

---

## What This Repo Tests

Four experiments comparing how the **observation error variance R** is formulated
in the EnKF update. R controls how strongly the DA pulls model states toward the
kriging-interpolated observation (Qkrig) at each hourly timestep:

| Folder | R Formula | Key idea |
|---|---|---|
| [folder1_variance_scaled_vrugt/](folder1_variance_scaled_vrugt/) | `(0.10·y)² + 0.001·σ²_krig` | Flow-scaled R — uncertainty grows with flow magnitude (Vrugt 2005) |
| [folder2_fixed_r_007/](folder2_fixed_r_007/) | `0.07` (constant) | Simple fixed R — high constant Kalman gain regardless of flow |
| [folder3_dynamic_vrugt_seeded/](folder3_dynamic_vrugt_seeded/) | `(0.10·y)² + 0.001·σ²_krig` | Same Vrugt formula with dynamic kriging variance and fixed RNG seed |
| [folder4_dynamic_variance_direct/](folder4_dynamic_variance_direct/) | `σ²_krig` directly | Raw kriging variance as R — no flow-magnitude term |

---

## Key Results

| Folder | Full KGE | Full NSE | Helene KGE | Helene NSE | Helene peak |
|---|---|---|---|---|---|
| F1 — Vrugt | +0.277 | **+0.555** | +0.213 | **+0.459** | 645.7 m³/s (34%) |
| F2 — Fixed R=0.07 | **+0.503** | +0.163 | **+0.439** | -0.009 | **1227.1 m³/s (65%)** |
| F3 — Vrugt seeded | +0.261 | +0.485 | +0.189 | +0.372 | 693.8 m³/s (37%) |
| F4 — Direct σ²_krig | +0.200 | +0.428 | +0.132 | +0.303 | 667.9 m³/s (35%) |

**USGS Helene peak: 1885.7 m³/s**

**Finding:** Fixed R=0.07 (F2) gives the best KGE and peak capture. Constant small R
keeps Kalman gain K = P/(P+R) near 1 throughout the flood, so the model tracks
observations aggressively at every timestep. The Vrugt formula inflates R at high flows
(R ≈ 35,000 at the Helene peak), collapsing the gain exactly when the DA update matters
most. F1 achieves the best NSE because it fits the overall hydrograph shape better —
F2's high constant gain overshoots at low flows, hurting NSE, but captures the peak
far better.

---

## Which Comparisons Are Controlled

**Not all four folders can be directly compared to each other.** Two design differences
exist between F1/F2 and F3/F4:

### 1 — Different RNG seeds

| Folder | RNG seed |
|---|---|
| F1, F2 | `hash(catchment_id) & 0x7fffffff` — unique per catchment |
| F3, F4 | `42` — fixed, same for all catchments |

Different seeds → different perturbation realizations → different analysis trajectories,
even with the same R formula. The gap between F1 (+0.277) and F3 (+0.261) is explained
by seed difference alone, not methodology.

### 2 — Different kriging variance during Helene

F1/F2 use static σ²_krig ≈ 4.17 throughout.  
F3/F4 use dynamic σ²_krig that drops to ~1.6 at the Helene peak.

The **Qkrig flow values are identical** in both datasets; only the variance changes.
For F3 (Vrugt), the 0.001 weight on σ²_krig makes the variance difference negligible.
For F4 (direct σ²_krig), the lower variance at Helene peak means stronger DA updates —
but this is offset by the absence of the flow-magnitude term entirely.

### Valid controlled comparisons

| Pair | Seed | Obs variance | R formula | Valid? |
|---|---|---|---|---|
| F1 vs F2 | both hash | both static ~4.17 | Vrugt vs Fixed 0.07 | ✅ **clean** |
| F3 vs F4 | both seed=42 | both dynamic | Vrugt vs Direct σ² | ✅ **clean** |
| F1 vs F3 | different | small difference | same (Vrugt) | ⚠️ seed confound |
| F2 vs F4 | different | larger difference | different | ⚠️ seed + obs confound |

To fully compare all four R formulas in a controlled way, F3 and F4 would need to
be re-run with the hash-based seed and the same obs dataset as F1/F2. See
[EXPERIMENTS.md](EXPERIMENTS.md) for the full completeness matrix.

---

## Scientific Approach

### EnKF update

At each hourly timestep t:

```
K(t)  = P(t) / (P(t) + R(t))          Kalman gain
x̂(t)  = x(t) + K(t) · (y_obs(t) − H·x(t))   state update
```

where P is ensemble state variance, y_obs is Qkrig, and R is the observation error
variance. Each experiment differs only in how R is computed.

### Vrugt 2005 heteroscedastic R (F1, F3)

```
R(t) = (α · y_obs(t))² + scale · σ²_krig(t)
       α = 0.10   (10% relative error on the observation)
       scale = 0.001
```

R scales with flow magnitude: small at low flow (gain near 1, aggressive DA) and
large at peak flow. At Helene peak (~3.2 mm/h per catchment): R ≈ 0.1, K collapses.

### Fixed R (F2)

```
R = 0.07  (constant)
```

Kalman gain K = P/(P+0.07) stays near 1 throughout the simulation including the peak.
The DA update is equally aggressive at low and high flows.

### Direct kriging variance (F4)

```
R(t) = σ²_krig(t)
```

σ²_krig is the spatial kriging interpolation uncertainty — typically ~4.17 (static obs)
or ~1.6–2.5 at Helene peak (dynamic obs). This gives weaker DA than Fixed R=0.07 at
low flows and stronger DA at peak when σ²_krig drops.

---

## Folder Structure

Each experiment is self-contained with identical pipeline structure:

```
da_methods/
├── README.md                     ← this file
├── EXPERIMENTS.md                ← results table, completeness matrix, server paths
├── compare_all_folders.py        ← cross-experiment comparison plots
├── build_pptx.py                 ← builds helene_da_results.pptx
├── pptx_figures/                 ← all presentation figures (tracked)
│
├── folder1_variance_scaled_vrugt/    R = (0.10·y)² + 0.001·σ²_krig, hash seed
├── folder2_fixed_r_007/              R = 0.07 constant, hash seed
├── folder3_dynamic_vrugt_seeded/     R = (0.10·y)² + 0.001·σ²_krig, seed=42
└── folder4_dynamic_variance_direct/  R = σ²_krig directly, seed=42

Each folder/:
├── 1_calibrate/    calibrate_catchment_cfe_da_v2.py (shared params, R not involved)
├── 2_assimilation/ run_perturbation_da_on.py, run_crossed_ensemble.py,
│                   run_lead_time_forecast_sweep.py, batch scripts
├── 3_routing/      run_route.py, run_route_crossed_ensemble.py,
│                   route_lead_time_forecasts.py
└── 4_evaluation/
    ├── 4a_error_decay/      forecast skill decay vs lead time
    ├── 4b_ensemble_vs_obs/  600-member ensemble vs USGS
    └── 4c_timeseries/       reconstructed timeseries + spaghetti plots
```

---

## Running an Experiment (F3 as example)

```bash
# Stage A — perturbation arms (2a/2b) for all 21 catchments
bash folder3_dynamic_vrugt_seeded/2_assimilation/batch_run_all_f3.sh arms

# Stage B — 18hr forecast cycles (2c)
bash folder3_dynamic_vrugt_seeded/2_assimilation/batch_run_all_f3.sh forecast

# Stage C — 600-member crossed ensemble (2d)
bash folder3_dynamic_vrugt_seeded/2_assimilation/batch_run_all_f3.sh ensemble

# Route analysis trajectory + ensemble + forecast cycles (3_routing/)
# See folder README for exact run_route.py invocations

# Evaluate (4_evaluation/4c_timeseries/)
bash folder3_dynamic_vrugt_seeded/4_evaluation/4c_timeseries/run_4c_f3.sh
```

F4 works identically — replace `f3` with `f4` and use `batch_run_all_f4.sh`.

---

## Ensemble Traceability

Any of the 600 members can be fully reconstructed:

```python
member_k → forcing_draw_i = k // 20   # which of 30 met perturbation draws
         → hydro_draw_j   = k % 20    # which of 20 initial-state draws
```

Four provenance files are saved per catchment:

| File | Contents |
|---|---|
| `{cat-id}_da_snapshots.parquet` | DA-analyzed state at every issue time — enables restart |
| `{cat-id}_member_manifest.csv` | Decoder: member_k → (forcing_draw_i, hydro_draw_j) |
| `{cat-id}_hydro_draw_states.parquet` | Exact perturbed initial states for each of 20 hydro draws |
| `{cat-id}_forcing_draw_sequences.parquet` | Exact P/PET scale factors for each of 30 forcing draws |

---

## Server Paths (dualearth1)

All paths under `/mnt/disk2/` unless noted.

| Resource | Path |
|---|---|
| Retro forcing | `suma_helen_poster/nwm_retro_catchment_forcings/` |
| Test forcing 1 | `/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings/` |
| Test forcing 2 | `/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings/` |
| Obs (F1/F2, static σ²) | `1400_sites_helene/catchment_ts_03463300_with_variance/` |
| Obs (F3/F4, dynamic σ²) | `1400_sites_helene/catchment_ts_03463300_spliced_dyn_helene/` |
| BMI config | `suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json` |
| Param bounds | `suma_helen_poster/run_gpu/CFE_parameter_bounds.json` |
| CFE source | `suma_helen_poster/cfe_py/` |
| GPKG | `/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg` |
| USGS hourly obs | `suma_helen_poster/03463300_usgs_hourly_2018_2024.csv` |
| F1 analysis | `suma_helen_poster/da_results/v2_true_enkf_vrugt/` |
| F2 analysis | `suma_helen_poster/da_results/v2_fixed_r007_analysis/` |
| F3 analysis | `1400_sites_helene/da_results_dynamic_vrugt_seeded/` |
| F4 analysis | `1400_sites_helene/da_results_dynamic_novrugt_seeded/` |

See [EXPERIMENTS.md](EXPERIMENTS.md) for full server paths including forecast cycles,
ensemble directories, and routed output paths for all 4 folders.
