# Test: 1 Held-Out Gauge (03463300)

## What this experiment tests

Same four DA R-formula experiments as `test_20pct_heldout_gauges/`, but with a
fundamentally better Qkrig observation field:

| | 20pct holdout | **This experiment** |
|---|---|---|
| Gauges withheld from kriging | ~20% of all South Toe gauges | Only gauge 03463300 (the outlet we evaluate) |
| Qkrig quality near outlet | Degraded — many local gauges missing | Much better — only 1 gauge missing |
| Expected DA performance | Baseline | **Should improve** — EnKF gets better observations |

The hypothesis: with only one gauge held out, the kriging field is far more
accurate near the South Toe outlet. The EnKF should therefore correct CFE states
more reliably, improving routed streamflow at 03463300.

---

## Steps to run (in order)

### Step 1 — Get new Qkrig observations (PENDING DATA)
New per-catchment Qkrig time series with only gauge 03463300 held out.

```
# TODO: update when data is available
NEW_OBS_DIR = /mnt/disk2/???/catchment_ts_1gauge_heldout_with_variance/
```

### Step 2 — Re-calibrate CFE (uses shared script)
Run `1_distributed_cfe/calibrate_catchment_cfe_da_v2.py` for all 21 catchments
using the new obs dir above.

```bash
python ../../1_distributed_cfe/calibrate_catchment_cfe_da_v2.py \
    --cat-id cat-XXXXXXX \
    --obs-dir <NEW_OBS_DIR> \           # <-- changes from 20pct version
    --forcing-dir <NWM_RETRO_DIR> \
    --cfe-dir /mnt/disk2/suma_helen_poster/cfe_py \
    --config-file <NEW_CONFIG_FILE> \   # <-- new per-catchment config from re-cal
    --param-bounds /mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json \
    --N 1000
```

Outputs per catchment: `best_params.json`, `cat-XXXXXXX_cal_results.csv`

### Step 3 — T-route routing (uses shared script)
Same as 20pct version. Run `2_troute_routing/run_route.py` using the new
calibration outputs from Step 2.

### Step 4 — Run DA experiments (F1–F4)
Scripts are copied from `test_20pct_heldout_gauges/`. Before running, update
these paths in each folder's `2_assimilation/run_perturbation_da_on.py`:

```python
# Lines to update in run_perturbation_da_on.py for each of F1–F4:
--obs-dir      <NEW_OBS_DIR>               # 1-gauge holdout Qkrig
--config-file  <NEW_PER_CATCHMENT_CONFIG>  # from re-calibration
--out-dir      <NEW_OUTPUT_DIR>            # separate from 20pct results
```

Static and dynamic variance paths for F1/F3/F4 will also need updating
when the new calibration results are available.

### Step 5 — Evaluation
Same evaluation scripts as `test_20pct_heldout_gauges/`. Compare results
against 20pct to quantify the improvement from using better Qkrig observations.

---

## PENDING — waiting on

- [ ] New Qkrig obs dir path (1 gauge held out)
- [ ] New per-catchment config files from re-calibration
- [ ] New output directory paths on server
- [ ] Static/dynamic variance paths for F3/F4

Once these are provided, update paths in each `2_assimilation/run_perturbation_da_on.py`
and the corresponding batch shell scripts.
