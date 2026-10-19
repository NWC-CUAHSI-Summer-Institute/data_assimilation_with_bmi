#!/usr/bin/env bash
# Open-loop CFE run (no DA) — 20% gauge holdout, all 21 catchments.
# Runs calibrate_catchment_cfe_da_v2.py WITHOUT --enkf-enabled.
# Params staged from da_results_dynamic_novrugt_seeded (20% holdout calibration).
# Obs: catchment_ts_03463300_dynamic_variance (03463300 in Qkrig, 22 others withheld).
#
# Usage:
#   bash batch_run_openloop_20pct.sh > ~/logs/openloop_20pct.log 2>&1

set -euo pipefail

PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
PROD_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/1_distributed_cfe/calibrate_catchment_cfe_da_v2.py

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

PARAM_SRC=/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded
OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
FORCING_DIR=/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results/openloop_20pct_heldout

mkdir -p "$OUT_DIR"
echo "Open-loop 20pct (no DA) — out: $OUT_DIR"
echo ""

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    echo "==============================="
    echo "=== $CAT ==="
    echo "==============================="

    CAT_OUT="$OUT_DIR/$CAT"
    OUT_FILE="$CAT_OUT/${CAT}_test_results.csv"

    PARAMS="$PARAM_SRC/$CAT/${CAT}_best_params.json"
    if [ ! -f "$PARAMS" ]; then
        echo "  WARNING: No best_params for $CAT at $PARAMS — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    if [ -f "$OUT_FILE" ]; then
        echo "  Already exists — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    mkdir -p "$CAT_OUT"
    cp "$PARAMS" "$CAT_OUT/"

    "$PYTHON" "$PROD_SCRIPT" \
        --cat-id            "$CAT"          \
        --forcing-dir       "$FORCING_DIR"  \
        --obs-dir           "$OBS_DIR"      \
        --cfe-dir           "$CFE_DIR"      \
        --config-file       "$CONFIG_FILE"  \
        --param-bounds      "$PARAM_BOUNDS" \
        --out-dir           "$OUT_DIR"      \
        --test-forcing-dir1 "$FORCING_DIR1" \
        --test-forcing-dir2 "$FORCING_DIR2" \
        --enkf-members 1                    \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== Open-loop 20pct done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
