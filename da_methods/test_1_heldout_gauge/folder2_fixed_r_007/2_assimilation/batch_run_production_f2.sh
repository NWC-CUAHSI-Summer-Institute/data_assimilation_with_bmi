#!/usr/bin/env bash
# Batch production per-member run — F2 Fixed R=0.07, 1 gauge holdout.
# Runs run_production_per_member.py for all 21 catchments with --hardcoded-r 0.07.
# Skips if output CSV already exists.
#
# Usage:
#   nohup bash batch_run_production_f2.sh > ~/logs/production_f2.log 2>&1 &

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT="$SCRIPT_DIR/run_production_per_member.py"
PROD_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/1_distributed_cfe/calibrate_catchment_cfe_da_v2.py

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

OBS_DIR=/home/svyas/catchment_ts_no_03463300_gapfilled
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

echo "F2 production per-member (R=0.07 fixed) — out: $OUT_DIR"
echo ""

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    echo "==============================="
    echo "=== $CAT ==="
    echo "==============================="

    CAT_OUT="$OUT_DIR/$CAT"
    OUT_FILE="$CAT_OUT/${CAT}_production_per_member.csv"

    if [ ! -f "$CAT_OUT/${CAT}_best_params.json" ]; then
        echo "  WARNING: No best_params for $CAT — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    if [ -f "$OUT_FILE" ]; then
        echo "  Production per-member CSV already exists — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    "$PYTHON" "$SCRIPT" \
        --cat-id            "$CAT"          \
        --forcing-dir       "$FORCING_DIR1" \
        --obs-dir           "$OBS_DIR"      \
        --cfe-dir           "$CFE_DIR"      \
        --config-file       "$CONFIG_FILE"  \
        --param-bounds      "$PARAM_BOUNDS" \
        --out-dir           "$OUT_DIR"      \
        --test-forcing-dir1 "$FORCING_DIR1" \
        --test-forcing-dir2 "$FORCING_DIR2" \
        --enkf-members 20                   \
        --hardcoded-r  0.07                 \
        --prod-script  "$PROD_SCRIPT"       \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== F2 production per-member done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
