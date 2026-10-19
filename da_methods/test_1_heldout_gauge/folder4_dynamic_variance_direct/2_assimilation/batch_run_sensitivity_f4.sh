#!/usr/bin/env bash
# Batch perturbation sensitivity analysis — F4 Direct variance, 1 gauge holdout.
# Runs run_perturbation_sensitivity.py for all 21 catchments × 3 sources
# (init, forcing, process). DA is OFF; only perturbation source differs.
# Skips if output CSV already exists.
#
# Usage:
#   nohup bash batch_run_sensitivity_f4.sh > ~/logs/sensitivity_f4.log 2>&1 &

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT="$SCRIPT_DIR/run_perturbation_sensitivity.py"

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

SOURCES=(init forcing process)

OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

echo "F4 sensitivity analysis — out: $OUT_DIR"
echo ""

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    for SRC in "${SOURCES[@]}"; do
        OUT_FILE="$OUT_DIR/$CAT/${CAT}_sensitivity_${SRC}.csv"

        if [ -f "$OUT_FILE" ]; then
            echo "  [$CAT/$SRC] already exists — skipping"
            SKIP=$((SKIP + 1))
            continue
        fi

        echo "==============================="
        echo "=== $CAT  source=$SRC ==="
        echo "==============================="

        "$PYTHON" "$SCRIPT" \
            --cat-id            "$CAT"          \
            --source            "$SRC"          \
            --forcing-dir       "$FORCING_DIR1" \
            --obs-dir           "$OBS_DIR"      \
            --cfe-dir           "$CFE_DIR"      \
            --config-file       "$CONFIG_FILE"  \
            --param-bounds      "$PARAM_BOUNDS" \
            --out-dir           "$OUT_DIR"      \
            --test-forcing-dir1 "$FORCING_DIR1" \
            --test-forcing-dir2 "$FORCING_DIR2" \
        && DONE=$((DONE + 1)) \
        || { echo "  FAILED: $CAT/$SRC"; FAIL=$((FAIL + 1)); }
    done
done

echo ""
echo "=== F4 sensitivity done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
