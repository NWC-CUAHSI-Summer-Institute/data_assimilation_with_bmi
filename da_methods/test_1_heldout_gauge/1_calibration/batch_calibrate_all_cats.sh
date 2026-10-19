#!/usr/bin/env bash
# Batch DDS calibration for all 21 catchments — 1 held-out gauge Qkrig obs.
# Runs one background job per catchment; logs go to $OUT_DIR/<cat-id>_cal.log.
# After all jobs finish, check logs for any non-zero exit.
#
# Usage (from server home dir, troute env active):
#   bash batch_calibrate_all_cats.sh

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
CALIBRATE=/mnt/disk2/suma_helen_poster/calibrate_catchment_nwm.py
FORCING_DIR=/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings
OBS_DIR=$HOME/catchment_ts_no_03463300_gapfilled
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
OUT_DIR=/mnt/disk2/suma_helen_poster/catchment_results_1gauge_heldout
TEST_FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
TEST_FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings
N=1000

mkdir -p "$OUT_DIR"
echo "Starting calibration for ${#CATS[@]} catchments"
echo "  obs-dir : $OBS_DIR"
echo "  out-dir : $OUT_DIR"
echo ""

for CAT in "${CATS[@]}"; do
    LOG="$OUT_DIR/${CAT}_cal.log"
    echo "  Launching $CAT  →  $LOG"
    "$PYTHON" "$CALIBRATE" \
        --cat-id            "$CAT"              \
        --forcing-dir       "$FORCING_DIR"      \
        --obs-dir           "$OBS_DIR"          \
        --cfe-dir           "$CFE_DIR"          \
        --config-file       "$CONFIG_FILE"      \
        --param-bounds      "$PARAM_BOUNDS"     \
        --out-dir           "$OUT_DIR"          \
        --test-forcing-dir1 "$TEST_FORCING_DIR1" \
        --test-forcing-dir2 "$TEST_FORCING_DIR2" \
        --N "$N"                                \
        > "$LOG" 2>&1 &
done

echo ""
echo "All ${#CATS[@]} jobs launched. Waiting for completion..."
wait
echo ""
echo "All catchments done. Check for failures:"
grep -l "Error\|Traceback\|error" "$OUT_DIR"/cat-*_cal.log 2>/dev/null || echo "  No errors found in logs."
