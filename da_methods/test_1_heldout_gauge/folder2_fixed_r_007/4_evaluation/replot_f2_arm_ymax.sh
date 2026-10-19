#!/usr/bin/env bash
# Regenerate F2 2b hydro-state arm plots with y-axis capped at 2200 m³/s
# and the routed open-loop baseline overlaid.
#
# Usage:
#   bash replot_f2_arm_ymax.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/test_1_heldout_gauge/folder1_variance_scaled_vrugt/4_evaluation/plot_da_perturbation_arms.py
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OL_CSV=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/openloop/routed_Q_test.csv
F2_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

if [ ! -f "$OL_CSV" ]; then
    echo "ERROR: routed open-loop not found: $OL_CSV"
    exit 1
fi

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    HYDRO_CSV="$F2_DIR/$CAT/${CAT}_da_hydro_arm.csv"

    if [ ! -f "$HYDRO_CSV" ]; then
        echo "  [$CAT] No hydro arm CSV — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    echo "  [$CAT] plotting (ymax=2200)..."
    $TROUTE "$SCRIPT" \
        --arm-dir   "$F2_DIR"   \
        --cat-id    "$CAT"      \
        --usgs-csv  "$USGS_CSV" \
        --ol-ts-csv "$OL_CSV"   \
        --ymax      2200        \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== F2 replot done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
