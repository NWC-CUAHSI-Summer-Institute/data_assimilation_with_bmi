#!/usr/bin/env bash
# Regenerate 2b hydro-state arm plots for F1, F2, F4, F5 with the routed
# open-loop baseline overlaid on each panel.
#
# Uses openloop/routed_Q_test.csv (Q_routed_m3s column, already in m³/s,
# correct timing via T-route) — same line on every catchment plot.
#
# Overwrites existing _2b_hydro_arm_helene.png files in each catchment subfolder.
#
# Usage:
#   bash run_hydro_arm_plots_openloop.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/test_1_heldout_gauge/folder1_variance_scaled_vrugt/4_evaluation/plot_da_perturbation_arms.py
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OL_CSV=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/openloop/routed_Q_test.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

declare -A EXPS
EXPS["F1"]=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt
EXPS["F2"]=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007
EXPS["F4"]=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct
EXPS["F5"]=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct

if [ ! -f "$OL_CSV" ]; then
    echo "ERROR: routed open-loop not found: $OL_CSV"
    exit 1
fi
echo "Open-loop: $OL_CSV"

SKIP=0; DONE=0; FAIL=0

for EXP in F1 F2 F4 F5; do
    ARM_DIR="${EXPS[$EXP]}"
    echo ""
    echo "======================================="
    echo "=== $EXP  arm-dir: $ARM_DIR ==="
    echo "======================================="

    for CAT in "${CATS[@]}"; do
        HYDRO_CSV="$ARM_DIR/$CAT/${CAT}_da_hydro_arm.csv"

        if [ ! -f "$HYDRO_CSV" ]; then
            echo "  [$EXP $CAT] No hydro arm CSV — skipping"
            SKIP=$((SKIP + 1)); continue
        fi

        echo "  [$EXP $CAT] plotting..."
        $TROUTE "$SCRIPT" \
            --arm-dir   "$ARM_DIR"  \
            --cat-id    "$CAT"      \
            --usgs-csv  "$USGS_CSV" \
            --ol-ts-csv "$OL_CSV"   \
        && DONE=$((DONE + 1)) \
        || { echo "  FAILED: $EXP $CAT"; FAIL=$((FAIL + 1)); }
    done
done

echo ""
echo "=== Hydro arm plots (with routed open loop) done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
