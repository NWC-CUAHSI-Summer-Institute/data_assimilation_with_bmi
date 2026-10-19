#!/usr/bin/env bash
# Generate hydro-state arm plots (2b_hydro_arm_helene.png) for F1, F2, F4, F5
# of the 1-gauge holdout experiment — all 21 catchments per experiment.
#
# Reads _da_hydro_arm.csv and _da_forcing_arm.csv from each experiment's dir.
# Saves PNGs alongside the arm CSVs in each catchment subfolder.
#
# Usage:
#   bash run_hydro_arm_plots.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/test_1_heldout_gauge/folder1_variance_scaled_vrugt/4_evaluation/plot_da_perturbation_arms.py
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

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
            continue
        fi
        PNG="$ARM_DIR/$CAT/${CAT}_2b_hydro_arm_helene.png"
        if [ -f "$PNG" ]; then
            echo "  [$EXP $CAT] Already exists — skipping"
            continue
        fi
        echo "  [$EXP $CAT] plotting..."
        $TROUTE "$SCRIPT" \
            --arm-dir  "$ARM_DIR" \
            --cat-id   "$CAT"     \
            --usgs-csv "$USGS_CSV"
    done
done

echo ""
echo "=== Hydro arm plots done ==="
