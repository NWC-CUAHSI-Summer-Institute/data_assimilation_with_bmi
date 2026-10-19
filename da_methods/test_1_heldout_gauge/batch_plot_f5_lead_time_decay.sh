#!/usr/bin/env bash
# Plot lead-time error decay curves for F5 — all 21 catchments.
# Run AFTER batch_run_f5_lead_time_sweep.sh has completed.
#
# Usage:
#   bash batch_plot_f5_lead_time_decay.sh

set -euo pipefail

PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/test_1_heldout_gauge/folder5_rekrig_variance_direct/4_evaluation/4a_error_decay/plot_lead_time_decay.py
F5_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    DA_CSV="$F5_DIR/$CAT/${CAT}_lead_time_forecasts_da.csv"
    if [ ! -f "$DA_CSV" ]; then
        echo "  [$CAT] No lead-time CSV — skipping (run sweep first)"
        SKIP=$((SKIP + 1)); continue
    fi

    echo "  [$CAT] plotting..."
    "$PYTHON" "$SCRIPT" \
        --cat-id       "$CAT"   \
        --leadtime-dir "$F5_DIR" \
        --da-dir       "$F5_DIR" \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== F5 lead-time decay plots done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
