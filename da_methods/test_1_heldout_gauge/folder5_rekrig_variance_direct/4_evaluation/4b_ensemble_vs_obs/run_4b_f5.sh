#!/usr/bin/env bash
# F5 Re-kriging dynamic variance (1 gauge holdout) — 4b ensemble plots.
#
# Runs:
#   1. run_4b_crossed_f5.sh  — per-catchment crossed ensemble PNGs
#   2. plot_routed_ensemble_vs_usgs.py — F5 ensemble envelope vs USGS at outlet
#   3. plot_routed_ensemble_combined.py — F1 Vrugt vs F5 re-kriging comparison
#
# Must be run AFTER:
#   - crossed ensemble parquets exist (Kunal has done this)
#   - routed_crossed_ensemble.parquet exists (Kunal has done this)
#   - routed_Q_test.csv exists (run deterministic routing first)
#
# Usage:
#   bash run_4b_f5.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

F1_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt
F5_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

# 1. Per-catchment crossed ensemble plots
echo "[F5-4b] Per-catchment crossed ensemble plots..."
for CAT in "${CATS[@]}"; do
    PQ="$F5_DIR/$CAT/${CAT}_crossed_ensemble.parquet"
    if [ ! -f "$PQ" ]; then
        echo "  [$CAT] No parquet — skipping"
        continue
    fi
    echo "  [$CAT] crossed ensemble..."
    $TROUTE "$SCRIPT_DIR/plot_crossed_ensemble.py" \
        --ensemble-dir "$F5_DIR" \
        --cat-id       "$CAT"    \
        --usgs-csv     "$USGS_CSV"
done

# 2. Routed ensemble vs USGS (F5 only)
echo "[F5-4b] Routed ensemble vs USGS (F5 only)..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_vs_usgs.py" \
    --routed-pq "$F5_DIR/routed_crossed_ensemble.parquet" \
    --usgs-csv  "$USGS_CSV" \
    --label     "F5 Re-kriging variance — 1 gauge holdout" \
    --out-dir   "$F5_DIR"

# 3. Combined: F1 Vrugt vs F5 re-kriging
echo "[F5-4b] Combined comparison: F1 Vrugt vs F5 Re-kriging..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_combined.py" \
    --vrugt-csv   "$F1_DIR/routed_Q_test.csv" \
    --novrugt-csv "$F5_DIR/routed_Q_test.csv" \
    --ensemble-pq "$F5_DIR/routed_crossed_ensemble.parquet" \
    --out-dir     "$F5_DIR"

echo "[F5-4b] Done."
ls "$F5_DIR"/*.png 2>/dev/null || echo "  (no root PNGs yet)"
