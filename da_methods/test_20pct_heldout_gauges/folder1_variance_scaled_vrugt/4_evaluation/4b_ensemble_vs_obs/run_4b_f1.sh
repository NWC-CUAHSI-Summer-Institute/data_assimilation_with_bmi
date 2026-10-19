#!/usr/bin/env bash
# F1 variance-scaled Vrugt (20pct gauge holdout) — 4b routed ensemble plots.
#
# Runs:
#   1. plot_routed_ensemble_vs_usgs.py  — F1 ensemble envelope vs USGS at outlet
#
# Requires routed_Q_test.csv and routed_crossed_ensemble.parquet in CROSSED_DIR.
#
# Usage:
#   bash run_4b_f1.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

CROSSED_DIR=/mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble_vrugt
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OUT_DIR="$CROSSED_DIR"

echo "[F1-4b] Routed ensemble vs USGS (F1 only)..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_vs_usgs.py" \
    --routed-pq "$CROSSED_DIR/routed_crossed_ensemble.parquet" \
    --usgs-csv  "$USGS_CSV" \
    --label     "F1 Variance-scaled Vrugt — 20pct gauge holdout" \
    --out-dir   "$OUT_DIR"

echo "[F1-4b] Done. Outputs in: $OUT_DIR"
ls "$OUT_DIR"/*.png 2>/dev/null || echo "  (no PNGs in out dir)"
