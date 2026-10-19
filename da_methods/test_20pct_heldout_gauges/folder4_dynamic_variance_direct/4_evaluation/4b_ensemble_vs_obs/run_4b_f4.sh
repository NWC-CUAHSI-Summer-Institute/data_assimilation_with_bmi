#!/usr/bin/env bash
# F4 dynamic variance direct (20pct gauge holdout) — 4b routed ensemble plots.
#
# Runs:
#   1. plot_routed_ensemble_vs_usgs.py  — F4 ensemble envelope vs USGS at outlet
#   2. plot_routed_ensemble_combined.py — F1 Vrugt vs F4 dynamic variance direct comparison
#
# Requires routed_Q_test.csv (both folders) and routed_crossed_ensemble.parquet
# (CROSSED_DIR) to be present.
#
# Usage:
#   bash run_4b_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

F1_DA_DIR=/mnt/disk2/suma_helen_poster/da_results/v2_perturbation_da_on
F4_DA_DIR=/mnt/disk2/1400_sites_helene/da_arms_dynamic_novrugt_seeded
CROSSED_DIR=/mnt/disk2/1400_sites_helene/da_crossed_dynamic_novrugt_seeded
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OUT_DIR="$CROSSED_DIR"

echo "[F4-4b] Routed ensemble vs USGS (F4 only)..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_vs_usgs.py" \
    --routed-pq "$CROSSED_DIR/routed_crossed_ensemble.parquet" \
    --usgs-csv  "$USGS_CSV" \
    --label     "F4 Dynamic variance direct — 20pct gauge holdout" \
    --out-dir   "$OUT_DIR"

echo "[F4-4b] Combined comparison: F1 Vrugt vs F4 Dynamic variance direct..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_combined.py" \
    --vrugt-csv   "$F1_DA_DIR/routed_Q_test.csv" \
    --novrugt-csv "$F4_DA_DIR/routed_Q_test.csv" \
    --ensemble-pq "$CROSSED_DIR/routed_crossed_ensemble.parquet" \
    --out-dir     "$OUT_DIR"

echo "[F4-4b] Done. Outputs in: $OUT_DIR"
ls "$OUT_DIR"/*.png 2>/dev/null || echo "  (no PNGs in out dir)"
