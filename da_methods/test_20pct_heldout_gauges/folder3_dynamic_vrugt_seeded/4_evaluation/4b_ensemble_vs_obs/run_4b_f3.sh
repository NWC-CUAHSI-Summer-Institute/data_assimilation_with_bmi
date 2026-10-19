#!/usr/bin/env bash
# F3 dynamic Vrugt seeded (20pct gauge holdout) — 4b routed ensemble plots.
#
# Runs:
#   1. plot_routed_ensemble_vs_usgs.py  — F3 ensemble envelope vs USGS at outlet
#   2. plot_routed_ensemble_combined.py — F1 Vrugt vs F3 dynamic Vrugt seeded comparison
#
# Requires routed_Q_test.csv (both folders) and routed_crossed_ensemble.parquet
# (CROSSED_DIR) to be present.
#
# Usage:
#   bash run_4b_f3.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

F1_DA_DIR=/mnt/disk2/suma_helen_poster/da_results/v2_perturbation_da_on
F3_DA_DIR=/mnt/disk2/1400_sites_helene/da_arms_dynamic_vrugt_seeded
CROSSED_DIR=/mnt/disk2/1400_sites_helene/da_crossed_dynamic_vrugt_seeded
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OUT_DIR="$CROSSED_DIR"

echo "[F3-4b] Routed ensemble vs USGS (F3 only)..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_vs_usgs.py" \
    --routed-pq "$CROSSED_DIR/routed_crossed_ensemble.parquet" \
    --usgs-csv  "$USGS_CSV" \
    --label     "F3 Dynamic Vrugt seeded — 20pct gauge holdout" \
    --out-dir   "$OUT_DIR"

echo "[F3-4b] Combined comparison: F1 Vrugt vs F3 Dynamic Vrugt seeded..."
$TROUTE "$SCRIPT_DIR/plot_routed_ensemble_combined.py" \
    --vrugt-csv   "$F1_DA_DIR/routed_Q_test.csv" \
    --novrugt-csv "$F3_DA_DIR/routed_Q_test.csv" \
    --ensemble-pq "$CROSSED_DIR/routed_crossed_ensemble.parquet" \
    --out-dir     "$OUT_DIR"

echo "[F3-4b] Done. Outputs in: $OUT_DIR"
ls "$OUT_DIR"/*.png 2>/dev/null || echo "  (no PNGs in out dir)"
