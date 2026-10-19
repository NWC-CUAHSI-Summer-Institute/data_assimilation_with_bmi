#!/usr/bin/env bash
# F5 re-kriged variance (20% holdout) — gauge-level lead-time NSE evaluation.
# Run AFTER the lead-time forecast parquets have been routed via T-Route.
#
# Usage:
#   bash run_4a_f5.sh

set -euo pipefail

PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ROUTE_DIR=/mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct/routed
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct/figures

mkdir -p "$OUT_DIR"

echo "[F5-20pct] Lead-time NSE evaluation..."
$PYTHON "$SCRIPT_DIR/plot_lead_time_nse_gauge_f5.py" \
    --route-dir "$ROUTE_DIR" \
    --usgs-csv  "$USGS_CSV"  \
    --out-dir   "$OUT_DIR"

echo "[F5-20pct] Done. Figures in: $OUT_DIR"
