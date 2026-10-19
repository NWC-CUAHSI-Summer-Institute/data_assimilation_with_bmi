#!/usr/bin/env bash
# F4 — Dynamic variance direct: 4c reconstructed timeseries plots
#
# Reads routed_leadtime_da_full.parquet + routed_leadtime_openloop_full.parquet
# from the F4 forecast route dir and writes two PNGs there.
#
# Usage:
#   bash run_4c_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

ROUTE_DIR="/mnt/disk2/suma_helen_poster/da_results/da_forecast_dynamic_novrugt_seeded_routed"
USGS_CSV="/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

echo "[F4-4c] Reconstructed timeseries from: $ROUTE_DIR"
$TROUTE "$SCRIPT_DIR/plot_lead_time_reconstructed_timeseries.py" \
    --route-dir "$ROUTE_DIR" \
    --usgs-csv  "$USGS_CSV"

echo "[F4-4c] Done. Outputs:"
ls "$ROUTE_DIR"/lead_time_reconstructed_timeseries*.png 2>/dev/null || echo "  (no PNGs found)"
