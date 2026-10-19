#!/usr/bin/env bash
# F4 — Dynamic variance direct: 4c reconstructed timeseries plots
#
# Reads routed_leadtime_da_full.parquet + routed_leadtime_openloop_full.parquet
# from the F4 forecast route dir and writes two PNGs there.
#
# Usage:
#   bash run_4c_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

LEADTIME_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct
ROUTE_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct_leadtime_routed
OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "[F4-4c] Reconstructed timeseries from: $ROUTE_DIR"
$TROUTE "$SCRIPT_DIR/plot_lead_time_reconstructed_timeseries.py" \
    --route-dir "$ROUTE_DIR" \
    --usgs-csv  "$USGS_CSV"

echo "[F4-4c] Forecast spaghetti..."
$TROUTE "$SCRIPT_DIR/plot_forecast_spaghetti.py" \
    --route-dir "$ROUTE_DIR" \
    --usgs-csv  "$USGS_CSV"

echo "[F4-4c] Helene vs low-flow issue-time hydrograph (per catchment)..."
for CAT in "${CATS[@]}"; do
    DA_CSV="$LEADTIME_DIR/$CAT/${CAT}_lead_time_forecasts_da.csv"
    [ -f "$DA_CSV" ] || { echo "  [$CAT] no leadtime CSV — skipping"; continue; }
    echo "  [$CAT] helene hydrograph..."
    $TROUTE "$SCRIPT_DIR/plot_helene_issue_time_hydrograph.py" \
        --cat-id       "$CAT"          \
        --leadtime-dir "$LEADTIME_DIR" \
        --obs-dir      "$OBS_DIR"
done

echo "[F4-4c] Done. Outputs in: $ROUTE_DIR"
ls "$ROUTE_DIR"/*.png 2>/dev/null || echo "  (no PNGs in route dir)"
