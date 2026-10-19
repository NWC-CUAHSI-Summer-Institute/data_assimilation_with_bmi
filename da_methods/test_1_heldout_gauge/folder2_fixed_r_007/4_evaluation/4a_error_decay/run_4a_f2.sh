#!/usr/bin/env bash
# F2 Fixed R=0.07 (1 gauge holdout) — 4a lead-time decay plots.
#
# Runs 3 plot scripts per catchment:
#   plot_lead_time_decay.py          — pooled RMSE vs lead (catchment-level)
#   plot_lead_time_decay_by_regime.py — same split by flow regime
#   plot_lead_time_decay_gauge.py    — gauge-level (requires routed parquets)
#
# plot_lead_time_decay_gauge.py is skipped if the routed parquets don't exist.
#
# Usage:
#   bash run_4a_f2.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

LEADTIME_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007
DA_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007
OBS_DIR=/home/svyas/catchment_ts_no_03463300_gapfilled
ROUTE_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder2_fixed_r007_leadtime_routed
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "[F2-4a] Lead-time decay plots — leadtime dir: $LEADTIME_DIR"

for CAT in "${CATS[@]}"; do
    DA_CSV="$LEADTIME_DIR/$CAT/${CAT}_lead_time_forecasts_da.csv"
    if [ ! -f "$DA_CSV" ]; then
        echo "  [$CAT] No lead-time CSV — skipping"
        continue
    fi

    echo "  [$CAT] Plotting pooled decay..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay.py" \
        --cat-id       "$CAT"          \
        --leadtime-dir "$LEADTIME_DIR" \
        --da-dir       "$DA_DIR"       \
        --obs-dir      "$OBS_DIR"

    echo "  [$CAT] Plotting decay by regime..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay_by_regime.py" \
        --cat-id       "$CAT"          \
        --leadtime-dir "$LEADTIME_DIR" \
        --da-dir       "$DA_DIR"       \
        --obs-dir      "$OBS_DIR"
done

# Gauge-level decay (needs routed parquets)
if [ -f "$ROUTE_DIR/routed_leadtime_da_full.parquet" ]; then
    echo "[F2-4a] Gauge-level decay (routed)..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay_gauge.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"

    echo "[F2-4a] Fixed-target error plots..."
    $TROUTE "$SCRIPT_DIR/plot_forecast_error_fixed_target.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"

    echo "[F2-4a] Per-init error plots..."
    $TROUTE "$SCRIPT_DIR/plot_forecast_error_per_init.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"
else
    echo "[F2-4a] Routed parquets not found — skipping gauge-level and error plots."
    echo "         Run route_leadtime_f2.sh first, then re-run this script."
fi

echo "[F2-4a] Done."
