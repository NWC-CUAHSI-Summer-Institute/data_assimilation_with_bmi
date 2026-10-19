#!/usr/bin/env bash
# F5 re-kriged variance (1 gauge holdout) — 4a lead-time decay plots.
#
# Runs 3 plot scripts per catchment:
#   plot_lead_time_decay.py           -- pooled RMSE vs lead (catchment-level)
#   plot_lead_time_decay_by_regime.py -- same split by flow regime
#   plot_lead_time_decay_gauge.py     -- gauge-level (requires routed parquets)
#
# F5 has no _test_results.csv per catchment, so obs are read from the
# re-kriged obs directory (--krig-obs-dir).
#
# Usage:
#   bash run_4a_f5.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

F5_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct
KRIG_OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance_rekrig
ROUTE_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "[F5-4a] Lead-time decay plots — leadtime dir: $F5_DIR"

for CAT in "${CATS[@]}"; do
    DA_CSV="$F5_DIR/$CAT/${CAT}_lead_time_forecasts_da.csv"
    if [ ! -f "$DA_CSV" ]; then
        echo "  [$CAT] No lead-time CSV — skipping"
        continue
    fi

    echo "  [$CAT] Plotting pooled decay..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay.py" \
        --cat-id        "$CAT"           \
        --leadtime-dir  "$F5_DIR"        \
        --da-dir        "$F5_DIR"        \
        --krig-obs-dir  "$KRIG_OBS_DIR"

    echo "  [$CAT] Plotting decay by regime..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay_by_regime.py" \
        --cat-id        "$CAT"           \
        --leadtime-dir  "$F5_DIR"        \
        --da-dir        "$F5_DIR"        \
        --krig-obs-dir  "$KRIG_OBS_DIR"
done

if [ -f "$ROUTE_DIR/routed_leadtime_da_full.parquet" ]; then
    echo "[F5-4a] Gauge-level decay (routed)..."
    $TROUTE "$SCRIPT_DIR/plot_lead_time_decay_gauge.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"

    echo "[F5-4a] Fixed-target error plots..."
    $TROUTE "$SCRIPT_DIR/plot_forecast_error_fixed_target.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"

    echo "[F5-4a] Per-init error plots..."
    $TROUTE "$SCRIPT_DIR/plot_forecast_error_per_init.py" \
        --route-dir "$ROUTE_DIR" \
        --usgs-csv  "$USGS_CSV"
else
    echo "[F5-4a] Routed parquets not found — skipping gauge-level and error plots."
    echo "         Run route_leadtime_f5.sh (in 3_routing/) first, then re-run this script."
fi

echo "[F5-4a] Done."
