#!/usr/bin/env bash
# F5 re-kriged variance (20% holdout) — route lead-time forecast CSVs.
# Routes 21-catchment lead-time forecast CSVs through Muskingum-Cunge
# to gauge 03463300. Writes routed_leadtime_da_full.parquet and
# routed_leadtime_openloop_full.parquet for gauge-level NSE decay analysis.
#
# Run AFTER batch_run_f5_analysis_20pct.sh forecast completes.
#
# Usage:
#   bash route_leadtime_f5.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route_leadtime_forecasts.py

F5_DIR=/mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
OUT_DIR=$F5_DIR/routed

echo "[F5-20pct] Routing lead-time forecasts through T-Route..."
echo "  forecast-dir: $F5_DIR"
echo "  out-dir     : $OUT_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg         "$GPKG"   \
    --forecast-dir "$F5_DIR" \
    --out-dir      "$OUT_DIR"

echo "[F5-20pct] Lead-time routing done. Outputs in: $OUT_DIR"
ls "$OUT_DIR"/*.parquet 2>/dev/null || echo "  (no parquets found)"
