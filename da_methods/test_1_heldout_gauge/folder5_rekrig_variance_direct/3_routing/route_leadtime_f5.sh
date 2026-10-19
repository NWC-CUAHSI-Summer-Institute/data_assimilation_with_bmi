#!/usr/bin/env bash
# F5 re-kriging variance (1 gauge holdout) — route lead-time forecast CSVs.
# Routes all 21-catchment lead-time forecast CSVs through Muskingum-Cunge
# to gauge 03463300. Writes routed_leadtime_da_full.parquet and
# routed_leadtime_openloop_full.parquet for gauge-level error decay analysis.
#
# Run AFTER batch_run_leadtime_f5.sh finishes.
#
# Usage:
#   bash route_leadtime_f5.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route_leadtime_forecasts.py

F5_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed

echo "[F5] Routing lead-time forecasts through T-route..."
echo "  forecast-dir: $F5_DIR"
echo "  out-dir     : $OUT_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg         "$GPKG"   \
    --forecast-dir "$F5_DIR" \
    --out-dir      "$OUT_DIR"

echo "[F5] Lead-time routing done. Outputs in: $OUT_DIR"
ls "$OUT_DIR"/*.parquet 2>/dev/null || echo "  (no parquets found)"
