#!/usr/bin/env bash
# F4 Direct variance (1 gauge holdout) — crossed-ensemble T-route routing.
# Routes all 600 crossed-ensemble members through Muskingum-Cunge
# to gauge 03463300 and writes routed_crossed_ensemble.parquet.
#
# Must be run AFTER crossed ensemble parquets exist in F4_DIR.
#
# Usage:
#   bash route_ensemble_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route_crossed_ensemble.py

F4_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[F4] Crossed-ensemble T-route routing..."
echo "  ensemble-dir: $F4_DIR"
echo "  out-dir     : $F4_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg         "$GPKG"    \
    --ensemble-dir "$F4_DIR"  \
    --out-dir      "$F4_DIR"  \
    --usgs-csv     "$USGS_CSV"

echo "[F4] Ensemble routing done. Output: $F4_DIR/routed_crossed_ensemble.parquet"
