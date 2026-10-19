#!/usr/bin/env bash
# F1 Vrugt (1 gauge holdout) — crossed-ensemble T-route routing.
# Routes all 600 crossed-ensemble members through Muskingum-Cunge
# to gauge 03463300 and writes routed_crossed_ensemble.parquet.
#
# Must be run AFTER crossed ensemble parquets exist in F1_DIR.
#
# Usage:
#   bash route_ensemble_f1.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route_crossed_ensemble.py

F1_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[F1] Crossed-ensemble T-route routing..."
echo "  ensemble-dir: $F1_DIR"
echo "  out-dir     : $F1_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg         "$GPKG"    \
    --ensemble-dir "$F1_DIR"  \
    --out-dir      "$F1_DIR"  \
    --usgs-csv     "$USGS_CSV"

echo "[F1] Ensemble routing done. Output: $F1_DIR/routed_crossed_ensemble.parquet"
