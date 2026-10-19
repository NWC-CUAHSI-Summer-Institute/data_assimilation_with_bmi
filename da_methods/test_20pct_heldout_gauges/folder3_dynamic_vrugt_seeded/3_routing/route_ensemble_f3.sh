#!/usr/bin/env bash
# F3 dynamic Vrugt seeded (20pct gauge holdout) — crossed-ensemble T-route routing.
# Routes all 600 crossed-ensemble members through Muskingum-Cunge
# to gauge 03463300 and writes routed_crossed_ensemble.parquet.
#
# Must be run AFTER crossed ensemble parquets exist in CROSSED_DIR.
#
# Usage:
#   bash route_ensemble_f3.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route_crossed_ensemble.py

CROSSED_DIR=/mnt/disk2/1400_sites_helene/da_crossed_dynamic_vrugt_seeded
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[F3] Crossed-ensemble T-route routing..."
echo "  ensemble-dir: $CROSSED_DIR"
echo "  out-dir     : $CROSSED_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg         "$GPKG"         \
    --ensemble-dir "$CROSSED_DIR"  \
    --out-dir      "$CROSSED_DIR"  \
    --usgs-csv     "$USGS_CSV"

echo "[F3] Ensemble routing done. Output: $CROSSED_DIR/routed_crossed_ensemble.parquet"
