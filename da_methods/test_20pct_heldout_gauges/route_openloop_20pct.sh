#!/usr/bin/env bash
# Route open-loop 20% holdout outputs through T-route to gauge 03463300.
# Run AFTER batch_run_openloop_20pct.sh completes.
#
# Usage:
#   bash route_openloop_20pct.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route.py

OL_DIR=/mnt/disk2/suma_helen_poster/da_results/openloop_20pct_heldout
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[OL-20pct] T-route routing..."
echo "  da-dir : $OL_DIR"
echo "  out    : $OL_DIR/routed_Q_test.csv"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg     "$GPKG"    \
    --da-dir   "$OL_DIR"  \
    --out-dir  "$OL_DIR"  \
    --usgs-csv "$USGS_CSV"

echo "[OL-20pct] Done."
