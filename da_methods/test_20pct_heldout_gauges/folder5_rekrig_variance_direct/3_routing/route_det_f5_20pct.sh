#!/usr/bin/env bash
# F5 re-kriged variance (20% gauge holdout) — deterministic T-route routing.
# Routes per-catchment _test_results.csv through Muskingum-Cunge to gauge 03463300.
#
# Run AFTER batch_run_f5_20pct.sh completes.
#
# Usage:
#   bash route_det_f5_20pct.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route.py

F5_DIR=/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[F5-20pct] Deterministic T-route routing..."
echo "  da-dir : $F5_DIR"
echo "  out-dir: $F5_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg     "$GPKG"     \
    --da-dir   "$F5_DIR"   \
    --out-dir  "$F5_DIR"   \
    --usgs-csv "$USGS_CSV"

echo "[F5-20pct] Done. Output: $F5_DIR/routed_Q_test.csv"
