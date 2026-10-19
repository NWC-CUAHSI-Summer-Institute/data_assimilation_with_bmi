#!/usr/bin/env bash
# F2 fixed R=0.07 (20pct gauge holdout) — deterministic T-route routing.
# Routes per-catchment _test_results.csv (sim_mm_h column) through Muskingum-Cunge
# to gauge 03463300 and writes routed_Q_test.csv + KGE/NSE summary.
#
# Usage:
#   bash route_det_f2.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route.py

F2_DA_DIR=/mnt/disk2/suma_helen_poster/da_results/v2_fixed_r007_analysis
F2_OUT_DIR=/mnt/disk2/suma_helen_poster/da_results/fixed_r007_routed
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

mkdir -p "$F2_OUT_DIR"

echo "[F2] Deterministic T-route routing..."
echo "  da-dir : $F2_DA_DIR"
echo "  out-dir: $F2_OUT_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg    "$GPKG"        \
    --da-dir  "$F2_DA_DIR"  \
    --out-dir "$F2_OUT_DIR" \
    --usgs-csv "$USGS_CSV"

echo "[F2] Deterministic routing done. Output: $F2_OUT_DIR/routed_Q_test.csv"
