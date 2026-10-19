#!/usr/bin/env bash
# F3 dynamic Vrugt seeded (20pct gauge holdout) — deterministic T-route routing.
# Routes per-catchment _test_results.csv (sim_mm_h column) through Muskingum-Cunge
# to gauge 03463300 and writes routed_Q_test.csv + KGE/NSE summary.
#
# Usage:
#   bash route_det_f3.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route.py

F3_DIR=/mnt/disk2/1400_sites_helene/da_arms_dynamic_vrugt_seeded
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

echo "[F3] Deterministic T-route routing..."
echo "  da-dir : $F3_DIR"
echo "  out-dir: $F3_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg    "$GPKG"    \
    --da-dir  "$F3_DIR"  \
    --out-dir "$F3_DIR"  \
    --usgs-csv "$USGS_CSV"

echo "[F3] Deterministic routing done. Output: $F3_DIR/routed_Q_test.csv"
