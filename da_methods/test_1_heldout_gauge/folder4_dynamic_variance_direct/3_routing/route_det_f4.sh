#!/usr/bin/env bash
# F4 dynamic variance direct (1 gauge holdout) — deterministic T-route routing.
# Routes per-catchment _test_results.csv (sim_mm_h column) through Muskingum-Cunge
# to gauge 03463300 and writes routed_Q_test.csv + KGE/NSE summary.
#
# Usage:
#   bash route_det_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
ROUTE_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/2_troute_routing/run_route.py

F4_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder4_dynamic_variance_direct
GPKG=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv
KV_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance

echo "[F4] Deterministic T-route routing..."
echo "  da-dir : $F4_DIR"
echo "  out-dir: $F4_DIR"
echo "  kv-dir : $KV_DIR"

$TROUTE "$ROUTE_SCRIPT" \
    --gpkg    "$GPKG"    \
    --da-dir  "$F4_DIR"  \
    --out-dir "$F4_DIR"  \
    --usgs-csv "$USGS_CSV" \
    --kv-dir  "$KV_DIR"

echo "[F4] Deterministic routing done. Output: $F4_DIR/routed_Q_test.csv"
