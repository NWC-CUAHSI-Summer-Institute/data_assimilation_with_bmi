#!/usr/bin/env bash
# F4 dynamic variance direct (20pct gauge holdout) — 4b crossed ensemble plots.
# Runs plot_crossed_ensemble.py for all 21 catchments.
# Reads <DATA_DIR>/<cat>/<cat>_crossed_ensemble.parquet (already generated).
#
# Usage:
#   bash run_4b_crossed_f4.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

DATA_DIR=/mnt/disk2/1400_sites_helene/da_crossed_dynamic_novrugt_seeded
USGS_CSV=/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "[F4-4b-crossed] Crossed ensemble plots — data: $DATA_DIR"

for CAT in "${CATS[@]}"; do
    PQ="$DATA_DIR/$CAT/${CAT}_crossed_ensemble.parquet"
    if [ ! -f "$PQ" ]; then
        echo "  [$CAT] No crossed ensemble parquet — skipping"
        continue
    fi
    echo "  [$CAT] crossed ensemble..."
    $TROUTE "$SCRIPT_DIR/plot_crossed_ensemble.py" \
        --ensemble-dir "$DATA_DIR" \
        --cat-id       "$CAT"      \
        --usgs-csv     "$USGS_CSV"
done

echo "[F4-4b-crossed] Done."
