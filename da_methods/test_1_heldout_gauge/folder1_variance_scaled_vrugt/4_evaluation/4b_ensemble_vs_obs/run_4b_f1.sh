#!/usr/bin/env bash
# F1 Vrugt (1 gauge holdout) — 4b ensemble vs obs plots.
#
# Per-catchment: perturbation_category_shaded, per_member_factor_decomp,
#   production_ensemble_forecast.
# Multi-catchment (runs once): sensitivity_spaghetti, sensitivity_spread.
#
# Usage:
#   bash run_4b_f1.sh

set -euo pipefail

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EVAL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "[F1-4b] Per-catchment ensemble plots — data: $DATA_DIR"

for CAT in "${CATS[@]}"; do
    PROD_CSV="$DATA_DIR/$CAT/${CAT}_production_per_member.csv"
    SEN_CSV="$DATA_DIR/$CAT/${CAT}_sensitivity_init.csv"

    if [ ! -f "$PROD_CSV" ] || [ ! -f "$SEN_CSV" ]; then
        echo "  [$CAT] Missing prod or sensitivity CSV — skipping"
        continue
    fi

    echo "  [$CAT] perturbation category shaded..."
    $TROUTE "$SCRIPT_DIR/plot_perturbation_category_shaded.py"  --cat-id "$CAT"

    echo "  [$CAT] per-member factor decomposition..."
    $TROUTE "$SCRIPT_DIR/plot_per_member_factor_decomposition.py" --cat-id "$CAT"

    echo "  [$CAT] production ensemble forecast..."
    $TROUTE "$SCRIPT_DIR/plot_production_ensemble_forecast.py"  --cat-id "$CAT"
done

echo "[F1-4b] Multi-catchment sensitivity plots..."
$TROUTE "$EVAL_DIR/plot_perturbation_sensitivity_spaghetti.py"
$TROUTE "$EVAL_DIR/plot_perturbation_sensitivity_spread.py"

echo "[F1-4b] Done."
