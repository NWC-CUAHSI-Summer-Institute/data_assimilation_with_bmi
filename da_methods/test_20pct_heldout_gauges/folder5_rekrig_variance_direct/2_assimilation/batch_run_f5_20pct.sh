#!/usr/bin/env bash
# F5 re-kriged variance (20% gauge holdout) — DA assimilation for all 21 catchments.
#
# R(t) = sigma2_rekrig  (re-kriged per-hour variogram variance, --no-vrugt-r)
# Obs : /mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig/
# Params sourced from F4 DA results (shared calibration).
#
# Usage:
#   bash batch_run_f5_20pct.sh > ~/logs/f5_20pct.log 2>&1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
PROD_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/1_distributed_cfe/calibrate_catchment_cfe_da_v2.py

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

# Re-kriged obs (per-hour variogram refit, 03463300 withheld)
OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig

# Params source — F4 results carry best_params for all 21 catchments
PARAM_SRC=/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded

CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
FORCING_DIR=/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

OUT_DIR=/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct

mkdir -p "$OUT_DIR"
echo "[F5-20pct] Re-kriged variance DA — out: $OUT_DIR"
echo ""

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    echo "==============================="
    echo "=== $CAT ==="
    echo "==============================="

    CAT_OUT="$OUT_DIR/$CAT"
    OUT_FILE="$CAT_OUT/${CAT}_test_results.csv"

    # Stage best_params from F4 param source
    PARAMS="$PARAM_SRC/$CAT/${CAT}_best_params.json"
    if [ ! -f "$PARAMS" ]; then
        echo "  WARNING: No best_params for $CAT at $PARAMS — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    if [ -f "$OUT_FILE" ]; then
        echo "  Already exists — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    mkdir -p "$CAT_OUT"
    cp "$PARAMS" "$CAT_OUT/"

    "$PYTHON" "$PROD_SCRIPT" \
        --cat-id            "$CAT"          \
        --forcing-dir       "$FORCING_DIR"  \
        --obs-dir           "$OBS_DIR"      \
        --cfe-dir           "$CFE_DIR"      \
        --config-file       "$CONFIG_FILE"  \
        --param-bounds      "$PARAM_BOUNDS" \
        --out-dir           "$OUT_DIR"      \
        --test-forcing-dir1 "$FORCING_DIR1" \
        --test-forcing-dir2 "$FORCING_DIR2" \
        --enkf-enabled                      \
        --enkf-members 20                   \
        --no-vrugt-r                        \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== F5-20pct done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
