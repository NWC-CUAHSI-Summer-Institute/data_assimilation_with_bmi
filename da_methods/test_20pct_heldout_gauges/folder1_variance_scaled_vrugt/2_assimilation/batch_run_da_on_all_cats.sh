#!/bin/bash
# batch_run_da_on_all_cats.sh
# Runs run_perturbation_da_on.py for all 21 catchments.
# Stages best_params.json from v2_true_enkf_pn for each catchment.
# cat-1016300 is skipped (already done).

set -e

CATS=(
  cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
  cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
  cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
  cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results/v2_perturbation_da_on
BEST_PARAMS_SRC=/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_pn
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

SKIP=0
DONE=0
FAIL=0

for CAT in "${CATS[@]}"; do
    echo ""
    echo "==============================="
    echo "=== $CAT ==="
    echo "==============================="

    CAT_OUT=$OUT_DIR/$CAT
    mkdir -p "$CAT_OUT"

    # Stage best_params if not already present
    if [ ! -f "$CAT_OUT/${CAT}_best_params.json" ]; then
        SRC="$BEST_PARAMS_SRC/$CAT/${CAT}_best_params.json"
        if [ -f "$SRC" ]; then
            cp "$SRC" "$CAT_OUT/"
            echo "  Staged best_params from v2_true_enkf_pn"
        else
            echo "  WARNING: No best_params for $CAT — skipping"
            SKIP=$((SKIP + 1))
            continue
        fi
    else
        echo "  best_params already staged"
    fi

    # Skip if both arm CSVs already exist
    if [ -f "$CAT_OUT/${CAT}_da_forcing_arm.csv" ] && \
       [ -f "$CAT_OUT/${CAT}_da_hydro_arm.csv" ]; then
        echo "  Arms already complete — skipping"
        SKIP=$((SKIP + 1))
        continue
    fi

    python3 ~/run_perturbation_da_on.py \
        --cat-id "$CAT" \
        --forcing-dir "$FORCING_DIR1" \
        --obs-dir "$OBS_DIR" \
        --cfe-dir "$CFE_DIR" \
        --config-file "$CONFIG_FILE" \
        --param-bounds "$PARAM_BOUNDS" \
        --out-dir "$OUT_DIR" \
        --test-forcing-dir1 "$FORCING_DIR1" \
        --test-forcing-dir2 "$FORCING_DIR2" \
        --hardcoded-r 0.07 \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== Batch complete: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
