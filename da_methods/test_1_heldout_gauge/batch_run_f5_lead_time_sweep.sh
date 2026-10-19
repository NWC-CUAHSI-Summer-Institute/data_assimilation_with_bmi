#!/usr/bin/env bash
# Lead-time forecast sweep for F5 (re-kriged variance direct) — all 21 catchments.
# Runs run_lead_time_forecast_sweep.py with --direct-variance and --no-vrugt-r,
# generating _lead_time_forecasts_da.csv and _lead_time_forecasts_openloop.csv
# per catchment. Plot with batch_plot_f5_lead_time_decay.sh after this completes.
#
# Runtime: ~30-60 min per catchment (8760-hour test period × 2 trajectories).
#
# Usage:
#   nohup bash batch_run_f5_lead_time_sweep.sh > ~/logs/f5_lead_time_sweep.log 2>&1 &

set -euo pipefail

PYTHON=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/test_1_heldout_gauge/folder5_rekrig_variance_direct/2_assimilation/run_lead_time_forecast_sweep.py
PROD_SCRIPT=/home/svyas/cfe_DA_bmi/da_methods/1_distributed_cfe/calibrate_catchment_cfe_da_v2.py

OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance_rekrig
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
OUT_DIR=/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct
FORCING_DIR=/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings
FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

CATS=(
    cat-1016279 cat-1016280 cat-1016281 cat-1016282 cat-1016283
    cat-1016300
    cat-1016301 cat-1016302 cat-1016303 cat-1016304 cat-1016305
    cat-1016306 cat-1016307 cat-1016308 cat-1016309 cat-1016310
    cat-1016311 cat-1016312 cat-1016313 cat-1016314 cat-1016315
)

echo "F5 lead-time sweep — out: $OUT_DIR"
echo ""

SKIP=0; DONE=0; FAIL=0

for CAT in "${CATS[@]}"; do
    echo "==============================="
    echo "=== $CAT ==="
    echo "==============================="

    DA_CSV="$OUT_DIR/$CAT/${CAT}_lead_time_forecasts_da.csv"
    OL_CSV="$OUT_DIR/$CAT/${CAT}_lead_time_forecasts_openloop.csv"

    if [ -f "$DA_CSV" ] && [ -f "$OL_CSV" ]; then
        echo "  Already exists — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    if [ ! -f "$OUT_DIR/$CAT/${CAT}_best_params.json" ]; then
        echo "  WARNING: No best_params for $CAT — skipping"
        SKIP=$((SKIP + 1)); continue
    fi

    "$PYTHON" "$SCRIPT" \
        --cat-id            "$CAT"          \
        --forcing-dir       "$FORCING_DIR"  \
        --obs-dir           "$OBS_DIR"      \
        --cfe-dir           "$CFE_DIR"      \
        --config-file       "$CONFIG_FILE"  \
        --param-bounds      "$PARAM_BOUNDS" \
        --out-dir           "$OUT_DIR"      \
        --test-forcing-dir1 "$FORCING_DIR1" \
        --test-forcing-dir2 "$FORCING_DIR2" \
        --enkf-members      20              \
        --no-vrugt-r                        \
        --direct-variance                   \
        --prod-script       "$PROD_SCRIPT"  \
    && DONE=$((DONE + 1)) \
    || { echo "  FAILED: $CAT"; FAIL=$((FAIL + 1)); }
done

echo ""
echo "=== F5 lead-time sweep done: done=$DONE  skipped=$SKIP  failed=$FAIL ==="
