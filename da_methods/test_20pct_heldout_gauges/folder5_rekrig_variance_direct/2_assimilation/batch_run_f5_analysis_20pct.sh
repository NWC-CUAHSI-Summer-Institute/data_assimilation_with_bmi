#!/usr/bin/env bash
# F5 re-kriged variance (20% gauge holdout) — analysis pipeline.
#
# Stage A: perturbation arms  (run_perturbation_da_on.py, --direct-variance, seed=42)
# Stage B: lead time forecast (run_lead_time_forecast_sweep.py, --no-vrugt-r, seed=42)
#
# Run AFTER batch_run_f5_20pct.sh (production DA) completes.
#
# Usage:
#   bash batch_run_f5_analysis_20pct.sh arms
#   bash batch_run_f5_analysis_20pct.sh forecast
#   bash batch_run_f5_analysis_20pct.sh all

set -euo pipefail

STAGE="${1:-all}"

TROUTE=/home/svyas/miniconda3/envs/troute/bin/python3
SCRIPT_DIR=/home/svyas

DA_SRC=/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded
OBS_DIR=/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig
CFE_DIR=/mnt/disk2/suma_helen_poster/cfe_py
CONFIG_FILE=/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json
PARAM_BOUNDS=/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json
FORCING_DIR=/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings
TEST_FORCING_DIR1=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings
TEST_FORCING_DIR2=/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings

OUT_ARMS=/mnt/disk2/suma_helen_poster/da_results/da_arms_f5_rekrig_20pct
OUT_FORECAST=/mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct

LOG_DIR=$HOME/f5_20pct_logs
mkdir -p "$LOG_DIR"

CATS=$(ls -d "$DA_SRC"/cat-* 2>/dev/null | xargs -I{} basename {})

# ── Stage A: perturbation arms ────────────────────────────────────────────────
run_arms() {
    echo "[F5-20pct] Stage A: perturbation arms (re-kriged variance direct, seed=42)"
    for cat_id in $CATS; do
        params="$DA_SRC/$cat_id/${cat_id}_best_params.json"
        [ -f "$params" ] || { echo "  [skip] $cat_id — no params"; continue; }
        mkdir -p "$OUT_ARMS/$cat_id"
        cp "$params" "$OUT_ARMS/$cat_id/"
        log="$LOG_DIR/arms_${cat_id}.log"
        echo "  [run] $cat_id"
        $TROUTE "$SCRIPT_DIR/run_perturbation_da_on.py" \
            --cat-id            "$cat_id"           \
            --forcing-dir       "$FORCING_DIR"       \
            --obs-dir           "$OBS_DIR"           \
            --cfe-dir           "$CFE_DIR"           \
            --config-file       "$CONFIG_FILE"       \
            --param-bounds      "$PARAM_BOUNDS"      \
            --out-dir           "$OUT_ARMS"          \
            --test-forcing-dir1 "$TEST_FORCING_DIR1" \
            --test-forcing-dir2 "$TEST_FORCING_DIR2" \
            --direct-variance                        \
            --rng-seed 42                            \
            > "$log" 2>&1 && echo "  [ok ] $cat_id" || echo "  [FAIL] $cat_id — see $log"
    done
}

# ── Stage B: lead time forecast sweep ────────────────────────────────────────
run_forecast() {
    echo "[F5-20pct] Stage B: lead time forecast sweep (re-kriged variance, seed=42)"
    for cat_id in $CATS; do
        params="$DA_SRC/$cat_id/${cat_id}_best_params.json"
        [ -f "$params" ] || { echo "  [skip] $cat_id — no params"; continue; }
        mkdir -p "$OUT_FORECAST/$cat_id"
        cp "$params" "$OUT_FORECAST/$cat_id/"
        log="$LOG_DIR/forecast_${cat_id}.log"
        echo "  [run] $cat_id"
        $TROUTE "$SCRIPT_DIR/run_lead_time_forecast_sweep.py" \
            --cat-id            "$cat_id"           \
            --forcing-dir       "$FORCING_DIR"       \
            --obs-dir           "$OBS_DIR"           \
            --cfe-dir           "$CFE_DIR"           \
            --config-file       "$CONFIG_FILE"       \
            --param-bounds      "$PARAM_BOUNDS"      \
            --out-dir           "$OUT_FORECAST"      \
            --test-forcing-dir1 "$TEST_FORCING_DIR1" \
            --test-forcing-dir2 "$TEST_FORCING_DIR2" \
            --no-vrugt-r                             \
            --rng-seed 42                            \
            > "$log" 2>&1 && echo "  [ok ] $cat_id" || echo "  [FAIL] $cat_id — see $log"
    done
}

case "$STAGE" in
    arms)     run_arms ;;
    forecast) run_forecast ;;
    all)      run_arms; run_forecast ;;
    *) echo "Usage: $0 [arms|forecast|all]"; exit 1 ;;
esac

echo "[F5-20pct] Analysis done."
