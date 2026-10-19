#!/usr/bin/env bash
# Folder 4 — Dynamic variance direct: full 2_assimilation pipeline
#
# Runs 3 stages for all 21 catchments:
#   Stage A: 2a/2b perturbation arms  (run_perturbation_da_on.py, R=krig_var, seed=42)
#   Stage B: 2c 18hr forecast cycles  (run_lead_time_forecast_sweep.py, --no-vrugt-r, seed=42)
#   Stage C: 2d 600-member ensemble   (run_crossed_ensemble.py, --direct-variance, seed=42)
#
# Usage:
#   bash batch_run_all_f4.sh            # runs all 3 stages
#   bash batch_run_all_f4.sh arms       # stage A only
#   bash batch_run_all_f4.sh forecast   # stage B only
#   bash batch_run_all_f4.sh ensemble   # stage C only

set -euo pipefail

STAGE="${1:-all}"

SCRIPT_DIR="$(dirname "$0")"
TROUTE=/home/svyas/miniconda3/envs/troute/bin/python

# Paths
DA_SRC="/mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded"
OBS_DIR="/mnt/disk2/1400_sites_helene/catchment_ts_03463300_spliced_dyn_helene"
CFE_DIR="/mnt/disk2/suma_helen_poster/cfe_py"
CONFIG_FILE="/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
PARAM_BOUNDS="/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json"
FORCING_DIR="/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings"
TEST_FORCING_DIR1="/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings"
TEST_FORCING_DIR2="/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings"

OUT_ARMS="/mnt/disk2/1400_sites_helene/da_arms_dynamic_novrugt_seeded"
OUT_FORECAST="/mnt/disk2/1400_sites_helene/da_forecast_dynamic_novrugt_seeded"
OUT_ENSEMBLE="/mnt/disk2/1400_sites_helene/da_crossed_dynamic_novrugt_seeded"

LOG_DIR="$HOME/f4_logs"
mkdir -p "$LOG_DIR"

CATS=$(ls -d "$DA_SRC"/cat-* 2>/dev/null | xargs -I{} basename {})

# ── Stage A: 2a/2b perturbation arms ─────────────────────────────────────────
run_arms() {
  echo "[F4] Stage A: perturbation arms (R=krig_var direct, seed=42)"
  for cat_id in $CATS; do
    params="$DA_SRC/$cat_id/${cat_id}_best_params.json"
    [ -f "$params" ] || { echo "  [skip] $cat_id — no params"; continue; }
    mkdir -p "$OUT_ARMS/$cat_id"
    cp "$params" "$OUT_ARMS/$cat_id/"
    log="$LOG_DIR/arms_${cat_id}.log"
    echo "  [run] $cat_id"
    $TROUTE "$SCRIPT_DIR/run_perturbation_da_on.py" \
      --cat-id       "$cat_id" \
      --forcing-dir  "$FORCING_DIR" \
      --obs-dir      "$OBS_DIR" \
      --cfe-dir      "$CFE_DIR" \
      --config-file  "$CONFIG_FILE" \
      --param-bounds "$PARAM_BOUNDS" \
      --out-dir      "$OUT_ARMS" \
      --test-forcing-dir1 "$TEST_FORCING_DIR1" \
      --test-forcing-dir2 "$TEST_FORCING_DIR2" \
      --direct-variance \
      --rng-seed 42 \
      > "$log" 2>&1 && echo "  [ok ] $cat_id" || echo "  [FAIL] $cat_id — see $log"
  done
}

# ── Stage B: 2c 18hr forecast cycles ─────────────────────────────────────────
run_forecast() {
  echo "[F4] Stage B: 18hr forecast cycles (R=krig_var direct, seed=42)"
  for cat_id in $CATS; do
    params="$DA_SRC/$cat_id/${cat_id}_best_params.json"
    [ -f "$params" ] || { echo "  [skip] $cat_id — no params"; continue; }
    mkdir -p "$OUT_FORECAST/$cat_id"
    cp "$params" "$OUT_FORECAST/$cat_id/"
    log="$LOG_DIR/forecast_${cat_id}.log"
    echo "  [run] $cat_id"
    $TROUTE "$SCRIPT_DIR/run_lead_time_forecast_sweep.py" \
      --cat-id       "$cat_id" \
      --forcing-dir  "$FORCING_DIR" \
      --obs-dir      "$OBS_DIR" \
      --cfe-dir      "$CFE_DIR" \
      --config-file  "$CONFIG_FILE" \
      --param-bounds "$PARAM_BOUNDS" \
      --out-dir      "$OUT_FORECAST" \
      --test-forcing-dir1 "$TEST_FORCING_DIR1" \
      --test-forcing-dir2 "$TEST_FORCING_DIR2" \
      --no-vrugt-r \
      --rng-seed 42 \
      > "$log" 2>&1 && echo "  [ok ] $cat_id" || echo "  [FAIL] $cat_id — see $log"
  done
}

# ── Stage C: 2d 600-member crossed ensemble ───────────────────────────────────
run_ensemble() {
  echo "[F4] Stage C: 600-member crossed ensemble (R=krig_var direct, seed=42)"
  for cat_id in $CATS; do
    params="$DA_SRC/$cat_id/${cat_id}_best_params.json"
    [ -f "$params" ] || { echo "  [skip] $cat_id — no params"; continue; }
    mkdir -p "$OUT_ENSEMBLE/$cat_id"
    cp "$params" "$OUT_ENSEMBLE/$cat_id/"
    log="$LOG_DIR/ensemble_${cat_id}.log"
    echo "  [run] $cat_id"
    $TROUTE "$SCRIPT_DIR/run_crossed_ensemble.py" \
      --cat-id       "$cat_id" \
      --forcing-dir  "$FORCING_DIR" \
      --obs-dir      "$OBS_DIR" \
      --cfe-dir      "$CFE_DIR" \
      --config-file  "$CONFIG_FILE" \
      --out-dir      "$OUT_ENSEMBLE" \
      --test-forcing-dir1 "$TEST_FORCING_DIR1" \
      --test-forcing-dir2 "$TEST_FORCING_DIR2" \
      --direct-variance \
      --rng-seed 42 \
      > "$log" 2>&1 && echo "  [ok ] $cat_id" || echo "  [FAIL] $cat_id — see $log"
  done
}

case "$STAGE" in
  arms)     run_arms ;;
  forecast) run_forecast ;;
  ensemble) run_ensemble ;;
  all)      run_arms; run_forecast; run_ensemble ;;
  *) echo "Usage: $0 [arms|forecast|ensemble|all]"; exit 1 ;;
esac

echo "[F4] Done."
