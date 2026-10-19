#!/usr/bin/env bash
# Batch: run CFE+EnKF analysis (R=0.07) for all 21 catchments → v2_fixed_r007_analysis/
# Discovers catchments from the existing hardcoded-R lead-time forecast dir.
# After this completes, route with:
#   conda run -n troute python3 da_methods/3_routing/run_route.py \
#       --gpkg /mnt/disk2/suma_helen_poster/gauge_03463300_network.gpkg \
#       --da-dir /mnt/disk2/suma_helen_poster/da_results/v2_fixed_r007_analysis \
#       --out-dir /mnt/disk2/suma_helen_poster/da_results/fixed_r007_routed \
#       --usgs-csv /mnt/disk2/suma_helen_poster/da_results/usgs_03463300_test.csv \
#       --kv-dir /mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance

set -euo pipefail

SCRIPT="$(dirname "$0")/run_analysis_fixed_r.py"
HARDCODED_R_DIR="/mnt/disk2/suma_helen_poster/da_results/v2_lead_time_forecast_hardcoded_r"
OBS_DIR="/mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance"
CFE_DIR="/mnt/disk2/suma_helen_poster/cfe_py"
CONFIG_FILE="/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
OUT_DIR="/mnt/disk2/suma_helen_poster/da_results/v2_fixed_r007_analysis"
LOG_DIR="$HOME/fixed_r007_analysis_logs"

mkdir -p "$LOG_DIR"

done=0; skipped=0; failed=0

for cat_dir in "$HARDCODED_R_DIR"/cat-*; do
    cat_id=$(basename "$cat_dir")
    forcing_file="$cat_dir/${cat_id}_nwm_operational_combined.csv"
    params_file="$cat_dir/${cat_id}_best_params.json"
    out_csv="$OUT_DIR/$cat_id/${cat_id}_test_results.csv"

    if [ -f "$out_csv" ]; then
        echo "[skip] $cat_id — already done"
        skipped=$((skipped + 1))
        continue
    fi

    if [ ! -f "$forcing_file" ]; then
        echo "[warn] $cat_id — forcing not found: $forcing_file"
        failed=$((failed + 1))
        continue
    fi

    if [ ! -f "$params_file" ]; then
        echo "[warn] $cat_id — params not found: $params_file"
        failed=$((failed + 1))
        continue
    fi

    log_file="$LOG_DIR/${cat_id}.log"
    echo "[run ] $cat_id"
    python3 "$SCRIPT" \
        --cat-id       "$cat_id" \
        --forcing-file "$forcing_file" \
        --obs-dir      "$OBS_DIR" \
        --params-file  "$params_file" \
        --cfe-dir      "$CFE_DIR" \
        --config-file  "$CONFIG_FILE" \
        --out-dir      "$OUT_DIR" \
        --hardcoded-r  0.07 \
        > "$log_file" 2>&1 \
    && { echo "[ok  ] $cat_id"; done=$((done + 1)); } \
    || { echo "[FAIL] $cat_id — see $log_file"; failed=$((failed + 1)); }
done

echo ""
echo "Summary: done=$done  skipped=$skipped  failed=$failed"
