#!/bin/bash
# Alabama storm 2026 — catchment-level CFE calibration (2020–2022)
# 12 gauges: 02361000 02369800 02371500 02372250 02374500 02408540
#            02422500 02450250 02464000 02465493 02469800 03574500
#
# Run:
#   chmod +x ~/batch_calibrate_al_storm.sh
#   nohup ~/batch_calibrate_al_storm.sh > ~/logs_cal_al_storm/run.log 2>&1 &
set -euo pipefail

JOBS=$(nproc)
FORCING_DIR="/mnt/disk2/svyas_clean_run/nwm_retro_al_storm_2026"
OBS_DIR="/mnt/disk2/al_storm_2026/catchment_obs_r100"
CFE_DIR="/mnt/disk2/suma_helen_poster/cfe_py"
CONFIG_FILE="/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
PARAM_BOUNDS="/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json"
OUT_DIR="/mnt/disk2/svyas_clean_run/calibration_al_storm_2026"
LOG_DIR="$HOME/logs_cal_al_storm"
PYTHON="/home/svyas/miniconda3/envs/troute/bin/python3"
SCRIPT="$HOME/calibrate_catchment_heldin.py"
N_ITER=1000

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "=== Alabama storm 2026 — CFE calibration (2020–2022) ==="
echo "Cores: $JOBS  |  N iterations: $N_ITER"
echo "FORCING: $FORCING_DIR"
echo "OBS:     $OBS_DIR"
echo "OUT:     $OUT_DIR"

for f in "$CFE_DIR" "$CONFIG_FILE" "$PARAM_BOUNDS" "$FORCING_DIR" "$OBS_DIR"; do
    [[ -e "$f" ]] || { echo "ERROR: missing: $f"; exit 1; }
done
echo "All prerequisites found."

ALL_CATS=()
while IFS= read -r f; do
    ALL_CATS+=("$(basename "$f" .csv)")
done < <(ls "$FORCING_DIR"/*.csv)
TOTAL=${#ALL_CATS[@]}
echo "Total catchments: $TOTAL"

RUNNING=0; DONE=0; SKIPPED=0
pids=()

for cat in "${ALL_CATS[@]}"; do
    best_file="$OUT_DIR/$cat/${cat}_best_params.json"
    if [[ -f "$best_file" ]]; then (( SKIPPED++ )) || true; continue; fi
    if [[ ! -f "$FORCING_DIR/$cat.csv" ]]; then echo "  SKIP $cat — no forcing"; (( SKIPPED++ )) || true; continue; fi
    if [[ ! -f "$OBS_DIR/$cat.csv" ]];    then echo "  SKIP $cat — no obs";     (( SKIPPED++ )) || true; continue; fi

    while true; do
        live=()
        for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && live+=("$pid")
        done
        pids=("${live[@]}")
        (( ${#pids[@]} < JOBS )) && break
        sleep 2
    done

    (
        $PYTHON "$SCRIPT" \
            --cat-id       "$cat" \
            --forcing-dir  "$FORCING_DIR" \
            --obs-dir      "$OBS_DIR" \
            --cfe-dir      "$CFE_DIR" \
            --config-file  "$CONFIG_FILE" \
            --param-bounds "$PARAM_BOUNDS" \
            --out-dir      "$OUT_DIR" \
            --N            "$N_ITER" \
        > "$LOG_DIR/${cat}.log" 2>&1
        echo "  DONE: $cat"
    ) &
    pids+=($!)
    (( RUNNING++ )) || true; (( DONE++ )) || true
    (( DONE % 50 == 0 )) && echo "Launched $DONE / $TOTAL (skipped $SKIPPED)..."
done

wait
echo ""
echo "=== All Alabama calibrations complete ==="
echo "  Launched: $DONE  |  Skipped: $SKIPPED  |  Logs: $LOG_DIR"
