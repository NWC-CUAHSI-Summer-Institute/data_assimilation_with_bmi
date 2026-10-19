"""
Launch hardcoded-R (R=0.07) lead-time sweep for all 21 catchments.

- Discovers catchment IDs from the existing Vrugt-R output dir
- Pre-stages best_params.json from run3 into the hardcoded-R out dir
- Skips catchments that already completed
- Launches the sweep with bounded parallelism (each catchment runs ~60 BMI
  instances, so 4-5 concurrent is usually safe on this box)
- Each sweep's stdout/stderr goes to ~/leadtime_logs/sweep_<cat>.log

Run:
    python3 run_sweep_all_21_hardcoded_r.py --concurrent 5
    nohup python3 run_sweep_all_21_hardcoded_r.py --concurrent 5 \
        > ~/sweep_all_21_master.log 2>&1 &
    # then monitor:
    tail -f ~/sweep_all_21_master.log
    ls ~/leadtime_logs/

Each catchment ~6-12 hr. With concurrency=5, total wall clock ~2-4 batches
× 6-12 hr each ≈ overnight + a few hours.
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

EXISTING_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_lead_time_forecast"
OUT_DIR_ROOT = "/mnt/disk2/suma_helen_poster/da_results/v2_lead_time_forecast_hardcoded_r"
RUN3_DIR     = "/mnt/disk2/suma_helen_poster/catchment_results_range100_run3"
SWEEP_SCRIPT = os.path.expanduser("~/run_lead_time_forecast_sweep.py")
LOG_DIR      = os.path.expanduser("~/leadtime_logs")

COMMON_ARGS = [
    "--forcing-dir",  "/mnt/disk2/suma_helen_poster/nwm_retro_catchment_forcings",
    "--obs-dir",      "/mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance",
    "--cfe-dir",      "/mnt/disk2/suma_helen_poster/cfe_py",
    "--config-file",  "/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json",
    "--param-bounds", "/mnt/disk2/suma_helen_poster/run_gpu/CFE_parameter_bounds.json",
    "--out-dir",      OUT_DIR_ROOT,
    "--test-forcing-dir1",
        "/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings",
    "--test-forcing-dir2",
        "/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings",
    "--enkf-members", "20",
    "--hardcoded-r",  "0.07",
]


def stage_params(cat_id):
    out_cat_dir = Path(OUT_DIR_ROOT) / cat_id
    out_cat_dir.mkdir(parents=True, exist_ok=True)
    src = Path(RUN3_DIR) / cat_id / f"{cat_id}_best_params.json"
    dst = out_cat_dir / f"{cat_id}_best_params.json"
    if dst.exists():
        return True
    if not src.exists():
        print(f"  [warn] no run3 params for {cat_id} at {src}")
        return False
    dst.write_text(src.read_text())
    return True


def run_one(cat_id):
    if not stage_params(cat_id):
        return cat_id, -1, "missing params"
    log_path = Path(LOG_DIR) / f"sweep_{cat_id}.log"
    cmd = ["python3", SWEEP_SCRIPT, "--cat-id", cat_id, *COMMON_ARGS]
    with open(log_path, "w") as logf:
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
    return cat_id, result.returncode, str(log_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--concurrent', type=int, default=5,
                        help='Max concurrent sweeps. Each uses ~60 BMI instances + '
                             'Python overhead, so 5 ≈ 300 BMIs in flight. '
                             'Drop to 3 if memory is tight.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print plan and exit without launching.')
    args = parser.parse_args()

    os.makedirs(LOG_DIR, exist_ok=True)

    # Discover catchments
    catchments = sorted([
        d.name for d in Path(EXISTING_DIR).iterdir()
        if d.is_dir() and d.name.startswith('cat-')
    ])
    print(f"[plan] {len(catchments)} catchments discovered under {EXISTING_DIR}")

    # Filter to those not yet done
    to_run = []
    for c in catchments:
        done_marker = Path(OUT_DIR_ROOT) / c / f"{c}_lead_time_forecasts_da.csv"
        if done_marker.exists():
            print(f"  [skip] {c} — already complete")
        else:
            to_run.append(c)
    print(f"[plan] {len(to_run)} to run | concurrent={args.concurrent}")

    if args.dry_run:
        print("[dry-run] would launch:", to_run)
        return

    if not to_run:
        print("[plan] nothing to do — all catchments complete.")
        return

    # Launch with bounded parallelism. ThreadPoolExecutor is fine since each
    # task is a subprocess.run that releases the GIL while waiting.
    with ThreadPoolExecutor(max_workers=args.concurrent) as ex:
        futures = {ex.submit(run_one, c): c for c in to_run}
        completed = 0
        failed = []
        for fut in as_completed(futures):
            cat_id, rc, log = fut.result()
            completed += 1
            status = "OK" if rc == 0 else f"FAIL (rc={rc})"
            print(f"[done {completed}/{len(to_run)}] {cat_id} {status} — log: {log}",
                  flush=True)
            if rc != 0:
                failed.append(cat_id)

    print(f"\n[summary] {len(to_run) - len(failed)} succeeded, {len(failed)} failed")
    if failed:
        print(f"  failed catchments: {failed}")
        sys.exit(1)


if __name__ == "__main__":
    main()
