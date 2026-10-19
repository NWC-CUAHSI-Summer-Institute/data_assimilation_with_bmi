#!/usr/bin/env python3
"""
run_daily_update_al_storm.py -- Daily re-issue of Alabama storm forecast.

Runs all four steps for one new issue time:
  1. Download HRRR for the new day (if not already cached)
  2. Extend EnKF DA to the new issue time
  3. Run 18h ensemble forecast from the new HRRR init
  4. Route through T-Route for all 12 gauges → timestamped parquets

The resulting per-gauge parquets are named
  gauge_<id>_forecast18h_<YYYYMMDDHHH>.parquet
and are automatically picked up by plot_forecast_all_gauges.py to build the
multi-colored spaghetti fan plot (one fan per issue time, colored cool→warm).

USAGE
    # Re-issue at Oct 9 12z
    python run_daily_update_al_storm.py --issue-time 2026-10-09T12

    # Dry-run: show commands without executing
    python run_daily_update_al_storm.py --issue-time 2026-10-09T12 --dry-run

    # Skip individual steps (e.g. if HRRR already downloaded):
    python run_daily_update_al_storm.py --issue-time 2026-10-09T12 --skip-download

ADDING ISSUE TIMES BUILDS THE HELENE-STYLE MULTI-FAN PLOT
    Issue Oct  8 12z → blue fan   (first run, already done)
    Issue Oct  9 12z → cyan fan
    Issue Oct 10 12z → green fan
    Issue Oct 11 12z → orange fan  (storm making landfall)
    Issue Oct 12 12z → red fan
"""

import argparse
import subprocess
import sys
from pathlib import Path

# ── Server paths (must match the individual scripts) ──────────────────────────
HRRR_DIR    = "/mnt/disk2/svyas_clean_run/hrrr_al_storm_2026"
DA_DIR      = "/mnt/disk2/svyas_clean_run/da_al_storm_2026"
FC_DIR      = "/mnt/disk2/svyas_clean_run/forecast_al_storm_2026"
OUT_BASE    = "/home/svyas/qkrig_DA/results/Alabama_gauges"
SCRIPT_DIR  = Path(__file__).parent
ROUTE_DIR   = SCRIPT_DIR.parent / "2_troute_routing"

DOWNLOAD_SCRIPT = str(SCRIPT_DIR / "download_hrrr_al_storm.py")
DA_SCRIPT       = str(SCRIPT_DIR / "run_cfe_enkf_al_storm.py")
FORECAST_SCRIPT = str(SCRIPT_DIR / "run_forecast_al_storm.py")
ROUTE_SCRIPT    = str(ROUTE_DIR  / "route_forecast_al_storm.py")


def run(cmd, dry_run=False, step=""):
    print(f"\n{'─'*70}")
    print(f"STEP: {step}")
    print(f"CMD : {' '.join(str(c) for c in cmd)}")
    print('─'*70)
    if dry_run:
        print("[dry-run — not executed]")
        return 0
    result = subprocess.run([str(c) for c in cmd])
    if result.returncode != 0:
        print(f"\nERROR: step '{step}' exited with code {result.returncode}")
        sys.exit(result.returncode)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Daily re-issue of Alabama storm 18h ensemble forecast.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--issue-time", required=True,
                        help="Forecast issue time, e.g. 2026-10-09T12")
    parser.add_argument("--dry-run",        action="store_true")
    parser.add_argument("--skip-download",  action="store_true",
                        help="skip HRRR download (if already done for this day)")
    parser.add_argument("--skip-da",        action="store_true",
                        help="skip EnKF DA extension (use existing member states)")
    parser.add_argument("--skip-forecast",  action="store_true",
                        help="skip CFE ensemble forecast run")
    parser.add_argument("--skip-routing",   action="store_true",
                        help="skip T-Route routing")
    parser.add_argument("--n-members",  type=int, default=600)
    parser.add_argument("--hrrr-dir",   default=HRRR_DIR)
    parser.add_argument("--da-dir",     default=DA_DIR)
    parser.add_argument("--fc-dir",     default=FC_DIR)
    parser.add_argument("--out-base",   default=OUT_BASE)
    args = parser.parse_args()

    import pandas as pd
    issue = pd.Timestamp(args.issue_time)
    issue_date  = issue.strftime("%Y-%m-%dT%H")      # "2026-10-09T12"
    issue_ymdh  = issue.strftime("%Y%m%d%H")         # "2026100912"  (output tag)
    da_end_str  = issue.strftime("%Y-%m-%d %H:%M:%S")  # "2026-10-09 12:00:00"
    day_str     = issue.strftime("%Y%m%d")            # "20261009"

    print(f"\n{'═'*70}")
    print(f"Alabama storm daily update")
    print(f"  Issue time : {issue_date}")
    print(f"  DA cutoff  : {da_end_str}")
    print(f"  Output tag : {issue_ymdh}")
    print(f"{'═'*70}")

    # ── 1. Download HRRR for the new day ─────────────────────────────────────
    if not args.skip_download:
        run([
            sys.executable, DOWNLOAD_SCRIPT,
            "--out-dir", args.hrrr_dir,
            "--init-start", f"{day_str[:4]}-{day_str[4:6]}-{day_str[6:]}T00",
            "--init-end",   issue_date,
        ], dry_run=args.dry_run, step="1/4 HRRR download")
    else:
        print("\nStep 1/4: HRRR download — SKIPPED")

    # ── 2. Extend EnKF DA to new issue time ───────────────────────────────────
    if not args.skip_da:
        run([
            sys.executable, DA_SCRIPT,
            "--cat-id", "all",
            "--hrrr-dir", args.hrrr_dir,
            "--out-dir",  args.da_dir,
            "--test-end", da_end_str,
            "--rebuild-hrrr-cache",          # force re-concat with new day's data
            "--n-members", str(args.n_members),
        ], dry_run=args.dry_run, step="2/4 EnKF DA continuation")
    else:
        print("\nStep 2/4: EnKF DA — SKIPPED")

    # ── 3. Run 18h ensemble forecast ──────────────────────────────────────────
    if not args.skip_forecast:
        run([
            sys.executable, FORECAST_SCRIPT,
            "--cat-id", "all",
            "--da-dir",  args.da_dir,
            "--out-dir", args.fc_dir,
            "--forecast-init-time",  issue_date,
            "--forecast-state-time", da_end_str,
            "--out-tag", issue_ymdh,
            "--n-members", str(args.n_members),
        ], dry_run=args.dry_run, step="3/4 18h ensemble forecast")
    else:
        print("\nStep 3/4: Forecast — SKIPPED")

    # ── 4. Route through T-Route for all 12 gauges ───────────────────────────
    if not args.skip_routing:
        run([
            sys.executable, ROUTE_SCRIPT,
            "--fc-dir",  args.fc_dir,
            "--out-base", args.out_base,
            "--fc-tag",  issue_ymdh,
            "--out-tag", issue_ymdh,
        ], dry_run=args.dry_run, step="4/4 T-Route routing")
    else:
        print("\nStep 4/4: Routing — SKIPPED")

    print(f"\n{'═'*70}")
    print(f"Done. Issue {issue_date} complete.")
    print(f"  Routed parquets: {args.out_base}/<gauge>/gauge_<id>_forecast18h_{issue_ymdh}.parquet")
    print(f"  Re-run plot    : python plot_forecast_all_gauges.py")
    print(f"{'═'*70}")


if __name__ == "__main__":
    main()
