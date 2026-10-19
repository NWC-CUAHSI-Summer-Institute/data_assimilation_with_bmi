#!/usr/bin/env python3
"""
route_forecast_al_gauges.py -- Run route_forecast.py for all 12 Alabama gauges.

Thin wrapper: calls route_forecast.py once per gauge, substituting the correct
GPKG and output directory.  All other arguments (fc-dir, da-traj-dir, etc.)
are forwarded unchanged.

USAGE
  # route both DA and open-loop forecasts for all 12 gauges (20 procs each)
  python route_forecast_al_gauges.py \\
      --fc-dir        /mnt/disk2/svyas_clean_run/da_al_storm_2026_forecast_sweep \\
      --da-traj-dir   /mnt/disk2/svyas_clean_run/da_al_storm_2026_traj_da \\
      --ol-traj-dir   /mnt/disk2/svyas_clean_run/da_al_storm_2026_traj_ol \\
      --gpkg-dir      /mnt/disk2/al_storm_2026/gpkg \\
      --out-base      /home/svyas/qkrig_DA/results/Alabama_gauges \\
      --procs         20

  # dry-run
  python route_forecast_al_gauges.py --dry-run ...
"""

import argparse
import subprocess
import sys
from pathlib import Path

AL_GAUGES = [
    "02361000", "02369800", "02371500", "02372250", "02374500", "02408540",
    "02422500", "02450250", "02464000", "02465493", "02469800", "03574500",
]

ROUTE_SCRIPT = Path(__file__).parent / "route_forecast.py"


def main():
    parser = argparse.ArgumentParser(
        description="Route Alabama storm ensemble forecasts for all 12 gauges.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--fc-dir",      required=True,
                        help="forecast dir from run_cfe_enkf_forecast_hrrr_multihour_al_storm.py")
    parser.add_argument("--da-traj-dir", required=True,
                        help="DA member_q trajectories for warm-up")
    parser.add_argument("--ol-traj-dir", required=True,
                        help="open-loop member_q trajectories")
    parser.add_argument("--gpkg-dir",    default="/mnt/disk2/al_storm_2026/gpkg",
                        help="directory containing gage-<gauge>_subset.gpkg files")
    parser.add_argument("--procs",       type=int, default=20)
    parser.add_argument("--warmup-h",    type=int, default=48)
    parser.add_argument("--scenario",    default="both",
                        choices=["da", "openloop", "both"])
    parser.add_argument("--gauge",       default=None,
                        help="route one gauge only (default: all 12)")
    parser.add_argument("--dry-run",     action="store_true")
    args = parser.parse_args()

    gauges = [args.gauge] if args.gauge else AL_GAUGES
    failed = []

    for gauge in gauges:
        gpkg_matches = sorted(Path(args.gpkg_dir).glob(f"*{gauge}*.gpkg"))
        if not gpkg_matches:
            print(f"\n[{gauge}] SKIP — no GPKG in {args.gpkg_dir}")
            continue

        gpkg = str(gpkg_matches[0])
        # output lands in <fc-dir>/<cat>/routed/ per route_forecast.py convention
        routed_dir = str(Path(args.fc_dir) / "routed")

        cmd = [
            sys.executable, str(ROUTE_SCRIPT),
            "--fc-dir",      args.fc_dir,
            "--da-traj-dir", args.da_traj_dir,
            "--ol-traj-dir", args.ol_traj_dir,
            "--gpkg",        gpkg,
            "--procs",       str(args.procs),
            "--warmup-h",    str(args.warmup_h),
            "--scenario",    args.scenario,
        ]

        print(f"\n{'='*70}")
        print(f"Routing gauge {gauge}  |  GPKG: {Path(gpkg).name}")
        print(f"  cmd: {' '.join(cmd)}")
        if args.dry_run:
            print("  [dry-run — not executed]")
            continue

        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"  ERROR: route_forecast.py exited {result.returncode} for {gauge}")
            failed.append(gauge)
        else:
            # move outputs into gauge-specific subdir so gauges don't overwrite each other
            gauge_out = Path(routed_dir) / gauge
            gauge_out.mkdir(parents=True, exist_ok=True)
            for fname in ["routed_forecast_leadtime_da.parquet",
                          "routed_forecast_leadtime_openloop.parquet"]:
                src = Path(routed_dir) / fname
                if src.exists():
                    src.rename(gauge_out / fname)
            print(f"  Done -> {gauge_out}/")

    print(f"\n{'='*70}")
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"All {len(gauges)} gauge(s) routed. Output: {args.fc_dir}/routed/")
    print("Next: plot_forecast_all_gauges.py reads routed/routed_forecast_leadtime_da.parquet")


if __name__ == "__main__":
    main()
