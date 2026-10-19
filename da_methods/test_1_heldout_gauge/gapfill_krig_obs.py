"""
gapfill_krig_obs.py

Cleans per-catchment Qkrig CSVs:
  1. Drops all rows where date is 2019-01-01 (leading NaN day)
  2. Linear-interpolates any remaining isolated interior NaN hours
  3. Renames columns to match calibration script format: datetime, qkrig, variance
  4. Overwrites files in-place (or writes to --out-dir if given)

Output columns: datetime, qkrig (mm/h), variance

Usage:
    python gapfill_krig_obs.py --obs-dir /mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_with_variance
    python gapfill_krig_obs.py --obs-dir <in-dir> --out-dir <out-dir>   # non-destructive
"""

import argparse
from pathlib import Path

import pandas as pd


def gapfill(path_in: Path, path_out: Path) -> dict:
    df = pd.read_csv(path_in)
    df["time"] = pd.to_datetime(df["time"])

    n_total = len(df)

    # Step 1 — drop leading 2019-01-01 day
    df = df[df["time"].dt.date.astype(str) != "2019-01-01"].copy()
    n_after_drop = len(df)
    dropped = n_total - n_after_drop

    # Step 2 — linear-interpolate interior NaNs (limit=48 to avoid extrapolation)
    nan_before = df["qkrig_mm_hr"].isna().sum()
    df["qkrig_mm_hr"]   = df["qkrig_mm_hr"].interpolate(method="linear", limit=48, limit_direction="both")
    df["qkrig_variance"] = df["qkrig_variance"].interpolate(method="linear", limit=48, limit_direction="both")
    nan_after = df["qkrig_mm_hr"].isna().sum()

    # Step 3 — rename to match calibration script expected format
    df = df.reset_index(drop=True)
    df = df.rename(columns={
        "time": "datetime",
        "qkrig_mm_hr": "qkrig",
        "qkrig_variance": "variance",
    })[["datetime", "qkrig", "variance"]]

    path_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path_out, index=False)

    return {"dropped": dropped, "interpolated": nan_before - nan_after, "remaining_nan": nan_after}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs-dir", required=True)
    parser.add_argument("--out-dir", default=None,
                        help="Output directory. If omitted, overwrites files in-place.")
    args = parser.parse_args()

    obs_dir = Path(args.obs_dir)
    out_dir = Path(args.out_dir) if args.out_dir else obs_dir

    csvs = sorted(obs_dir.glob("cat-*.csv"))
    if not csvs:
        print(f"No cat-*.csv files found in {obs_dir}")
        return

    print(f"Processing {len(csvs)} catchment files...")
    total_interp = 0
    total_remaining = 0

    for csv in csvs:
        out_path = out_dir / csv.name
        stats = gapfill(csv, out_path)
        total_interp    += stats["interpolated"]
        total_remaining += stats["remaining_nan"]
        if stats["remaining_nan"] > 0:
            print(f"  WARNING {csv.name}: {stats['remaining_nan']} NaNs remain after interpolation")

    print(f"\nDone.")
    print(f"  Dropped 2019-01-01 rows : yes (24 rows per file)")
    print(f"  Total hours interpolated: {total_interp}")
    print(f"  Total NaNs remaining    : {total_remaining}")
    if total_remaining > 0:
        print("  WARNING: some NaNs could not be interpolated — check files above")
    else:
        print("  All series are now gap-free.")


if __name__ == "__main__":
    main()
