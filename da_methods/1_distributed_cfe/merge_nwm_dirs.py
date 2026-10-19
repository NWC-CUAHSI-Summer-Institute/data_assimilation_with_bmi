#!/usr/bin/env python3
"""
merge_nwm_dirs.py -- Merge two NWM operational forcing dirs into one.

Concatenates base + ext CSVs for each catchment, deduplicates by time
(keep='last' = most recent analysis cycle wins), sorts by time.

USAGE
  python merge_nwm_dirs.py \
      --base-dir  /mnt/disk2/al_storm_2026/nwm_operational/wy2027/forcings/current/forcings \
      --ext-dir   /mnt/disk2/al_storm_2026/nwm_operational/wy2027/isaias/forcings \
      --out-dir   /mnt/disk2/svyas_clean_run/nwm_merged_oct10 \
      --procs     20
"""
import argparse
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

G = {}


def merge_one(cat_id):
    base_dir = Path(G["base_dir"])
    ext_dir  = Path(G["ext_dir"])
    out_dir  = Path(G["out_dir"])

    out_csv = out_dir / f"{cat_id}.csv"
    if out_csv.exists():
        return cat_id, "skip"

    base_csv = base_dir / f"{cat_id}.csv"
    ext_csv  = ext_dir  / f"{cat_id}.csv"

    if not base_csv.exists():
        return cat_id, f"ERROR: base CSV not found: {base_csv}"

    df = pd.read_csv(base_csv)
    if ext_csv.exists():
        df = pd.concat([df, pd.read_csv(ext_csv)], ignore_index=True)

    df = (df.drop_duplicates(subset="time", keep="last")
            .sort_values("time")
            .reset_index(drop=True))
    df.to_csv(out_csv, index=False)
    return cat_id, f"ok  {df['time'].iloc[0]} → {df['time'].iloc[-1]}  ({len(df)} rows)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--ext-dir",  required=True)
    parser.add_argument("--out-dir",  required=True)
    parser.add_argument("--procs",    type=int, default=20)
    parser.add_argument("--force",    action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cats = sorted(p.stem for p in Path(args.base_dir).glob("cat-*.csv"))
    print(f"Merging {len(cats)} catchments | procs={args.procs}")

    if args.force:
        for f in out_dir.glob("cat-*.csv"):
            f.unlink()

    G["base_dir"] = args.base_dir
    G["ext_dir"]  = args.ext_dir
    G["out_dir"]  = args.out_dir

    errors = []
    with Pool(args.procs) as pool:
        for cat_id, msg in pool.imap_unordered(merge_one, cats):
            if msg.startswith("ERROR"):
                print(f"  {cat_id}: {msg}")
                errors.append(cat_id)

    done = sum(1 for f in out_dir.glob("cat-*.csv"))
    print(f"\nDone: {done} CSVs in {args.out_dir}")
    if errors:
        print(f"ERRORS: {', '.join(errors)}")


if __name__ == "__main__":
    main()
