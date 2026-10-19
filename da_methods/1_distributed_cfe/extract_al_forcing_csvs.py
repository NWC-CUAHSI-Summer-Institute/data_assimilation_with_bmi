#!/usr/bin/env python3
"""Extract per-catchment forcing CSVs from NWM retro forcings.nc files.

Reads:  /mnt/disk2/al_storm_2026/nwm_retro/gage-XXXXXXXX/forcings/forcings.nc
Writes: /mnt/disk2/svyas_clean_run/nwm_retro_al_storm_2026/{cat-ID}.csv

Output CSV columns match the existing heldin calibration format:
  time, UGRD_10maboveground, precip_rate, DSWRF_surface, TMP_2maboveground,
  SPFH_2maboveground, VGRD_10maboveground, DLWRF_surface, PRES_surface, APCP_surface
"""
import sys
from pathlib import Path
import datetime
import numpy as np
import xarray as xr

RETRO_DIR = Path("/mnt/disk2/al_storm_2026/nwm_retro")
OUT_DIR   = Path("/mnt/disk2/svyas_clean_run/nwm_retro_al_storm_2026")

COLS = [
    "UGRD_10maboveground", "precip_rate", "DSWRF_surface",
    "TMP_2maboveground", "SPFH_2maboveground", "VGRD_10maboveground",
    "DLWRF_surface", "PRES_surface", "APCP_surface",
]

OUT_DIR.mkdir(parents=True, exist_ok=True)

gauge_dirs = sorted(RETRO_DIR.glob("gage-*/forcings/forcings.nc"))
print(f"Found {len(gauge_dirs)} gauge NetCDF files")

total_written = 0
total_skipped = 0

for nc_path in gauge_dirs:
    gauge = nc_path.parts[-3]  # gage-02361000
    print(f"\n── {gauge} ──")
    ds = xr.open_dataset(str(nc_path))

    cat_ids  = ds["ids"].values          # shape: (catchment-id,)
    time_arr = ds["Time"].values[0, :]   # unix timestamps, shape: (time,)

    # decode unix timestamps to datetime strings once
    times = [
        datetime.datetime.utcfromtimestamp(int(t)).strftime("%Y-%m-%d %H:%M:%S")
        for t in time_arr
    ]

    n_cats = len(cat_ids)
    for i, cat_id in enumerate(cat_ids):
        cat_id = str(cat_id).strip()
        out_csv = OUT_DIR / f"{cat_id}.csv"

        if out_csv.exists():
            total_skipped += 1
            continue

        # build dataframe-like output directly to CSV for speed
        rows = ["time," + ",".join(COLS)]
        data = {col: ds[col].values[i, :] for col in COLS}

        for j, t in enumerate(times):
            vals = ",".join(f"{data[col][j]:.7g}" for col in COLS)
            rows.append(f"{t},{vals}")

        out_csv.write_text("\n".join(rows) + "\n")
        total_written += 1

        if (i + 1) % 50 == 0 or (i + 1) == n_cats:
            print(f"  {i+1}/{n_cats} catchments written")

    ds.close()

print(f"\n=== Done ===")
print(f"  Written: {total_written}  |  Skipped (already exist): {total_skipped}")
print(f"  Output: {OUT_DIR}")
