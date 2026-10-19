"""
Sanity check: NWM operational precip totals during Hurricane Helene across
all 21 catchments of USGS gauge 03463300 (South Toe River near Celo, NC).

Reads the precip_mm_h column from each catchment's <cat>_test_results.csv
(produced by the production EnKF run), sums hourly precip over the Helene
event window, and prints per-catchment totals + per-day breakdown.

Reference event totals for the basin (Helene, Sep 25-28 2024):
  - Asheville area:     ~300-350 mm 3-day total
  - Mt. Mitchell area:  500+ mm reported (South Toe basin is right here)
  - General western NC: 250-500+ mm common at elevation

If the NWM operational totals across the 21 catchments are < 200 mm typical,
the forcing is the dominant bottleneck and no state-space DA can manufacture
the missing runoff. That shifts the diagnosis from "DA tuning"
to "forcing quality / precip ensemble".
"""
import glob
import os
import pandas as pd
import numpy as np

DA_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"

EVENT_START = pd.Timestamp("2024-09-25 00:00:00")
EVENT_END   = pd.Timestamp("2024-09-28 00:00:00")

DAY_WINDOWS = [
    ("sep_25_mm", "2024-09-25 00:00:00", "2024-09-26 00:00:00"),
    ("sep_26_mm", "2024-09-26 00:00:00", "2024-09-27 00:00:00"),
    ("sep_27_mm", "2024-09-27 00:00:00", "2024-09-28 00:00:00"),
]


def main():
    cat_csvs = sorted(glob.glob(os.path.join(DA_DIR, "cat-*", "cat-*_test_results.csv")))
    if not cat_csvs:
        print(f"No catchment CSVs found under {DA_DIR}.")
        return

    rows = []
    for csv in cat_csvs:
        cat_id = os.path.basename(csv).replace("_test_results.csv", "")
        df = pd.read_csv(csv, parse_dates=["date"])
        if "precip_mm_h" not in df.columns:
            print(f"  WARN: {cat_id} missing precip_mm_h column — skipping")
            continue
        ev = df[(df["date"] >= EVENT_START) & (df["date"] < EVENT_END)]
        total = float(ev["precip_mm_h"].sum())
        peak  = float(ev["precip_mm_h"].max())
        peak_time = ev.loc[ev["precip_mm_h"].idxmax(), "date"] if len(ev) > 0 else pd.NaT
        row = {
            "cat_id": cat_id,
            "total_3day_mm": total,
            "peak_hourly_mm_h": peak,
            "peak_time": peak_time,
        }
        for label, s, e in DAY_WINDOWS:
            day = df[(df["date"] >= pd.Timestamp(s)) & (df["date"] < pd.Timestamp(e))]
            row[label] = float(day["precip_mm_h"].sum())
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("total_3day_mm", ascending=False)

    print("=" * 90)
    print("NWM operational precip totals across 21 catchments — Helene window (Sep 25-27 2024)")
    print("=" * 90)
    with pd.option_context("display.max_rows", None,
                            "display.float_format", "{:.1f}".format,
                            "display.width", 160):
        print(out.to_string(index=False))
    print("=" * 90)
    print(f"Cross-catchment 3-day total (Sep 25-27):")
    print(f"  min  : {out['total_3day_mm'].min():.1f} mm")
    print(f"  mean : {out['total_3day_mm'].mean():.1f} mm")
    print(f"  max  : {out['total_3day_mm'].max():.1f} mm")
    print()
    print(f"Peak hourly intensity across catchments:")
    print(f"  min  : {out['peak_hourly_mm_h'].min():.2f} mm/h")
    print(f"  mean : {out['peak_hourly_mm_h'].mean():.2f} mm/h")
    print(f"  max  : {out['peak_hourly_mm_h'].max():.2f} mm/h")
    print()
    print("Reference (observed event totals):")
    print("  - Western NC widely reported 250-500+ mm 3-day totals")
    print("  - South Toe basin sits in the high-rainfall corridor")
    print("  - Asheville (downstream, lower elevation) reported ~300-350 mm")
    print("=" * 90)

    out_csv = os.path.join(DA_DIR, "_helene_precip_totals_all_cats.csv")
    out.to_csv(out_csv, index=False)
    print(f"\nSaved per-catchment totals: {out_csv}")


if __name__ == "__main__":
    main()
