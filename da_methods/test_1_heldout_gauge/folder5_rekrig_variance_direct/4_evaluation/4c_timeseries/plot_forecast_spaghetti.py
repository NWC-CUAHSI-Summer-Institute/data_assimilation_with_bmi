"""
plot_forecast_spaghetti.py — F5 (re-kriged variance direct)

"All forecasts ending at this hour" — verification-time-fixed spaghetti view.

Window: Sep 24 18 UTC -> Sep 30 06 UTC.
For every initialization time whose 18-hour forecast overlaps this window:
    DA  : shaded min/max band across 20 members + ensemble mean line
    OL  : ensemble mean line only (no shading), dashed

Lines are colored by initialization DATE (7 colors, Sep 24-30).
USGS obs overlaid in black.

Inputs:
    <route-dir>/routed_leadtime_da_full.parquet
    <route-dir>/routed_leadtime_openloop_full.parquet
    Wide format: issue_time, lead_hour, member_00..member_19 (q in m3/s)

Obs:
    /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

Outputs:
    <out-dir>/forecast_spaghetti_helene.png
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

PLOT_START = pd.Timestamp("2024-09-24 18:00:00")
PLOT_END   = pd.Timestamp("2024-09-30 06:00:00")

HELENE_PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_PEAK_END   = pd.Timestamp("2024-09-28 00:00:00")

DATE_COLORS = {
    "2024-09-24": "#1f77b4",
    "2024-09-25": "#ff7f0e",
    "2024-09-26": "#2ca02c",
    "2024-09-27": "#d62728",
    "2024-09-28": "#9467bd",
    "2024-09-29": "#8c564b",
    "2024-09-30": "#e377c2",
}


def load_parquet(path):
    df = pd.read_parquet(path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    df["ens_mean"] = df[member_cols].mean(axis=1)
    df["ens_min"]  = df[member_cols].min(axis=1)
    df["ens_max"]  = df[member_cols].max(axis=1)
    df["valid_time"] = df["issue_time"] + pd.to_timedelta(df["lead_hour"], unit="h")
    return df


def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns
                    if c.lower() in ("datetime", "date", "time", "timestamp"))
    q_col    = next(c for c in df.columns
                    if "q" in c.lower() or "flow" in c.lower()
                    or "discharge" in c.lower())
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if "mm" in q_col.lower():
        series = series * MM_H_TO_M3_S
    return series


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--route-dir", default=DEFAULT_ROUTE_DIR)
    parser.add_argument("--out-dir",   default=None)
    parser.add_argument("--usgs-csv",  default=DEFAULT_USGS_CSV)
    parser.add_argument("--da-name",   default="routed_leadtime_da_full.parquet")
    parser.add_argument("--ol-name",   default="routed_leadtime_openloop_full.parquet")
    args = parser.parse_args()
    out_dir = args.out_dir or args.route_dir
    os.makedirs(out_dir, exist_ok=True)

    print("Loading parquets...")
    da  = load_parquet(os.path.join(args.route_dir, args.da_name))
    ol  = load_parquet(os.path.join(args.route_dir, args.ol_name))
    obs = load_usgs(args.usgs_csv)

    da = da[(da["valid_time"] >= PLOT_START) & (da["valid_time"] <= PLOT_END)]
    ol = ol[(ol["valid_time"] >= PLOT_START) & (ol["valid_time"] <= PLOT_END)]

    issue_times = sorted(da["issue_time"].unique())
    print(f"  {len(issue_times)} initialization times contribute to this window")

    fig, ax = plt.subplots(figsize=(18, 6))

    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color="salmon", alpha=0.12, zorder=0)
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color="firebrick", ha="left", va="top")

    for t0 in issue_times:
        date_str = str(pd.Timestamp(t0).date())
        color = DATE_COLORS.get(date_str, "gray")

        da_t = da[da["issue_time"] == t0].sort_values("valid_time")
        ol_t = ol[ol["issue_time"] == t0].sort_values("valid_time")

        if da_t.empty:
            continue

        ax.fill_between(da_t["valid_time"], da_t["ens_min"], da_t["ens_max"],
                        color=color, alpha=0.06, zorder=2)
        ax.plot(da_t["valid_time"], da_t["ens_mean"],
                color=color, lw=0.8, alpha=0.55, zorder=3)

        if not ol_t.empty:
            ax.plot(ol_t["valid_time"], ol_t["ens_mean"],
                    color=color, lw=0.6, alpha=0.30, linestyle="--", zorder=2)

    obs_w = obs.loc[PLOT_START:PLOT_END]
    ax.plot(obs_w.index, obs_w.values,
            color="black", lw=1.8, zorder=6, label="USGS obs")

    patches = [mpatches.Patch(color=c, label=f"Init {d}")
               for d, c in DATE_COLORS.items()]
    patches.append(plt.Line2D([0], [0], color="black", lw=1.8, label="USGS obs"))
    patches.append(plt.Line2D([0], [0], color="gray", lw=1.2, alpha=0.7,
                               label="DA — ens mean (solid) + spread (shaded)"))
    patches.append(plt.Line2D([0], [0], color="gray", lw=1.0, linestyle="--",
                               alpha=0.5, label="OL — ens mean (dashed)"))
    ax.legend(handles=patches, fontsize=8, loc="upper left",
              frameon=True, framealpha=0.9, ncol=2)

    ax.set_xlim(PLOT_START, PLOT_END)
    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("Discharge (m3/s)", fontsize=11)
    ax.set_title(
        "All 18-hour forecast trajectories ending in this window — USGS 03463300 — F5 (re-kriged variance)\n"
        "Sep 24 18 UTC -> Sep 30 06 UTC  |  DA: shaded band + mean  |  OL: mean dashed",
        fontsize=11,
    )
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    ax.grid(True, alpha=0.2, lw=0.4)
    plt.tight_layout()

    out_path = os.path.join(out_dir, "forecast_spaghetti_helene.png")
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
