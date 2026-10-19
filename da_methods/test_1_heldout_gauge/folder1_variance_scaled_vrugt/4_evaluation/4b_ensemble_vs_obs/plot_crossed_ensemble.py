"""
plot_crossed_ensemble.py  --  4b: 600-member crossed ensemble vs USGS obs.

Loads the crossed ensemble parquet produced by run_crossed_ensemble.py:
    <ensemble-dir>/<cat-id>/<cat-id>_crossed_ensemble.parquet
    Columns: issue_time, lead_hour, member_0000..member_0599

Projects all members to valid_time = issue_time + lead_hour hours.
Aggregates at each valid_time across ALL members from ALL issue times
to build a probabilistic envelope (5th/25th/50th/75th/95th percentiles).

Two output figures:
    4b_crossed_ensemble_helene.png
        -- Full Helene window (Sep 24-30) probabilistic envelope + USGS obs.
           Percentile shading in two layers (5-95 outer, 25-75 inner).
    4b_crossed_ensemble_peak.png
        -- Zoomed to Helene peak (Sep 26 12 UTC -> Sep 28 06 UTC),
           same shading + USGS obs + individual spaghetti from Sep 27 inits.

Usage:
    python3 plot_crossed_ensemble.py \\
        --ensemble-dir /mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble \\
        --cat-id cat-1016300 \\
        --usgs-csv /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv \\
        --out-dir /mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble/cat-1016300
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

DEFAULT_ENSEMBLE_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble"
DEFAULT_USGS_CSV     = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
DEFAULT_CAT_ID       = "cat-1016300"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S       = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

PLOT_START        = pd.Timestamp("2024-09-24 00:00:00")
PLOT_END          = pd.Timestamp("2024-09-28 18:00:00")   # cut trailing spike
PEAK_START        = pd.Timestamp("2024-09-26 12:00:00")
PEAK_END          = pd.Timestamp("2024-09-28 06:00:00")
HELENE_PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_PEAK_END   = pd.Timestamp("2024-09-28 00:00:00")


def load_ensemble(path):
    df = pd.read_parquet(path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    df["valid_time"] = df["issue_time"] + pd.to_timedelta(df["lead_hour"], unit="h")
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    return df, member_cols


def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns
                    if c.lower() in ("datetime", "date", "time", "timestamp"))
    q_col = next(c for c in df.columns
                 if "q" in c.lower() or "flow" in c.lower() or "discharge" in c.lower())
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if "mm" in q_col.lower():
        series = series * MM_H_TO_M3_S
    return series


def build_percentile_envelope(df, member_cols, t_start, t_end):
    """
    For each valid_time in [t_start, t_end], stack all member values
    from all issue times and compute percentiles.
    Returns DataFrame indexed by valid_time with columns p05, p25, p50, p75, p95.
    """
    mask = (df["valid_time"] >= t_start) & (df["valid_time"] <= t_end)
    sub = df[mask].copy()
    sub["valid_time"] = sub["valid_time"].dt.floor("1h")

    records = []
    for vt, grp in sub.groupby("valid_time"):
        vals = grp[member_cols].to_numpy(dtype=float).flatten()
        vals = vals[~np.isnan(vals)] * MM_H_TO_M3_S
        if len(vals) == 0:
            continue
        records.append({
            "valid_time": vt,
            "p05": np.percentile(vals, 5),
            "p25": np.percentile(vals, 25),
            "p50": np.percentile(vals, 50),
            "p75": np.percentile(vals, 75),
            "p95": np.percentile(vals, 95),
            "n":   len(vals),
        })
    return pd.DataFrame(records).set_index("valid_time").sort_index()


def _add_helene_band(ax, label=True):
    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color="salmon", alpha=0.12, zorder=0)
    if label:
        ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 0.97, "Helene peak",
                transform=ax.get_xaxis_transform(),
                fontsize=9, color="firebrick", ha="left", va="top")


def _format_xaxis(ax, major_interval=1, minor_hours=6):
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=major_interval))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=minor_hours))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    ax.grid(True, alpha=0.2, lw=0.4)


def plot_envelope(ax, env, obs_series, t_start, t_end,
                  ens_color="tab:blue", obs_color="black", helene_label=True):
    _add_helene_band(ax, label=helene_label)

    # Outer band: 5-95
    ax.fill_between(env.index, env["p05"], env["p95"],
                    color=ens_color, alpha=0.15, zorder=2,
                    label="5th-95th percentile")
    # Inner band: 25-75
    ax.fill_between(env.index, env["p25"], env["p75"],
                    color=ens_color, alpha=0.30, zorder=3,
                    label="25th-75th percentile")
    # Median
    ax.plot(env.index, env["p50"],
            color=ens_color, lw=2.2, zorder=5,
            label="Ensemble median (600 members)")

    # USGS obs
    obs_w = obs_series.loc[t_start:t_end]
    ax.plot(obs_w.index, obs_w.values,
            color=obs_color, lw=1.8, zorder=6, label="USGS obs")

    ax.set_xlim(t_start, t_end)
    ax.set_ylabel("Discharge (m3/s)", fontsize=11)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ensemble-dir", default=DEFAULT_ENSEMBLE_DIR)
    parser.add_argument("--cat-id",       default=DEFAULT_CAT_ID)
    parser.add_argument("--usgs-csv",     default=DEFAULT_USGS_CSV)
    parser.add_argument("--out-dir",      default=None)
    args = parser.parse_args()

    cat_dir  = os.path.join(args.ensemble_dir, args.cat_id)
    out_dir  = args.out_dir or cat_dir
    os.makedirs(out_dir, exist_ok=True)

    parquet_path = os.path.join(cat_dir, f"{args.cat_id}_crossed_ensemble.parquet")
    print(f"Loading crossed ensemble: {parquet_path}")
    df, member_cols = load_ensemble(parquet_path)
    obs = load_usgs(args.usgs_csv)
    print(f"  {df['issue_time'].nunique()} issue times, "
          f"{len(member_cols)} members per row")

    # ------------------------------------------------------------------ #
    # Figure 4b-full: full Helene window probabilistic envelope            #
    # ------------------------------------------------------------------ #
    print("Building percentile envelope (full window)...")
    env_full = build_percentile_envelope(df, member_cols, PLOT_START, PLOT_END)
    print(f"  Valid times with data: {len(env_full)}, "
          f"median members per timestep: {env_full['n'].median():.0f}")

    fig, ax = plt.subplots(figsize=(17, 7))
    plot_envelope(ax, env_full, obs, PLOT_START, PLOT_END)
    ax.set_xlabel("Date", fontsize=11)
    ax.set_title(
        f"4b -- 600-member crossed ensemble (30 met x 20 hydro-state) | {args.cat_id}\n"
        "Shaded: 5-95th (light) and 25-75th (dark) percentile across all members "
        "and issue times | Sep 24-30 2024 | USGS 03463300",
        fontsize=11,
    )
    ax.legend(fontsize=10, loc="upper left", frameon=True, framealpha=0.9)
    _format_xaxis(ax)
    plt.tight_layout()
    out1 = os.path.join(out_dir, f"{args.cat_id}_4b_crossed_ensemble_helene.png")
    plt.savefig(out1, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out1}")

    # ------------------------------------------------------------------ #
    # Figure 4b-peak: zoomed to Helene peak window                        #
    # ------------------------------------------------------------------ #
    print("Building percentile envelope (peak window)...")
    env_peak = build_percentile_envelope(df, member_cols, PEAK_START, PEAK_END)

    fig, (ax_main, ax_zoom) = plt.subplots(
        2, 1, figsize=(14, 11),
        gridspec_kw={"height_ratios": [1, 1.4]}, sharex=False)

    # Top panel: full window context (smaller)
    plot_envelope(ax_main, env_full, obs, PLOT_START, PLOT_END, helene_label=False)
    ax_main.set_title("Full Helene window context (Sep 24 - Sep 28 18 UTC)", fontsize=10)
    _format_xaxis(ax_main)
    # Shade the zoom region on context panel
    ax_main.axvspan(PEAK_START, PEAK_END, color="gold", alpha=0.18, zorder=1)
    ax_main.text(PEAK_START + pd.Timedelta(hours=1), 0.97, "zoom",
                 transform=ax_main.get_xaxis_transform(),
                 fontsize=8, color="goldenrod", ha="left", va="top")

    # Bottom panel: peak zoom
    plot_envelope(ax_zoom, env_peak, obs, PEAK_START, PEAK_END, helene_label=True)
    ax_zoom.set_xlabel("Date", fontsize=11)
    ax_zoom.set_title(
        "Helene peak zoom (Sep 26 12 UTC - Sep 28 06 UTC)", fontsize=10)
    _format_xaxis(ax_zoom, major_interval=1, minor_hours=3)
    ax_zoom.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %HUTC"))
    ax_zoom.legend(fontsize=10, loc="lower right", frameon=True, framealpha=0.9)

    fig.suptitle(
        f"4b -- 600-member crossed ensemble vs USGS 03463300 | {args.cat_id}\n"
        "30 met forcing draws x 20 DA analysis state draws | R = (0.10·Q)² + 0.001·σ²_krig  (Vrugt)",
        fontsize=12, y=1.01,
    )
    plt.tight_layout()
    out2 = os.path.join(out_dir, f"{args.cat_id}_4b_crossed_ensemble_peak.png")
    plt.savefig(out2, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out2}")


if __name__ == "__main__":
    main()
