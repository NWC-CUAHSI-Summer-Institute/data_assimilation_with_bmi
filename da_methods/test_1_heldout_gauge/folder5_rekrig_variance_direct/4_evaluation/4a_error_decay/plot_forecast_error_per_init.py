"""
plot_forecast_error_per_init.py — F5 re-kriged variance direct (1 gauge holdout).

Error decay by initialization time for the Helene window (Sep 24-28 2024).

For each initialization time t0 in the Helene window:
    error[lead] = ensemble_mean(q_gauge_m3s at t0+lead) - USGS_obs(t0+lead)

Plotted as:
    DA  : thin colored lines (one per init time, colored by date) + thick mean
    OL  : thin gray dashed lines + thick gray dashed mean

x-axis: forecast lead hour (1 -> 18)
y-axis: signed error (m³/s), positive = forecast too high

The expected signal: DA error is small at lead 1 (just assimilated), grows
and converges toward the OL error curve by lead 18.

Inputs:
    <route-dir>/routed_leadtime_da_full.parquet
    <route-dir>/routed_leadtime_openloop_full.parquet

Obs:
    /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

Outputs:
    <out-dir>/forecast_error_per_init_helene.png   (signed error)
    <out-dir>/forecast_mae_per_lead_helene.png     (mean absolute error per lead)
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

HELENE_INIT_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_INIT_END   = pd.Timestamp("2024-09-28 23:00:00")

DATE_COLORS = {
    "2024-09-24": "#1f77b4",
    "2024-09-25": "#ff7f0e",
    "2024-09-26": "#2ca02c",
    "2024-09-27": "#d62728",
    "2024-09-28": "#9467bd",
}


def load_parquet(path):
    df = pd.read_parquet(path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    df["ens_mean"] = df[member_cols].mean(axis=1)
    df["valid_time"] = df["issue_time"] + pd.to_timedelta(df["lead_hour"], unit="h")
    return df[["issue_time", "lead_hour", "valid_time", "ens_mean"]]


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


def build_error_table(df, obs_series, init_start, init_end):
    df = df[(df["issue_time"] >= init_start) & (df["issue_time"] <= init_end)].copy()
    df["obs"]   = df["valid_time"].map(obs_series)
    df["error"] = df["ens_mean"] - df["obs"]
    return df.dropna(subset=["obs", "error"])


def plot_error(ax, err_df, color, alpha_thin, lw_thin, linestyle, label_prefix):
    leads = sorted(err_df["lead_hour"].unique())
    all_curves = []

    for t0, grp in err_df.groupby("issue_time"):
        date_str = str(pd.Timestamp(t0).date())
        c = DATE_COLORS.get(date_str, color)
        grp_sorted = grp.sort_values("lead_hour")
        curve = grp_sorted.set_index("lead_hour")["error"].reindex(leads).values
        ax.plot(leads, curve,
                color=c, lw=lw_thin, alpha=alpha_thin,
                linestyle=linestyle, zorder=2)
        all_curves.append(curve)

    if all_curves:
        mean_curve = np.nanmean(np.array(all_curves), axis=0)
        ax.plot(leads, mean_curve,
                color=color, lw=2.8, alpha=0.95,
                linestyle=linestyle, zorder=5,
                label=f"{label_prefix} — mean across all init times")


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
    da = load_parquet(os.path.join(args.route_dir, args.da_name))
    ol = load_parquet(os.path.join(args.route_dir, args.ol_name))
    obs = load_usgs(args.usgs_csv)
    print(f"  DA issue_times: {da['issue_time'].nunique()}  "
          f"OL issue_times: {ol['issue_time'].nunique()}")

    da_err = build_error_table(da, obs, HELENE_INIT_START, HELENE_INIT_END)
    ol_err = build_error_table(ol, obs, HELENE_INIT_START, HELENE_INIT_END)
    print(f"  DA init times in Helene window: {da_err['issue_time'].nunique()}")

    leads = sorted(da_err["lead_hour"].unique())

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.5, zorder=1)

    plot_error(ax, ol_err, color="tab:gray",   alpha_thin=0.12, lw_thin=0.7,
               linestyle="--", label_prefix="Open-loop")
    plot_error(ax, da_err, color="tab:purple", alpha_thin=0.18, lw_thin=0.8,
               linestyle="-",  label_prefix="DA")

    patches = [mpatches.Patch(color=c, label=f"Init {d}")
               for d, c in DATE_COLORS.items()]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + patches,
              fontsize=8, loc="upper left", frameon=True, framealpha=0.9, ncol=2)

    ax.set_xlabel("Forecast lead hour", fontsize=11)
    ax.set_ylabel("Error: forecast − USGS obs (m³/s)", fontsize=11)
    ax.set_xticks(np.arange(1, 19))
    ax.set_title(
        "Forecast error vs lead time — per initialization time, Helene window\n"
        "F5 (re-kriged variance) | DA (purple solid) vs Open-loop (gray dashed) | "
        "Sep 24–28 2024 | USGS 03463300",
        fontsize=11,
    )
    ax.grid(True, alpha=0.25, lw=0.4)
    plt.tight_layout()
    out1 = os.path.join(out_dir, "forecast_error_per_init_helene.png")
    plt.savefig(out1, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out1}")

    fig, ax = plt.subplots(figsize=(12, 6))

    def mean_abs_error_by_lead(err_df):
        return err_df.groupby("lead_hour")["error"].apply(
            lambda x: float(np.nanmean(np.abs(x)))
        )

    da_mae = mean_abs_error_by_lead(da_err)
    ol_mae = mean_abs_error_by_lead(ol_err)

    ax.plot(da_mae.index, da_mae.values,
            color="tab:purple", lw=2.4, marker="o", label="DA — mean |error|")
    ax.plot(ol_mae.index, ol_mae.values,
            color="tab:gray",   lw=2.4, marker="s", linestyle="--",
            label="Open-loop — mean |error|")

    ax.set_xlabel("Forecast lead hour", fontsize=11)
    ax.set_ylabel("Mean |error| vs USGS obs (m³/s)", fontsize=11)
    ax.set_xticks(np.arange(1, 19))
    ax.set_title(
        "Mean absolute forecast error vs lead time — Helene window\n"
        "F5 (re-kriged variance) | DA (purple) vs Open-loop (gray) | "
        "Sep 24–28 2024 | USGS 03463300",
        fontsize=11,
    )
    ax.legend(fontsize=10, loc="upper left", frameon=True, framealpha=0.9)
    ax.grid(True, alpha=0.25, lw=0.4)
    plt.tight_layout()
    out2 = os.path.join(out_dir, "forecast_mae_per_lead_helene.png")
    plt.savefig(out2, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out2}")


if __name__ == "__main__":
    main()
