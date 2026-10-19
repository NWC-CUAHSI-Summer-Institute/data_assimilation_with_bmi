"""
plot_forecast_error_fixed_target.py  —  4a: error decay, fixed-target-time view.
F1 variance-scaled Vrugt (1 gauge holdout).

For each target verification time T in the Helene peak window:
    Collect all forecasts that verify AT T:
        issue_time = T - lead_hour*1h, for lead in 1..18
    error[lead] = ensemble_mean(q at T, initialized T-lead) - USGS_obs(T)

This gives the correct operational picture:
    - lead 1  = initialized 1 hr before T  (DA just ran  -> small error)
    - lead 18 = initialized 18 hr before T (DA long ago  -> error ~ open loop)

Two panels:
    Top : signed error (m³/s) vs lead hour, one curve per target time
    Bot : same for open-loop
Plus a summary panel: mean across all target times, DA vs OL.

Outputs:
    <out-dir>/error_fixed_target_helene.png     (per-target spaghetti, DA vs OL)
    <out-dir>/error_fixed_target_mean.png       (mean across targets, DA vs OL)
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt_leadtime_routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

# Target verification times: hourly through the Helene peak window
TARGET_START = pd.Timestamp("2024-09-26 18:00:00")
TARGET_END   = pd.Timestamp("2024-09-28 06:00:00")


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


def build_fixed_target_errors(df, obs_series, target_times):
    df_idx = df.set_index(["issue_time", "lead_hour"])["ens_mean"]
    results = {}
    for T in target_times:
        obs_val = obs_series.get(T, np.nan)
        if np.isnan(obs_val):
            continue
        curve = {}
        for lead in range(1, 19):
            t0 = T - pd.Timedelta(hours=lead)
            try:
                q_fc = df_idx.loc[(t0, lead)]
                curve[lead] = float(q_fc) - obs_val
            except KeyError:
                curve[lead] = np.nan
        results[T] = curve
    return results


def plot_spaghetti(ax, error_dict, color_da, label_prefix, linestyle="-", lw=0.9, alpha=0.35):
    leads = list(range(1, 19))
    all_curves = []
    target_times = sorted(error_dict.keys())
    cmap = cm.get_cmap("plasma", len(target_times))

    for i, T in enumerate(target_times):
        curve = [error_dict[T].get(l, np.nan) for l in leads]
        ax.plot(leads, curve,
                color=cmap(i), lw=lw, alpha=alpha,
                linestyle=linestyle, zorder=2)
        all_curves.append(curve)

    if all_curves:
        mean_curve = np.nanmean(np.array(all_curves), axis=0)
        ax.plot(leads, mean_curve,
                color=color_da, lw=2.8, alpha=0.95,
                linestyle=linestyle, zorder=5,
                label=f"{label_prefix} — mean across targets")
    return all_curves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--route-dir", default=DEFAULT_ROUTE_DIR)
    parser.add_argument("--out-dir",   default=None)
    parser.add_argument("--usgs-csv",  default=DEFAULT_USGS_CSV)
    parser.add_argument("--da-name",   default="routed_leadtime_da_full.parquet")
    parser.add_argument("--ol-name",   default="routed_leadtime_openloop_full.parquet")
    parser.add_argument("--target-start", default=str(TARGET_START))
    parser.add_argument("--target-end",   default=str(TARGET_END))
    args = parser.parse_args()
    out_dir = args.out_dir or args.route_dir
    os.makedirs(out_dir, exist_ok=True)

    print("Loading parquets...")
    da = load_parquet(os.path.join(args.route_dir, args.da_name))
    ol = load_parquet(os.path.join(args.route_dir, args.ol_name))
    obs = load_usgs(args.usgs_csv)

    target_times = pd.date_range(args.target_start, args.target_end, freq="1h")
    print(f"  Target verification times: {len(target_times)} "
          f"({target_times[0]} → {target_times[-1]})")

    da_errors = build_fixed_target_errors(da, obs, target_times)
    ol_errors = build_fixed_target_errors(ol, obs, target_times)
    print(f"  Targets with obs: DA={len(da_errors)}  OL={len(ol_errors)}")

    leads = list(range(1, 19))

    fig, (ax_da, ax_ol) = plt.subplots(2, 1, figsize=(13, 10), sharex=True, sharey=True)

    ax_da.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.4)
    plot_spaghetti(ax_da, da_errors, "tab:blue", "DA", linestyle="-")
    ax_da.set_ylabel("Error: forecast − USGS obs (m³/s)", fontsize=10)
    ax_da.set_title("DA — error at each lead for fixed target times (Helene peak window)", fontsize=11)
    ax_da.grid(True, alpha=0.2)
    ax_da.legend(fontsize=9)

    ax_ol.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.4)
    plot_spaghetti(ax_ol, ol_errors, "tab:gray", "Open-loop", linestyle="--")
    ax_ol.set_ylabel("Error: forecast − USGS obs (m³/s)", fontsize=10)
    ax_ol.set_xlabel("Forecast lead hour (hours before target)", fontsize=11)
    ax_ol.set_title("Open-loop — error at each lead for fixed target times", fontsize=11)
    ax_ol.set_xticks(leads)
    ax_ol.grid(True, alpha=0.2)
    ax_ol.legend(fontsize=9)

    fig.suptitle(
        "Forecast error vs lead time — fixed verification time, Helene peak window\n"
        f"F1 (variance-scaled Vrugt) | Target: {TARGET_START.strftime('%b %d %H UTC')} → "
        f"{TARGET_END.strftime('%b %d %H UTC')}  |  USGS 03463300",
        fontsize=12, y=1.01,
    )
    plt.tight_layout()
    out1 = os.path.join(out_dir, "error_fixed_target_helene.png")
    plt.savefig(out1, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out1}")

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.4)

    def mean_curve(error_dict):
        arr = np.array([
            [error_dict[T].get(l, np.nan) for l in leads]
            for T in sorted(error_dict.keys())
        ])
        return np.nanmean(arr, axis=0), np.nanstd(arr, axis=0)

    da_mean, da_std = mean_curve(da_errors)
    ol_mean, ol_std = mean_curve(ol_errors)

    ax.fill_between(leads, da_mean - da_std, da_mean + da_std,
                    color="tab:blue", alpha=0.15, zorder=2)
    ax.fill_between(leads, ol_mean - ol_std, ol_mean + ol_std,
                    color="tab:gray",  alpha=0.15, zorder=2)
    ax.plot(leads, da_mean, color="tab:blue", lw=2.6, marker="o",
            zorder=5, label="DA — mean error (±1 std shaded)")
    ax.plot(leads, ol_mean, color="tab:gray",  lw=2.6, marker="s",
            linestyle="--", zorder=5, label="Open-loop — mean error (±1 std shaded)")

    ax.set_xlabel("Forecast lead hour (hours before target verification time)", fontsize=11)
    ax.set_ylabel("Mean error: forecast − USGS obs (m³/s)", fontsize=11)
    ax.set_xticks(leads)
    ax.set_title(
        "Mean forecast error vs lead time — fixed verification time, Helene peak window\n"
        f"F1 (variance-scaled Vrugt) | Target: {TARGET_START.strftime('%b %d %H UTC')} → "
        f"{TARGET_END.strftime('%b %d %H UTC')}  |  Lead 1 = init 1 hr before target",
        fontsize=11,
    )
    ax.legend(fontsize=10, loc="lower right", frameon=True, framealpha=0.92)
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    out2 = os.path.join(out_dir, "error_fixed_target_mean.png")
    plt.savefig(out2, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out2}")


if __name__ == "__main__":
    main()
