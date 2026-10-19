"""
plot_da_perturbation_arms.py  —  2a/2b: ensemble spread WITH DA on.

Loads the two arm CSVs produced by run_perturbation_da_on.py:
    <arm-dir>/<cat-id>/<cat-id>_da_forcing_arm.csv   (30 members)
    <arm-dir>/<cat-id>/<cat-id>_da_hydro_arm.csv     (20 members)

Plots:
    Panel 2a — Forcing arm: 30-member spaghetti over Helene window
               (spread = met forcing uncertainty with DA-corrected initial states)
    Panel 2b — Hydro-state arm: 20-member spaghetti over Helene window
               (spread = initial state uncertainty with deterministic forcing)
               Optional: --ol-csv overlays the open-loop grand median as a
               thick dashed gray line for comparison.
    Panel 2c — Comparison: median ± spread envelope, both arms + USGS obs

All trajectories projected to valid_time = issue_time + lead_hour hours.
Colored by initialization date (Sep 24-30). USGS obs in black.

Outputs:
    <out-dir>/<cat-id>/<cat-id>_2a_forcing_arm_helene.png
    <out-dir>/<cat-id>/<cat-id>_2b_hydro_arm_helene.png
    <out-dir>/<cat-id>/<cat-id>_2ab_arms_comparison.png

Run on server (troute env):
    python3 plot_da_perturbation_arms.py \\
        --cat-id cat-1016300 \\
        --arm-dir /mnt/disk2/suma_helen_poster/da_results/da_arms_f5_rekrig_20pct \\
        --ol-csv  /mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct/routed/routed_leadtime_openloop_full.parquet
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

DEFAULT_ARM_DIR  = "/mnt/disk2/suma_helen_poster/da_results/da_arms_f5_rekrig_20pct"
DEFAULT_USGS_CSV = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
DEFAULT_CAT_ID   = "cat-1016300"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S       = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

PLOT_START        = pd.Timestamp("2024-09-24 00:00:00")
PLOT_END          = pd.Timestamp("2024-09-30 06:00:00")
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


def load_arm(path):
    df = pd.read_csv(path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    df["valid_time"] = df["issue_time"] + pd.to_timedelta(df["lead_hour"], unit="h")
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    return df, member_cols


def load_openloop(path):
    """Load open-loop lead-time file (CSV or Parquet) and return grand median
    indexed by valid_time.

    Accepts:
      - routed_leadtime_openloop_full.parquet  (T-route output, m³/s, recommended)
      - cat-*_lead_time_forecasts_openloop.csv (unrouted CFE mm/h, single catchment)

    Returns a pd.Series (valid_time → median q in m³/s) clipped to plot window.
    """
    if path.endswith(".parquet"):
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)

    df["issue_time"] = pd.to_datetime(df["issue_time"])
    if "valid_time" in df.columns:
        df["valid_time"] = pd.to_datetime(df["valid_time"])
    elif "lead_hour" in df.columns:
        df["valid_time"] = (df["issue_time"]
                            + pd.to_timedelta(df["lead_hour"], unit="h"))
    else:
        raise ValueError("Open-loop file must have valid_time or lead_hour column")

    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    mask = (df["valid_time"] >= PLOT_START) & (df["valid_time"] <= PLOT_END)
    df = df[mask].copy()
    if df.empty:
        return pd.Series(dtype=float)

    vals = df[member_cols].to_numpy(dtype=float)
    # Only convert if values are clearly in mm/h (unrouted CFE output).
    # Routed parquet is already in m³/s — do not convert.
    if path.endswith(".csv") and np.nanmedian(vals[vals > 0]) < 5:
        vals = vals * MM_H_TO_M3_S

    df["q_grand_median"] = np.nanmedian(vals, axis=1)
    series = (df.groupby("valid_time")["q_grand_median"]
                .median()
                .sort_index())
    return series


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


def _add_helene_band(ax):
    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color="salmon", alpha=0.12, zorder=0)
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color="firebrick", ha="left", va="top")


def _format_xaxis(ax):
    ax.set_xlim(PLOT_START, PLOT_END)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    ax.grid(True, alpha=0.2, lw=0.4)


def plot_spaghetti_arm(ax, df, member_cols, arm_color, obs_series,
                       title, label_stem, lw_thin=0.7, alpha_thin=0.35):
    _add_helene_band(ax)

    issue_times = sorted(df["issue_time"].unique())
    all_means = []

    for t0 in issue_times:
        date_str = str(pd.Timestamp(t0).date())
        color = DATE_COLORS.get(date_str, arm_color)
        sub = df[df["issue_time"] == t0].sort_values("valid_time")
        mask = (sub["valid_time"] >= PLOT_START) & (sub["valid_time"] <= PLOT_END)
        sub_w = sub[mask]
        if sub_w.empty:
            continue
        vt_w = sub_w["valid_time"].values

        mem_vals = sub_w[member_cols].to_numpy(dtype=float)
        if "mm" not in label_stem.lower():
            mem_vals = mem_vals * MM_H_TO_M3_S

        ax.fill_between(vt_w,
                        np.nanmin(mem_vals, axis=1),
                        np.nanmax(mem_vals, axis=1),
                        color=color, alpha=0.06, zorder=2)
        ax.plot(vt_w, np.nanmedian(mem_vals, axis=1),
                color=color, lw=lw_thin, alpha=alpha_thin + 0.1, zorder=3)

        all_means.append(
            pd.Series(np.nanmedian(mem_vals, axis=1), index=vt_w))

    if all_means:
        full_idx = pd.date_range(PLOT_START, PLOT_END, freq="1h")
        stacked = pd.concat(all_means, axis=1).reindex(full_idx)
        grand_mean = stacked.mean(axis=1)
        ax.plot(grand_mean.index, grand_mean.values,
                color=arm_color, lw=2.4, alpha=0.95, zorder=5,
                label=f"{label_stem} — mean across all forecasts")

    obs_w = obs_series.loc[PLOT_START:PLOT_END]
    ax.plot(obs_w.index, obs_w.values,
            color="black", lw=1.8, zorder=6, label="USGS obs")

    ax.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax.set_title(title, fontsize=11)
    _format_xaxis(ax)


def compute_envelope(df, member_cols):
    records = []
    for t0, grp in df.groupby("issue_time"):
        mask = (grp["valid_time"] >= PLOT_START) & (grp["valid_time"] <= PLOT_END)
        sub = grp[mask]
        if sub.empty:
            continue
        vals = sub[member_cols].to_numpy(dtype=float) * MM_H_TO_M3_S
        for i, row in enumerate(sub.itertuples()):
            records.append({
                "valid_time": row.valid_time,
                "q_min":  np.nanmin(vals[i]),
                "q_med":  np.nanmedian(vals[i]),
                "q_max":  np.nanmax(vals[i]),
            })
    if not records:
        return pd.DataFrame(columns=["valid_time", "q_min", "q_med", "q_max"])

    env_df = pd.DataFrame(records)
    env_df = (env_df.groupby("valid_time")
              .agg(q_min=("q_min", "min"),
                   q_med=("q_med", "mean"),
                   q_max=("q_max", "max"))
              .reset_index()
              .sort_values("valid_time"))
    return env_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm-dir",  default=DEFAULT_ARM_DIR)
    parser.add_argument("--cat-id",   default=DEFAULT_CAT_ID)
    parser.add_argument("--usgs-csv", default=DEFAULT_USGS_CSV)
    parser.add_argument("--out-dir",  default=None)
    parser.add_argument(
        "--ol-csv", default=None,
        help=(
            "Path to open-loop lead-time forecast (CSV or Parquet). "
            "When provided, the grand median is overlaid on panel 2b "
            "as a thick dashed black line."
        ),
    )
    args = parser.parse_args()

    cat_dir = os.path.join(args.arm_dir, args.cat_id)
    out_dir = args.out_dir or cat_dir
    os.makedirs(out_dir, exist_ok=True)

    forcing_path = os.path.join(cat_dir, f"{args.cat_id}_da_forcing_arm.csv")
    hydro_path   = os.path.join(cat_dir, f"{args.cat_id}_da_hydro_arm.csv")

    print(f"Loading arm CSVs for {args.cat_id}...")
    df_fa, fa_cols = load_arm(forcing_path)
    df_ha, ha_cols = load_arm(hydro_path)
    obs = load_usgs(args.usgs_csv)
    print(f"  Forcing arm: {df_fa['issue_time'].nunique()} issue times, "
          f"{len(fa_cols)} members")
    print(f"  Hydro arm:   {df_ha['issue_time'].nunique()} issue times, "
          f"{len(ha_cols)} members")

    # ── Figure 2a: Forcing arm ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(17, 6))
    plot_spaghetti_arm(
        ax, df_fa, fa_cols,
        arm_color="tab:blue", obs_series=obs,
        title=(f"2a — Forcing arm (DA on, {len(fa_cols)} members): "
               "met forcing uncertainty | Sep 24-30 2024\n"
               f"{args.cat_id}  |  Shaded = member min-max  |  "
               "Line = median per init time  |  Thick = grand mean"),
        label_stem="Forcing arm",
    )
    patches = [mpatches.Patch(color=c, label=f"Init {d}")
               for d, c in DATE_COLORS.items()]
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + patches,
              fontsize=8, loc="upper left", frameon=True, framealpha=0.9, ncol=2)
    plt.tight_layout()
    out_a = os.path.join(out_dir, f"{args.cat_id}_2a_forcing_arm_helene.png")
    plt.savefig(out_a, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_a}")

    # ── Figure 2b: Hydro-state arm ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(17, 6))
    plot_spaghetti_arm(
        ax, df_ha, ha_cols,
        arm_color="tab:green", obs_series=obs,
        title=(f"2b — Hydro-state arm (DA on, {len(ha_cols)} members): "
               "initial state uncertainty | Sep 24-30 2024\n"
               f"{args.cat_id}  |  Shaded = member min-max  |  "
               "Line = median per init time  |  Thick = grand mean"),
        label_stem="Hydro-state arm",
    )
    if args.ol_csv:
        print(f"Loading open-loop: {args.ol_csv}")
        ol_series = load_openloop(args.ol_csv)
        if not ol_series.empty:
            ax.plot(
                ol_series.index, ol_series.values,
                color="black", lw=2.5, ls="--", alpha=0.95, zorder=7,
                label="Open loop (no DA) — grand median",
            )
        else:
            print("  Warning: open-loop yielded no data in plot window.")

    patches = [mpatches.Patch(color=c, label=f"Init {d}")
               for d, c in DATE_COLORS.items()]
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + patches,
              fontsize=8, loc="upper left", frameon=True, framealpha=0.9, ncol=2)
    plt.tight_layout()
    out_b = os.path.join(out_dir, f"{args.cat_id}_2b_hydro_arm_helene.png")
    plt.savefig(out_b, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_b}")

    # ── Figure 2ab: Comparison ────────────────────────────────────────────────
    print("Computing aggregated envelopes for comparison plot...")
    env_fa = compute_envelope(df_fa, fa_cols)
    env_ha = compute_envelope(df_ha, ha_cols)

    fig, ax = plt.subplots(figsize=(17, 7))
    _add_helene_band(ax)

    if not env_fa.empty:
        ax.fill_between(env_fa["valid_time"], env_fa["q_min"], env_fa["q_max"],
                        color="tab:blue", alpha=0.18, zorder=2,
                        label=f"Forcing arm spread (N={len(fa_cols)} members)")
        ax.plot(env_fa["valid_time"], env_fa["q_med"],
                color="tab:blue", lw=2.0, zorder=4,
                label="Forcing arm — grand median")

    if not env_ha.empty:
        ax.fill_between(env_ha["valid_time"], env_ha["q_min"], env_ha["q_max"],
                        color="tab:green", alpha=0.18, zorder=2,
                        label=f"Hydro-state arm spread (N={len(ha_cols)} members)")
        ax.plot(env_ha["valid_time"], env_ha["q_med"],
                color="tab:green", lw=2.0, zorder=4,
                label="Hydro-state arm — grand median")

    obs_w = obs.loc[PLOT_START:PLOT_END]
    ax.plot(obs_w.index, obs_w.values,
            color="black", lw=1.8, zorder=6, label="USGS obs")

    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_xlabel("Date", fontsize=11)
    ax.set_title(
        f"2a vs 2b — Forcing arm (blue) vs Hydro-state arm (green) | "
        f"DA on | {args.cat_id}\n"
        "Shaded = full member spread (min-max aggregated across all issue times). "
        "Lines = grand median.",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper left", frameon=True, framealpha=0.9)
    _format_xaxis(ax)
    plt.tight_layout()
    out_c = os.path.join(out_dir, f"{args.cat_id}_2ab_arms_comparison.png")
    plt.savefig(out_c, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_c}")


if __name__ == "__main__":
    main()
