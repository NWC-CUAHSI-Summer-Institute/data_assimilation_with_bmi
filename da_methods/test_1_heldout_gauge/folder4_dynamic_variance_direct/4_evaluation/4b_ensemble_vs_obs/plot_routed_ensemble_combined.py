#!/usr/bin/env python3
"""
plot_routed_ensemble_combined.py
Combined poster plot: 600-member routed ensemble envelope +
deterministic DA lines (Vrugt/no-Vrugt) vs USGS gauge obs at 03463300.

Inputs:
  --vrugt-csv    : routed_Q_test.csv from vrugt_dynamic_routed/
  --novrugt-csv  : routed_Q_test.csv from novrugt_dynamic_routed/
  --ensemble-pq  : routed_crossed_ensemble.parquet from run_route_crossed_ensemble.py
  --out-dir      : output directory

The ensemble parquet columns: issue_time, lead_hour, member_0000..member_0599 (Q in m3/s)

Usage:
    python3 plot_routed_ensemble_combined.py \
        --vrugt-csv   /mnt/disk2/suma_helen_poster/da_results/vrugt_dynamic_routed/routed_Q_test.csv \
        --novrugt-csv /mnt/disk2/suma_helen_poster/da_results/novrugt_dynamic_routed/routed_Q_test.csv \
        --ensemble-pq /mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble_routed/routed_crossed_ensemble.parquet \
        --out-dir     /mnt/disk2/suma_helen_poster/da_results/comparison_plots
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

HELENE_START = pd.Timestamp("2024-09-25 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-29 00:00:00")
ZOOM_START   = pd.Timestamp("2024-09-24 00:00:00")
ZOOM_END     = pd.Timestamp("2024-09-28 18:00:00")  # cut trailing ensemble bump

COLOR_OBS     = "black"
COLOR_VRUGT   = "#1f77b4"   # blue
COLOR_NOVRUGT = "#d62728"   # red
COLOR_ENS     = "#1f77b4"   # same blue family as Vrugt


def load_det_csv(path):
    df = pd.read_csv(path, parse_dates=["date"])
    return df.set_index("date").sort_index()


def build_envelope(pq_path, zoom_start, zoom_end):
    """
    Load routed ensemble parquet and build percentile envelope by valid_time.
    Returns DataFrame indexed by valid_time with columns p05, p25, p50, p75, p95.
    """
    df = pd.read_parquet(pq_path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    member_cols = [c for c in df.columns if c.startswith("member_")]

    records = []
    for _, row in df.iterrows():
        t0      = row["issue_time"]
        lead    = int(row["lead_hour"])
        vt      = t0 + pd.Timedelta(hours=lead)
        vals    = row[member_cols].values.astype(float)
        records.append({"valid_time": vt, "vals": vals})

    # Group by valid_time (floor to hour), stack all member values across issue_times
    from collections import defaultdict
    groups = defaultdict(list)
    for rec in records:
        vt = rec["valid_time"].floor("1h")
        groups[vt].extend(rec["vals"].tolist())

    rows = []
    for vt in sorted(groups):
        if vt < zoom_start or vt > zoom_end:
            continue
        v = np.array(groups[vt])
        v = v[~np.isnan(v)]
        if len(v) == 0:
            continue
        rows.append({
            "valid_time": vt,
            "p05": np.percentile(v, 5),
            "p25": np.percentile(v, 25),
            "p50": np.percentile(v, 50),
            "p75": np.percentile(v, 75),
            "p95": np.percentile(v, 95),
        })

    return pd.DataFrame(rows).set_index("valid_time")


def compute_kge(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.std(o) == 0:
        return np.nan
    r = np.corrcoef(o, s)[0, 1]
    return 1.0 - np.sqrt((r-1)**2 + (np.std(s)/np.std(o)-1)**2 + (np.mean(s)/np.mean(o)-1)**2)


def _format_xaxis(ax, locator, fmt, rotation=30):
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=rotation, ha="right", fontsize=9)


def _add_helene_band(ax):
    ax.axvspan(HELENE_START, HELENE_END, color="gold", alpha=0.15, zorder=0)


def plot_zoom(ax, dates, obs, sim_v, sim_nv, kge_v, kge_nv, env=None):
    mask = (dates >= ZOOM_START) & (dates <= ZOOM_END)
    dz   = dates[mask]

    # Ensemble envelope (bottom layer)
    if env is not None:
        ax.fill_between(env.index, env["p05"], env["p95"],
                        color=COLOR_ENS, alpha=0.18, label="Ensemble 5-95th pct")
        ax.fill_between(env.index, env["p25"], env["p75"],
                        color=COLOR_ENS, alpha=0.30, label="Ensemble 25-75th pct")
        ax.plot(env.index, env["p50"],
                color=COLOR_ENS, lw=1.2, linestyle="--", alpha=0.7, label="Ensemble median")

    # Deterministic DA lines
    ax.plot(dz, sim_nv[mask], color=COLOR_NOVRUGT, lw=1.4, linestyle="--",
            label=f"DA -- constant R  KGE={kge_nv:.3f}", zorder=3)
    ax.plot(dz, sim_v[mask],  color=COLOR_VRUGT,   lw=1.6,
            label=f"DA -- dynamic R (Vrugt)  KGE={kge_v:.3f}", zorder=4)

    # USGS obs on top
    ax.plot(dz, obs[mask], color=COLOR_OBS, lw=1.6, label="USGS obs (gauge 03463300)", zorder=5)

    peak_usgs = np.nanmax(obs[mask])
    ax.axhline(peak_usgs, color=COLOR_OBS, lw=0.7, linestyle=":", alpha=0.5)
    # Place annotation to the right of the peak to avoid legend overlap
    ax.text(pd.Timestamp("2024-09-28 06:00:00"), peak_usgs * 1.02,
            f"USGS peak {peak_usgs:.0f} m³/s", fontsize=8.5, color="black",
            alpha=0.8, ha="right")

    _add_helene_band(ax)
    _format_xaxis(ax, mdates.DayLocator(interval=1), "%b %d")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vrugt-csv",   required=True)
    parser.add_argument("--novrugt-csv", required=True)
    parser.add_argument("--ensemble-pq", required=True,
                        help="routed_crossed_ensemble.parquet from run_route_crossed_ensemble.py")
    parser.add_argument("--out-dir",     required=True)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    vrugt   = load_det_csv(args.vrugt_csv)
    novrugt = load_det_csv(args.novrugt_csv)

    obs   = vrugt["Q_usgs_m3s"].values
    dates = vrugt.index
    sim_v  = vrugt["Q_routed_m3s"].values
    sim_nv = novrugt["Q_routed_m3s"].reindex(dates).values

    kge_v  = compute_kge(obs, sim_v)
    kge_nv = compute_kge(obs, sim_nv)

    print("Building ensemble envelope...")
    env = build_envelope(args.ensemble_pq, ZOOM_START, ZOOM_END)
    print(f"  Envelope: {len(env)} timesteps, "
          f"peak p95={env['p95'].max():.1f} m3/s, peak p50={env['p50'].max():.1f} m3/s")

    peak_usgs = np.nanmax(obs[(dates >= ZOOM_START) & (dates <= ZOOM_END)])
    in_env = env["p95"].max() >= peak_usgs
    print(f"  USGS peak {peak_usgs:.1f} m3/s {'IS' if in_env else 'IS NOT'} within p95 envelope")

    # ------------------------------------------------------------------ #
    # Figure 1: Helene zoom only -- main poster panel
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(11, 5.5))
    plot_zoom(ax, dates, obs, sim_v, sim_nv, kge_v, kge_nv, env)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_title(
        "Hurricane Helene -- Probabilistic discharge forecast at gauge 03463300 (South Toe River)\n"
        "600-member CFE ensemble + DA (Muskingum routing)",
        fontsize=11)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    out1 = os.path.join(args.out_dir, "helene_ensemble_vs_usgs.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out1}")

    # ------------------------------------------------------------------ #
    # Figure 2: two-panel (full period + Helene zoom)
    # ------------------------------------------------------------------ #
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(14, 9),
        gridspec_kw={"height_ratios": [1, 1.5]})

    # Top: full test period (deterministic only, no ensemble for readability)
    ax_top.plot(dates, obs,    color=COLOR_OBS,     lw=0.8, label="USGS obs", zorder=3)
    ax_top.plot(dates, sim_v,  color=COLOR_VRUGT,   lw=0.9,
                label=f"DA dynamic R (Vrugt)  KGE={kge_v:.3f}", zorder=2)
    ax_top.plot(dates, sim_nv, color=COLOR_NOVRUGT, lw=0.9, linestyle="--",
                label=f"DA constant R          KGE={kge_nv:.3f}", zorder=2)
    ax_top.axvspan(ZOOM_START, ZOOM_END, color="gray", alpha=0.12, zorder=0)
    ax_top.axvspan(HELENE_START, HELENE_END, color="gold", alpha=0.15, zorder=0)
    ax_top.text(HELENE_START + pd.Timedelta(days=0.5),
                np.nanmax(obs) * 0.88, "Helene",
                fontsize=9, color="goldenrod", fontweight="bold")
    ax_top.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax_top.set_title(
        "Routed discharge at gauge 03463300 (South Toe River) -- Test period\n"
        "DA with dynamic Vrugt R vs. constant R -- Muskingum routing",
        fontsize=11)
    ax_top.legend(fontsize=9, loc="upper left")
    ax_top.grid(True, alpha=0.25)
    _format_xaxis(ax_top, mdates.MonthLocator(interval=2), "%Y-%m")

    # Bottom: Helene zoom with ensemble
    plot_zoom(ax_bot, dates, obs, sim_v, sim_nv, kge_v, kge_nv, env)
    ax_bot.set_xlabel("Date (UTC)", fontsize=10)
    ax_bot.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax_bot.set_title("Hurricane Helene window -- 600-member ensemble envelope", fontsize=10)
    ax_bot.legend(fontsize=9, loc="lower right")
    ax_bot.grid(True, alpha=0.25)

    plt.tight_layout()
    out2 = os.path.join(args.out_dir, "helene_ensemble_twopanel.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out2}")


if __name__ == "__main__":
    main()
