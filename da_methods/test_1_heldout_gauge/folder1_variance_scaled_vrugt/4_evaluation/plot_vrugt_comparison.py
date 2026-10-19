#!/usr/bin/env python3
"""
plot_vrugt_comparison.py
Compare routed discharge: dynamic Vrugt R vs. constant R (no-Vrugt)
against USGS gauge obs at 03463300 (South Toe River Near Celo, NC).

Reads routed_Q_test.csv from both run directories (produced by run_route.py
with --obs-csv). Produces:
  - Two-panel figure: full test period + Helene zoom
  - Single-panel Helene zoom only

Usage:
    python3 plot_vrugt_comparison.py \
        --vrugt-csv  /mnt/disk2/suma_helen_poster/da_results/vrugt_dynamic_routed/routed_Q_test.csv \
        --novrugt-csv /mnt/disk2/suma_helen_poster/da_results/novrugt_dynamic_routed/routed_Q_test.csv \
        --out-dir /mnt/disk2/suma_helen_poster/da_results/comparison_plots
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

HELENE_START = pd.Timestamp("2024-09-25 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-29 00:00:00")
ZOOM_START   = pd.Timestamp("2024-09-24 00:00:00")
ZOOM_END     = pd.Timestamp("2024-09-30 00:00:00")

COLOR_OBS    = "black"
COLOR_VRUGT  = "#1f77b4"   # blue
COLOR_NOVRUGT = "#d62728"  # red


def load_csv(path):
    df = pd.read_csv(path, parse_dates=["date"])
    df = df.set_index("date").sort_index()
    return df


def compute_kge(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.std(o) == 0:
        return np.nan
    r = np.corrcoef(o, s)[0, 1]
    return 1.0 - np.sqrt((r - 1)**2 + (np.std(s) / np.std(o) - 1)**2 + (np.mean(s) / np.mean(o) - 1)**2)


def compute_nse(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.sum((o - np.mean(o))**2) == 0:
        return np.nan
    return 1.0 - np.sum((o - s)**2) / np.sum((o - np.mean(o))**2)


def _add_helene_band(ax):
    ax.axvspan(HELENE_START, HELENE_END,
               color="gold", alpha=0.18, zorder=0, label="_nolegend_")


def _format_xaxis(ax, locator, fmt, rotation=30):
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=rotation, ha="right", fontsize=9)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vrugt-csv",   required=True)
    parser.add_argument("--novrugt-csv", required=True)
    parser.add_argument("--out-dir",     required=True)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    vrugt   = load_csv(args.vrugt_csv)
    novrugt = load_csv(args.novrugt_csv)

    # USGS obs from vrugt CSV (same timestamps, both runs used same --obs-csv)
    obs = vrugt["Q_usgs_m3s"].values
    dates = vrugt.index

    sim_v  = vrugt["Q_routed_m3s"].values
    sim_nv = novrugt["Q_routed_m3s"].reindex(dates).values

    kge_v  = compute_kge(obs, sim_v)
    nse_v  = compute_nse(obs, sim_v)
    kge_nv = compute_kge(obs, sim_nv)
    nse_nv = compute_nse(obs, sim_nv)

    peak_usgs = np.nanmax(obs)
    peak_v    = np.nanmax(sim_v)
    peak_nv   = np.nanmax(sim_nv)

    print(f"USGS peak:      {peak_usgs:.1f} m3/s")
    print(f"Vrugt peak:     {peak_v:.1f} m3/s   KGE={kge_v:.3f}  NSE={nse_v:.3f}")
    print(f"No-Vrugt peak:  {peak_nv:.1f} m3/s   KGE={kge_nv:.3f}  NSE={nse_nv:.3f}")

    # ------------------------------------------------------------------ #
    # Figure 1: two-panel (full period + Helene zoom)
    # ------------------------------------------------------------------ #
    fig, (ax_full, ax_zoom) = plt.subplots(
        2, 1, figsize=(14, 9),
        gridspec_kw={"height_ratios": [1, 1.4]})

    # -- Top: full period --
    ax_full.plot(dates, obs,    color=COLOR_OBS,     lw=0.8, label="USGS obs (gauge 03463300)", zorder=3)
    ax_full.plot(dates, sim_v,  color=COLOR_VRUGT,   lw=0.9, label=f"Vrugt (dynamic R)   KGE={kge_v:.3f}", zorder=2)
    ax_full.plot(dates, sim_nv, color=COLOR_NOVRUGT, lw=0.9, linestyle="--",
                 label=f"No-Vrugt (const R)  KGE={kge_nv:.3f}", zorder=2)
    _add_helene_band(ax_full)
    ax_full.text(HELENE_START + pd.Timedelta(hours=12), ax_full.get_ylim()[1] * 0.85,
                 "Helene", fontsize=9, color="goldenrod", fontweight="bold")
    ax_full.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax_full.set_title(
        "Routed discharge at gauge 03463300 (South Toe River)\n"
        "DA with dynamic Vrugt R vs. constant R -- Muskingum routing",
        fontsize=11)
    ax_full.legend(fontsize=9, loc="upper left")
    ax_full.grid(True, alpha=0.25)
    _format_xaxis(ax_full, mdates.MonthLocator(interval=2), "%Y-%m")

    # Add zoom indicator
    ax_full.axvspan(ZOOM_START, ZOOM_END, color="gray", alpha=0.10, zorder=0)
    ax_full.annotate("", xy=(ZOOM_END, ax_full.get_ylim()[1] * 0.5),
                     xytext=(ZOOM_START, ax_full.get_ylim()[1] * 0.5),
                     arrowprops=dict(arrowstyle="<->", color="gray", lw=1.0))

    # -- Bottom: Helene zoom --
    mask_zoom = (dates >= ZOOM_START) & (dates <= ZOOM_END)
    dz = dates[mask_zoom]
    ax_zoom.plot(dz, obs[mask_zoom],    color=COLOR_OBS,     lw=1.2, label="USGS obs (gauge 03463300)", zorder=3)
    ax_zoom.plot(dz, sim_v[mask_zoom],  color=COLOR_VRUGT,   lw=1.4,
                 label=f"Vrugt (dynamic R)   KGE={kge_v:.3f}  NSE={nse_v:.3f}  peak={peak_v:.0f} m³/s", zorder=2)
    ax_zoom.plot(dz, sim_nv[mask_zoom], color=COLOR_NOVRUGT, lw=1.4, linestyle="--",
                 label=f"No-Vrugt (const R)  KGE={kge_nv:.3f}  NSE={nse_nv:.3f}  peak={peak_nv:.0f} m³/s", zorder=2)
    _add_helene_band(ax_zoom)
    ax_zoom.axhline(peak_usgs, color=COLOR_OBS, lw=0.7, linestyle=":", alpha=0.6)
    ax_zoom.text(ZOOM_START + pd.Timedelta(hours=3), peak_usgs * 1.01,
                 f"USGS peak {peak_usgs:.0f} m³/s", fontsize=8, color="black", alpha=0.8)
    ax_zoom.set_xlabel("Date (UTC)", fontsize=10)
    ax_zoom.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax_zoom.set_title("Hurricane Helene window (Sep 24-30, 2024)", fontsize=10)
    ax_zoom.legend(fontsize=9, loc="lower right")
    ax_zoom.grid(True, alpha=0.25)
    _format_xaxis(ax_zoom, mdates.DayLocator(interval=1), "%b %d")

    plt.tight_layout()
    out1 = os.path.join(args.out_dir, "vrugt_vs_novrugt_twopanel.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out1}")

    # ------------------------------------------------------------------ #
    # Figure 2: Helene zoom only (poster-ready single panel)
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(dz, obs[mask_zoom],    color=COLOR_OBS,     lw=1.4, label="USGS obs (gauge 03463300)", zorder=3)
    ax.plot(dz, sim_v[mask_zoom],  color=COLOR_VRUGT,   lw=1.6,
            label=f"DA -- dynamic R (Vrugt)   KGE={kge_v:.3f}", zorder=2)
    ax.plot(dz, sim_nv[mask_zoom], color=COLOR_NOVRUGT, lw=1.6, linestyle="--",
            label=f"DA -- constant R           KGE={kge_nv:.3f}", zorder=2)
    _add_helene_band(ax)
    ax.axhline(peak_usgs, color=COLOR_OBS, lw=0.8, linestyle=":", alpha=0.5)
    ax.text(ZOOM_START + pd.Timedelta(hours=3), peak_usgs * 1.02,
            f"USGS peak {peak_usgs:.0f} m³/s", fontsize=9, color="black", alpha=0.75)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_title(
        "Hurricane Helene -- Routed discharge at gauge 03463300 (South Toe River)\n"
        "CFE + DA (Muskingum routing): dynamic vs. constant observation error variance",
        fontsize=11)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(True, alpha=0.25)
    _format_xaxis(ax, mdates.DayLocator(interval=1), "%b %d")
    plt.tight_layout()
    out2 = os.path.join(args.out_dir, "vrugt_vs_novrugt_helene_zoom.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out2}")


if __name__ == "__main__":
    main()
