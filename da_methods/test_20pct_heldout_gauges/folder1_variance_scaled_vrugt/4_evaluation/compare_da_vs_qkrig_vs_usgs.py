"""
compare_da_vs_qkrig_vs_usgs.py

Three-way comparison at gauge 03463300:
    DA-routed Q  vs  routed Qkrig  vs  USGS obs

Reads the routed_Q_test.csv produced by run_route.py (which already contains
Q_routed_m3s, Q_usgs_m3s, and Q_krig_m3s columns) and prints KGE/NSE/peak
for the full test period and the Helene window separately.

Also saves a two-panel comparison figure.

Usage:
    python3 compare_da_vs_qkrig_vs_usgs.py \
        --vrugt-csv /mnt/disk2/suma_helen_poster/da_results/vrugt_dynamic_routed/routed_Q_test.csv \
        --out-dir   /mnt/disk2/suma_helen_poster/da_results/comparison_plots
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-29 23:00:00")


def kge(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.std(o) == 0:
        return np.nan
    r = np.corrcoef(o, s)[0, 1]
    return 1.0 - np.sqrt((r-1)**2 + (np.std(s)/np.std(o)-1)**2 + (np.mean(s)/np.mean(o)-1)**2)


def nse(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    o, s = obs[mask], sim[mask]
    if len(o) < 2:
        return np.nan
    denom = np.sum((o - o.mean())**2)
    return 1.0 - np.sum((o - s)**2) / denom if denom > 0 else np.nan


def peak_ratio(obs, sim):
    return np.nanmax(sim) / np.nanmax(obs)


def print_stats(label, obs, sim, dates, window_name="full period"):
    mask = np.isfinite(obs) & np.isfinite(sim)
    print(f"  {label:35s}  KGE={kge(obs,sim):+.3f}  NSE={nse(obs,sim):+.3f}  "
          f"peak_sim={np.nanmax(sim):.1f}  peak_obs={np.nanmax(obs):.1f}  "
          f"ratio={peak_ratio(obs,sim):.2f}x  [{window_name}]")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vrugt-csv", required=True,
                        help="routed_Q_test.csv from vrugt_dynamic_routed/")
    parser.add_argument("--out-dir",   required=True)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.vrugt_csv, parse_dates=["date"]).set_index("date").sort_index()

    if "Q_usgs_m3s" not in df.columns:
        raise ValueError("CSV missing Q_usgs_m3s — run run_route.py with --obs-csv")
    if "Q_krig_m3s" not in df.columns:
        raise ValueError("CSV missing Q_krig_m3s — need run_route.py with --kv-dir or --obs-csv")

    obs  = df["Q_usgs_m3s"].values
    da   = df["Q_routed_m3s"].values
    krig = df["Q_krig_m3s"].values
    dates = df.index

    helene = (dates >= HELENE_START) & (dates <= HELENE_END)

    print("\n" + "="*85)
    print("THREE-WAY COMPARISON  —  gauge 03463300 (South Toe River)")
    print("="*85)

    print("\n[FULL TEST PERIOD]")
    print_stats("DA-routed   vs USGS", obs, da,   dates, "full")
    print_stats("Qkrig-routed vs USGS", obs, krig, dates, "full")
    print_stats("DA-routed   vs Qkrig-routed", krig, da, dates, "full")

    print("\n[HELENE WINDOW  Sep 24–29]")
    print_stats("DA-routed   vs USGS", obs[helene], da[helene],   dates[helene], "Helene")
    print_stats("Qkrig-routed vs USGS", obs[helene], krig[helene], dates[helene], "Helene")
    print_stats("DA-routed   vs Qkrig-routed", krig[helene], da[helene], dates[helene], "Helene")

    print("\n[PEAK VALUES  (Helene window)]")
    print(f"  USGS peak        : {np.nanmax(obs[helene]):.1f} m³/s")
    print(f"  DA-routed peak   : {np.nanmax(da[helene]):.1f} m³/s  "
          f"({np.nanmax(da[helene])/np.nanmax(obs[helene])*100:.0f}% of USGS)")
    print(f"  Qkrig-routed peak: {np.nanmax(krig[helene]):.1f} m³/s  "
          f"({np.nanmax(krig[helene])/np.nanmax(obs[helene])*100:.0f}% of USGS)")
    print(f"  DA vs Qkrig peak : DA is "
          f"{'higher' if np.nanmax(da[helene]) > np.nanmax(krig[helene]) else 'lower'} "
          f"by {abs(np.nanmax(da[helene])-np.nanmax(krig[helene])):.1f} m³/s")
    print("="*85 + "\n")

    # ---- Plot ----
    fig, axes = plt.subplots(2, 1, figsize=(15, 9),
                             gridspec_kw={"height_ratios": [1, 1.6]})

    # Top: full period
    ax = axes[0]
    ax.plot(dates, obs,  color="black",      lw=0.8, label="USGS obs", zorder=4)
    ax.plot(dates, krig, color="#e6820e",    lw=0.8, linestyle="--",
            label=f"Qkrig-routed  KGE={kge(obs,da):.3f}", zorder=2)
    ax.plot(dates, da,   color="#1f77b4",    lw=0.9,
            label=f"DA-routed  KGE={kge(obs,da):.3f}", zorder=3)
    ax.axvspan(HELENE_START, HELENE_END, color="gold", alpha=0.15, zorder=0)
    ax.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax.set_title("Full test period — DA-routed vs Qkrig-routed vs USGS obs", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right", fontsize=8)

    # Bottom: Helene zoom
    ax = axes[1]
    dh = dates[helene]
    kge_da_h   = kge(obs[helene], da[helene])
    kge_krig_h = kge(obs[helene], krig[helene])
    nse_da_h   = nse(obs[helene], da[helene])
    nse_krig_h = nse(obs[helene], krig[helene])

    ax.plot(dh, obs[helene],  color="black",   lw=1.8, zorder=4, label="USGS obs")
    ax.plot(dh, krig[helene], color="#e6820e", lw=1.4, linestyle="--", zorder=2,
            label=f"Qkrig-routed  KGE={kge_krig_h:.3f}  NSE={nse_krig_h:.3f}")
    ax.plot(dh, da[helene],   color="#1f77b4", lw=1.6, zorder=3,
            label=f"DA-routed      KGE={kge_da_h:.3f}  NSE={nse_da_h:.3f}")

    peak_usgs = np.nanmax(obs[helene])
    ax.axhline(peak_usgs, color="black", lw=0.6, linestyle=":", alpha=0.5)
    ax.text(HELENE_END - pd.Timedelta(hours=6), peak_usgs * 1.02,
            f"USGS peak {peak_usgs:.0f} m³/s", fontsize=8.5, ha="right", color="black")

    ax.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax.set_xlabel("Date (UTC)", fontsize=10)
    ax.set_title("Helene window — does DA add value beyond Qkrig?", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right", fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(args.out_dir, "da_vs_qkrig_vs_usgs.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
