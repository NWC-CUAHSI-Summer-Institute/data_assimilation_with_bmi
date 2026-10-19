"""
plot_routed_gauge_check_f5.py

Gauge-level evaluation for F5 (re-kriged variance, 20% holdout).

Reads routed_Q_test.csv (T-Route output from Step 2) and plots
DA-routed streamflow vs USGS observed discharge at gauge 03463300.
Optionally overlays the open-loop routed timeseries (--ol-csv).
Units are m³/s — directly comparable to the USGS gauge reading.

Two panels:
  Panel 1 : full test period  (Oct 2023 – Oct 2024)
  Panel 2 : Helene zoom       (Sep 20 – Oct 5 2024)

NSE and KGE reported in title for DA and open loop.

Run on server:
    python3 plot_routed_gauge_check_f5.py
    python3 plot_routed_gauge_check_f5.py --ol-csv /path/to/routed_Q_openloop.csv
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_ROUTED_CSV = "/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct/routed_Q_test.csv"
DEFAULT_OL_CSV     = "/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct/routed_Q_openloop.csv"
DEFAULT_OUT_DIR    = "/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct"

HELENE_ZOOM_START = pd.Timestamp("2024-09-20")
HELENE_ZOOM_END   = pd.Timestamp("2024-10-05")
HELENE_PEAK_START = pd.Timestamp("2024-09-24")
HELENE_PEAK_END   = pd.Timestamp("2024-09-30")

DA_COLOR  = "tab:purple"
OL_COLOR  = "tab:gray"
OBS_COLOR = "black"


def nse(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2:
        return float("nan")
    denom = ((o - o.mean()) ** 2).sum()
    return float(1.0 - ((o - s) ** 2).sum() / denom) if denom > 0 else float("nan")


def kge(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2:
        return float("nan")
    r = float(np.corrcoef(o, s)[0, 1])
    alpha = float(s.std() / o.std()) if o.std() > 0 else float("nan")
    beta  = float(s.mean() / o.mean()) if o.mean() > 0 else float("nan")
    return 1.0 - float(np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def load_routed_csv(path):
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["date"])
    return df.set_index("date").sort_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--routed-csv", default=DEFAULT_ROUTED_CSV)
    parser.add_argument("--ol-csv",     default=DEFAULT_OL_CSV)
    parser.add_argument("--out-dir",    default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── DA routed timeseries ──────────────────────────────────────────
    df      = load_routed_csv(args.routed_csv)
    sim_col = next(c for c in df.columns if "routed" in c.lower() or "sim" in c.lower())
    obs_col = next(c for c in df.columns if "usgs" in c.lower() or "obs" in c.lower())
    sim   = df[sim_col].values.astype(float)
    obs   = df[obs_col].values.astype(float)
    dates = df.index

    # ── Open loop routed timeseries (optional) ────────────────────────
    ol_sim, ol_dates = None, None
    if args.ol_csv and os.path.exists(args.ol_csv):
        df_ol    = load_routed_csv(args.ol_csv)
        ol_col   = next(c for c in df_ol.columns
                        if "routed" in c.lower() or "sim" in c.lower() or "ol" in c.lower())
        ol_sim   = df_ol[ol_col].values.astype(float)
        ol_dates = df_ol.index
        print(f"Open loop loaded: {args.ol_csv}")
    else:
        print("No open loop CSV found — plotting DA only.")

    # ── Metrics ───────────────────────────────────────────────────────
    nse_full   = nse(obs, sim)
    kge_full   = kge(obs, sim)
    helene_mask = (dates >= HELENE_ZOOM_START) & (dates <= HELENE_ZOOM_END)
    nse_helene  = nse(obs[helene_mask], sim[helene_mask])

    ol_nse_full = ol_nse_helene = ol_kge_full = float("nan")
    if ol_sim is not None:
        ol_obs_aligned = np.interp(
            ol_dates.astype(np.int64), dates.astype(np.int64), obs,
            left=np.nan, right=np.nan)
        ol_nse_full    = nse(ol_obs_aligned, ol_sim)
        ol_kge_full    = kge(ol_obs_aligned, ol_sim)
        ol_h_mask      = (ol_dates >= HELENE_ZOOM_START) & (ol_dates <= HELENE_ZOOM_END)
        ol_nse_helene  = nse(ol_obs_aligned[ol_h_mask], ol_sim[ol_h_mask])

    # ── Title ─────────────────────────────────────────────────────────
    title_lines = [
        f"Gauge 03463300  |  F5 routed DA output vs USGS (re-kriged σ², 20% holdout)",
        f"DA:  NSE = {nse_full:.3f}  KGE = {kge_full:.3f} (full year)   NSE = {nse_helene:.3f} (Helene zoom)",
    ]
    if ol_sim is not None:
        title_lines.append(
            f"OL:  NSE = {ol_nse_full:.3f}  KGE = {ol_kge_full:.3f} (full year)   "
            f"NSE = {ol_nse_helene:.3f} (Helene zoom)")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 9))
    fig.suptitle("\n".join(title_lines), fontsize=11)

    # ── Panel 1: full period ──────────────────────────────────────────
    if ol_sim is not None:
        ax1.plot(ol_dates, ol_sim, color=OL_COLOR, lw=1.0, linestyle='--',
                 alpha=0.85, label="Open loop (no DA)  [m³/s]", zorder=2)
    ax1.plot(dates, sim, color=DA_COLOR, lw=1.2, label="F5 DA routed  [m³/s]", zorder=3)
    ax1.plot(dates, obs, color=OBS_COLOR, lw=0.8, alpha=0.7,
             label="USGS obs  [m³/s]", zorder=4)
    ax1.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax1.set_title("Full test period: Oct 2023 – Oct 2024", fontsize=10)
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(True, alpha=0.2)

    # ── Panel 2: Helene zoom ──────────────────────────────────────────
    d_h   = dates[helene_mask]
    sim_h = sim[helene_mask]
    obs_h = obs[helene_mask]

    ax2.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
                color="#ffcccc", alpha=0.20, zorder=0)
    if ol_sim is not None:
        ol_h_mask = (ol_dates >= HELENE_ZOOM_START) & (ol_dates <= HELENE_ZOOM_END)
        ax2.plot(ol_dates[ol_h_mask], ol_sim[ol_h_mask],
                 color=OL_COLOR, lw=2.0, linestyle='--',
                 label="Open loop (no DA)  [m³/s]", zorder=3)
    ax2.plot(d_h, sim_h, color=DA_COLOR, lw=2.0,
             label="F5 DA routed  [m³/s]", zorder=4)
    ax2.plot(d_h, obs_h, color=OBS_COLOR, lw=1.5,
             label="USGS obs  [m³/s]", zorder=5)
    ax2.set_ylabel("Discharge (m³/s)", fontsize=10)
    ax2.set_title("Helene zoom: Sep 20 – Oct 5 2024  (pink = Sep 24–30)", fontsize=10)
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = os.path.join(args.out_dir, "cat-1016300_routed_gauge_check_f5.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")
    print(f"DA  NSE={nse_full:.3f}  KGE={kge_full:.3f} (full year)  NSE={nse_helene:.3f} (Helene)")
    if ol_sim is not None:
        print(f"OL  NSE={ol_nse_full:.3f}  KGE={ol_kge_full:.3f} (full year)  "
              f"NSE={ol_nse_helene:.3f} (Helene)")


if __name__ == "__main__":
    main()
