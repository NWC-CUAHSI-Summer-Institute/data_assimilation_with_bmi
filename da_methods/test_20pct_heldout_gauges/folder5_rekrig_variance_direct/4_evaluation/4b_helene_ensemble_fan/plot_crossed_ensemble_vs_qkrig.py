"""
plot_crossed_ensemble_vs_qkrig.py

Plots the 600-member crossed ensemble (30 forcing × 20 hydro-state) for
cat-1016300 against Qkrig — the spatially interpolated observation that
EnKF actually assimilated — rather than the USGS gauge.

This is the catchment-level view: shows how well the ensemble brackets
the observation the DA was targeting, before any routing.

Reads:
  <crossed-dir>/cat-1016300/cat-1016300_crossed_ensemble.parquet
      Columns: issue_time, lead_hour, member_0000..member_0599  (mm/h)

  <da-dir>/cat-1016300/cat-1016300_test_results.csv
      Columns: date, sim_mm_h, obs_mm_h (Qkrig), precip_mm_h

Output:
  <out-dir>/cat-1016300_crossed_600member_vs_qkrig_helene.png

Run on server (troute env):
    python3 plot_crossed_ensemble_vs_qkrig.py
"""

import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

DEFAULT_CROSSED_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble"
DEFAULT_DA_DIR      = "/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct"
DEFAULT_OUT_DIR     = None   # falls back to --crossed-dir

CAT_ID = "cat-1016300"

HELENE_START      = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END        = pd.Timestamp("2024-09-30 23:00:00")
HELENE_PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_PEAK_END   = pd.Timestamp("2024-09-28 00:00:00")

ENSEMBLE_COLOR = "#4878CF"   # blue
QKRIG_COLOR    = "#d62728"   # red — Qkrig obs


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--crossed-dir', default=DEFAULT_CROSSED_DIR)
    parser.add_argument('--da-dir',      default=DEFAULT_DA_DIR)
    parser.add_argument('--out-dir',     default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    out_dir = args.out_dir or os.path.join(args.crossed_dir, CAT_ID)
    os.makedirs(out_dir, exist_ok=True)

    # ── Load crossed ensemble parquet ─────────────────────────────────────────
    pq_path = os.path.join(args.crossed_dir, CAT_ID, f"{CAT_ID}_crossed_ensemble.parquet")
    if not os.path.exists(pq_path):
        raise FileNotFoundError(f"Not found: {pq_path}")
    print(f"Loading crossed ensemble: {pq_path}")
    df = pd.read_parquet(pq_path)
    df['issue_time'] = pd.to_datetime(df['issue_time'])

    member_cols = sorted(c for c in df.columns if c.startswith('member_'))
    n_members   = len(member_cols)
    print(f"  {n_members} members | {df['issue_time'].nunique()} issue times")

    # ── Filter to Helene issue times ──────────────────────────────────────────
    helene_mask = (df['issue_time'] >= HELENE_START) & (df['issue_time'] <= HELENE_END)
    df_h = df[helene_mask].copy()
    df_h['valid_time'] = (df_h['issue_time'] +
                          pd.to_timedelta(df_h['lead_hour'], unit='h'))
    print(f"  Helene issue times: {df_h['issue_time'].nunique()}")

    # ── Aggregate: pool all member values per valid_time ──────────────────────
    q_by_vt = defaultdict(list)
    for _, row in df_h.iterrows():
        vt   = row['valid_time']
        vals = np.maximum(row[member_cols].values.astype(float), 0.0)
        q_by_vt[vt].extend(vals.tolist())

    times  = sorted(q_by_vt.keys())
    median = np.array([np.median(q_by_vt[t])          for t in times])
    q5     = np.array([np.percentile(q_by_vt[t],  5)  for t in times])
    q95    = np.array([np.percentile(q_by_vt[t], 95)  for t in times])
    q25    = np.array([np.percentile(q_by_vt[t], 25)  for t in times])
    q75    = np.array([np.percentile(q_by_vt[t], 75)  for t in times])
    print(f"  Peak ensemble median: {median.max():.3f} mm/h")

    # ── Load Qkrig (obs_mm_h from DA test results) ────────────────────────────
    qkrig_path = os.path.join(args.da_dir, CAT_ID, f"{CAT_ID}_test_results.csv")
    if not os.path.exists(qkrig_path):
        raise FileNotFoundError(f"Qkrig CSV not found: {qkrig_path}")
    print(f"Loading Qkrig: {qkrig_path}")
    df_da  = pd.read_csv(qkrig_path)
    df_da['date'] = pd.to_datetime(df_da['date'])
    df_da  = df_da.set_index('date').sort_index()
    qkrig  = df_da['obs_mm_h'].astype(float)

    obs_mask = ((qkrig.index >= HELENE_START) &
                (qkrig.index <= HELENE_END + pd.Timedelta(hours=24)))
    print(f"  Qkrig peak (Helene): {qkrig[obs_mask].max():.3f} mm/h")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(17, 6))

    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color='salmon', alpha=0.12, zorder=0, label='_nolegend_')
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color='firebrick', ha='left', va='top')

    # Outer band: 5th–95th
    ax.fill_between(times, q5, q95,
                    color=ENSEMBLE_COLOR, alpha=0.18, linewidth=0, zorder=2,
                    label=f"600-member spread [5th–95th]")
    # Inner band: 25th–75th
    ax.fill_between(times, q25, q75,
                    color=ENSEMBLE_COLOR, alpha=0.32, linewidth=0, zorder=3,
                    label=f"600-member spread [25th–75th]")
    # Grand median
    ax.plot(times, median,
            color=ENSEMBLE_COLOR, lw=2.2, zorder=4,
            label="Ensemble median (600 members)")

    # Qkrig
    ax.plot(qkrig[obs_mask].index, qkrig[obs_mask].values,
            color=QKRIG_COLOR, lw=2.2, zorder=6,
            label="Qkrig obs (DA target)")

    ax.set_ylabel("Streamflow (mm/h)", fontsize=11)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_title(
        f"600-member crossed ensemble vs Qkrig — {CAT_ID}  |  catchment level (mm/h)\n"
        f"30 forcing × 20 hydro-state members  |  F5 re-kriged σ²  |  20% gauge holdout  |  "
        f"Aggregated across all Helene issue times (Sep 24–30 2024)",
        fontsize=10)
    ax.set_xlim(HELENE_START, HELENE_END + pd.Timedelta(hours=24))
    ax.legend(fontsize=10, loc='upper left', framealpha=0.92)
    ax.grid(True, alpha=0.22)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha='right', fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{CAT_ID}_crossed_600member_vs_qkrig_helene.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
