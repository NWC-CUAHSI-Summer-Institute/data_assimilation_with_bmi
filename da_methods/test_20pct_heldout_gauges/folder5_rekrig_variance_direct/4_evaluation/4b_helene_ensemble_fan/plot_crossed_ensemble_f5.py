"""
plot_crossed_ensemble_f5.py

Plots the 600-member crossed ensemble (30 forcing × 20 hydro-state members)
routed to gauge 03463300 during Hurricane Helene (Sep 24-30 2024).

Reads the already-routed parquet produced by run_route_crossed_ensemble.py
(Q in m³/s at gauge, 600 member columns). Aggregates across all Helene
issue times to show total forecast uncertainty — both forcing and initial
state uncertainty active simultaneously.

Inputs:
  <crossed-dir>/routed_crossed_ensemble.parquet   (m³/s, 600 members)
  USGS obs CSV

Output:
  <out-dir>/cat-1016300_2ab_crossed_600member_routed.png

Run on server (troute env):
    python3 plot_crossed_ensemble_f5.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

DEFAULT_CROSSED_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_crossed_ensemble_routed"
DEFAULT_USGS_CSV    = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
DEFAULT_OUT_DIR     = None   # falls back to --crossed-dir

WATERSHED_AREA_KM2 = 113.18
HELENE_START       = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END         = pd.Timestamp("2024-09-30 23:00:00")
HELENE_PEAK_START  = pd.Timestamp("2024-09-26 12:00:00")
HELENE_PEAK_END    = pd.Timestamp("2024-09-28 00:00:00")

CROSSED_COLOR = "#7b4fa6"   # purple — combined uncertainty


def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns
                    if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    q_col    = next(c for c in df.columns
                    if any(k in c.lower() for k in ('q', 'flow', 'discharge')))
    df[date_col] = pd.to_datetime(df[date_col])
    obs = df.set_index(date_col)[q_col].astype(float)
    if 'mm' in q_col.lower():
        obs = obs * WATERSHED_AREA_KM2 * 1000.0 / 3600.0
    return obs


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--crossed-dir', default=DEFAULT_CROSSED_DIR)
    parser.add_argument('--usgs-csv',    default=DEFAULT_USGS_CSV)
    parser.add_argument('--out-dir',     default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    out_dir = args.out_dir or args.crossed_dir
    os.makedirs(out_dir, exist_ok=True)

    # ── Load routed crossed ensemble ──────────────────────────────────────────
    pq_path = os.path.join(args.crossed_dir, "routed_crossed_ensemble.parquet")
    if not os.path.exists(pq_path):
        raise FileNotFoundError(
            f"Routed crossed ensemble not found: {pq_path}\n"
            f"Run run_route_crossed_ensemble.py first.")
    print(f"Loading routed crossed ensemble from {pq_path} ...")
    df = pd.read_parquet(pq_path)
    df['issue_time'] = pd.to_datetime(df['issue_time'])
    print(f"  Shape: {df.shape}  |  issue times: {df['issue_time'].nunique()}")

    member_cols = sorted(c for c in df.columns if c.startswith('member_'))
    n_members   = len(member_cols)
    print(f"  {n_members} member columns (expect 600)")

    # ── Filter to Helene window issue times ───────────────────────────────────
    helene_mask = (df['issue_time'] >= HELENE_START) & (df['issue_time'] <= HELENE_END)
    df_h = df[helene_mask].copy()
    print(f"  Helene issue times: {df_h['issue_time'].nunique()}")

    # ── Compute valid_time and aggregate across all issue times + members ─────
    df_h['valid_time'] = (df_h['issue_time'] +
                          pd.to_timedelta(df_h['lead_hour'], unit='h'))

    # For each valid_time: pool all 600 member values from all issue times
    # that have a forecast window covering that valid_time
    from collections import defaultdict
    q_by_vt = defaultdict(list)
    for _, row in df_h.iterrows():
        vt   = row['valid_time']
        vals = row[member_cols].values.astype(float)
        vals = np.maximum(vals, 0.0)
        q_by_vt[vt].extend(vals.tolist())

    times  = sorted(q_by_vt.keys())
    median = np.array([np.median(q_by_vt[t])       for t in times])
    q5     = np.array([np.percentile(q_by_vt[t],  5) for t in times])
    q95    = np.array([np.percentile(q_by_vt[t], 95) for t in times])
    q25    = np.array([np.percentile(q_by_vt[t], 25) for t in times])
    q75    = np.array([np.percentile(q_by_vt[t], 75) for t in times])

    print(f"  Valid times aggregated: {len(times)}")
    print(f"  Peak grand median: {median.max():.1f} m³/s")

    # ── USGS obs ──────────────────────────────────────────────────────────────
    print("Loading USGS obs...")
    usgs     = load_usgs(args.usgs_csv)
    obs_mask = ((usgs.index >= HELENE_START) &
                (usgs.index <= HELENE_END + pd.Timedelta(hours=24)))

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(17, 6))

    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color='salmon', alpha=0.12, zorder=0, label='_nolegend_')
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color='firebrick', ha='left', va='top')

    # Outer band: 5th–95th
    ax.fill_between(times, q5, q95,
                    color=CROSSED_COLOR, alpha=0.18, linewidth=0, zorder=2,
                    label=f"600-member spread [5th–95th]")
    # Inner band: 25th–75th
    ax.fill_between(times, q25, q75,
                    color=CROSSED_COLOR, alpha=0.30, linewidth=0, zorder=3,
                    label=f"600-member spread [25th–75th]")
    # Grand median
    ax.plot(times, median,
            color=CROSSED_COLOR, lw=2.4, zorder=4,
            label="600-member grand median")

    # USGS obs
    ax.plot(usgs[obs_mask].index, usgs[obs_mask].values,
            color='black', lw=2.6, zorder=6, label="USGS obs")

    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_title(
        f"600-member crossed ensemble (30 forcing × 20 hydro-state)  |  DA on  |  "
        f"All 21 catchments T-Route routed to gauge 03463300\n"
        f"F5 re-kriged σ²  |  20% gauge holdout  |  Sep 24–30 2024  |  "
        f"Aggregated across all Helene issue times",
        fontsize=10)
    ax.set_xlim(HELENE_START, HELENE_END + pd.Timedelta(hours=24))
    ax.legend(fontsize=10, loc='upper left', framealpha=0.92)
    ax.grid(True, alpha=0.22)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha='right', fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "cat-1016300_2ab_crossed_600member_routed.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
