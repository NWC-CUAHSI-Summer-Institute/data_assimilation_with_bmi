"""
plot_krig_obs_comparison.py

Compare Qkrig pseudo-observations used in two experiments at cat-1016300:
  - 20% holdout: catchment_ts_03463300_with_variance
      -> gauge 03463300 is HELD OUT; kriging built from ~80% of gauges
  - 1-gauge holdout: catchment_ts_no_03463300_gapfilled
      -> gauge 03463300 is HELD OUT; kriging built from all other gauges (~7,299)

In both cases gauge 03463300 is withheld from the kriging network.
The difference: 20% holdout removes ~20% of all gauges (sparser network),
while 1-gauge holdout removes only one gauge (denser network).

Overlays USGS gauge 03463300 observations for reference.

Output: krig_obs_comparison_full.png  (full test period)
        krig_obs_comparison_helene.png (Helene window zoom)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ── Paths ──────────────────────────────────────────────────────────────────
OBS_20PCT   = "/mnt/disk2/1400_sites_helene/catchment_ts_03463300_with_variance/cat-1016300.csv"
OBS_1GAUGE  = "/home/svyas/catchment_ts_no_03463300_gapfilled/cat-1016300.csv"
USGS_CSV    = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
OUT_DIR     = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

HELENE_START = pd.Timestamp("2024-09-24")
HELENE_END   = pd.Timestamp("2024-09-30")
PLOT_START   = pd.Timestamp("2019-01-01")
PLOT_END     = pd.Timestamp("2024-10-31")


def load_krig(path, label):
    df = pd.read_csv(path)
    time_col = next(c for c in df.columns if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    df[time_col] = pd.to_datetime(df[time_col])
    df = df.set_index(time_col).sort_index()
    q_col = next(c for c in df.columns if 'qkrig' in c.lower() and 'var' not in c.lower())
    print(f"  [{label}] loaded {len(df)} rows, q_col='{q_col}'")
    return df[q_col].astype(float)


def load_usgs(path):
    df = pd.read_csv(path)
    date_col = next(c for c in df.columns if c.lower() in ('datetime','date','time','timestamp'))
    q_col    = next(c for c in df.columns if 'q' in c.lower() or 'flow' in c.lower())
    df[date_col] = pd.to_datetime(df[date_col])
    s = df.set_index(date_col)[q_col].astype(float).sort_index()
    if 'mm' in q_col.lower():
        s = s * MM_H_TO_M3_S
    print(f"  [USGS] loaded {len(s)} rows, q_col='{q_col}'")
    return s


def make_plot(q20, q1g, usgs, t_start, t_end, suffix, title):
    mask20 = (q20.index >= t_start)  & (q20.index <= t_end)
    mask1g = (q1g.index >= t_start)  & (q1g.index <= t_end)
    usgs_m = (usgs.index >= t_start) & (usgs.index <= t_end)

    # Compute difference on common index
    common = q20.index.intersection(q1g.index)
    common = common[(common >= t_start) & (common <= t_end)]
    diff = q20.reindex(common) - q1g.reindex(common)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(14, 8),
                                   gridspec_kw={'height_ratios': [3, 1]},
                                   sharex=True)

    # ── Top panel: both Qkrig series + USGS ──────────────────────────────
    # Draw 1-gauge first, then 20% on top with dashed so both visible
    ax.plot(q1g.index[mask1g], q1g.values[mask1g],
            color='tomato', lw=1.2, alpha=0.9, zorder=2,
            label='Qkrig — 1-gauge holdout (03463300 removed; ~7,299 gauges)')
    ax.plot(q20.index[mask20], q20.values[mask20],
            color='steelblue', lw=1.4, alpha=0.9, linestyle='--', zorder=3,
            label='Qkrig — 20% holdout (03463300 removed; ~80% of gauges)')
    ax.plot(usgs.index[usgs_m], usgs.values[usgs_m] / MM_H_TO_M3_S,
            color='black', lw=1.6, zorder=4,
            label='USGS obs 03463300 (converted to mm/h)')

    if t_end > HELENE_START:
        hs = max(t_start, HELENE_START)
        he = min(t_end,   HELENE_END)
        ax.axvspan(hs, he, color='salmon', alpha=0.12, zorder=0, label='Helene window')

    ax.set_ylabel('q (mm/h)', fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.2)

    # ── Bottom panel: difference (20pct − 1gauge) ─────────────────────────
    ax2.plot(common, diff.values, color='purple', lw=0.8, alpha=0.85)
    ax2.axhline(0, color='black', lw=0.7, linestyle='--')
    if t_end > HELENE_START:
        hs = max(t_start, HELENE_START)
        he = min(t_end,   HELENE_END)
        ax2.axvspan(hs, he, color='salmon', alpha=0.12, zorder=0)
    ax2.set_ylabel('Difference\n(20% − 1-gauge)\nmm/h', fontsize=9)
    ax2.set_xlabel('Date', fontsize=11)
    ax2.grid(True, alpha=0.2)

    if suffix == 'helene':
        ax2.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    else:
        ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=30, ha='right')

    plt.tight_layout()
    out = os.path.join(OUT_DIR, f'krig_obs_comparison_{suffix}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading Qkrig obs...')
    q20 = load_krig(OBS_20PCT,  '20pct holdout')
    q1g = load_krig(OBS_1GAUGE, '1-gauge holdout')
    print('Loading USGS obs...')
    usgs = load_usgs(USGS_CSV)

    print('Plotting full test period...')
    make_plot(q20, q1g, usgs, PLOT_START, PLOT_END, 'full',
              'Qkrig pseudo-observations at cat-1016300\n'
              '20% holdout vs 1-gauge holdout — full test period (Oct 2023–Oct 2024)')

    print('Plotting Helene window...')
    make_plot(q20, q1g, usgs, HELENE_START, HELENE_END, 'helene',
              'Qkrig pseudo-observations at cat-1016300\n'
              '20% holdout vs 1-gauge holdout — Hurricane Helene (Sep 24–30, 2024)')

    # Print summary stats for the Helene window
    h20  = q20[(q20.index  >= HELENE_START) & (q20.index  <= HELENE_END)]
    h1g  = q1g[(q1g.index  >= HELENE_START) & (q1g.index  <= HELENE_END)]
    husgs= usgs[(usgs.index >= HELENE_START) & (usgs.index <= HELENE_END)]
    print('\n── Helene window peak (mm/h) ────────────────')
    print(f'  Qkrig 20% holdout : {h20.max():.3f} mm/h')
    print(f'  Qkrig 1-gauge     : {h1g.max():.3f} mm/h')
    print(f'  USGS (converted)  : {(husgs.max()/MM_H_TO_M3_S):.3f} mm/h')


if __name__ == '__main__':
    main()
