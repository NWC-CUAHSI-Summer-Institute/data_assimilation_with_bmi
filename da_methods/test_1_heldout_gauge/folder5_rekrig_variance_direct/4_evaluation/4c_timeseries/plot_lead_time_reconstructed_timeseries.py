"""
Reconstructed time series at gauge 03463300 — F5 (re-kriged variance direct).

For each target hour t, pool ALL forecasts that land on t (across many
issue_times x many lead_hours x 20 members). Compute median + 5th/95th
percentile envelope. Plot DA and open-loop reconstructions overlaid on USGS
obs, with the Helene window shaded.

Input parquets (from route_lead_time_forecasts.py):
  <route-dir>/routed_leadtime_da_full.parquet
  <route-dir>/routed_leadtime_openloop_full.parquet
  Wide format: issue_time, lead_hour, member_00..member_19 (q in m3/s)

Obs:
  /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

Output:
  <out-dir>/lead_time_reconstructed_timeseries.png
  <out-dir>/lead_time_reconstructed_timeseries_helene.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

DA_COLOR  = "tab:purple"
OL_COLOR  = "tab:gray"
OBS_COLOR = "black"

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

PLOT_START_DEFAULT = pd.Timestamp("2024-09-10 00:00:00")
PLOT_END_DEFAULT   = pd.Timestamp("2024-10-10 23:00:00")

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-30 23:00:00")


def load_parquet_long(path):
    df = pd.read_parquet(path)
    df['issue_time'] = pd.to_datetime(df['issue_time'])
    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    if not member_cols:
        raise ValueError(f"No member columns in {path}")
    long = df.melt(
        id_vars=['issue_time', 'lead_hour'],
        value_vars=member_cols,
        var_name='member',
        value_name='q_m3s',
    )
    long['valid_time'] = (long['issue_time']
                          + pd.to_timedelta(long['lead_hour'], unit='h'))
    return long


def reconstruct(long_df):
    grouped = long_df.groupby('valid_time')['q_m3s'].agg(
        median='median',
        p05=lambda s: s.quantile(0.05),
        p95=lambda s: s.quantile(0.95),
        count='count',
    ).reset_index()
    return grouped


def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next((c for c in df.columns
                     if c.lower() in ('datetime', 'date', 'time', 'timestamp')), None)
    q_col    = next((c for c in df.columns
                     if 'q' in c.lower() or 'flow' in c.lower()
                     or 'discharge' in c.lower()), None)
    if not date_col or not q_col:
        raise ValueError(f"Couldn't id date/Q columns in {usgs_csv}: {df.columns.tolist()}")
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if 'mm' in q_col.lower():
        series = series * MM_H_TO_M3_S
        print(f"  Converted obs from mm/h (column '{q_col}') to m3/s "
              f"using area = {WATERSHED_AREA_KM2} km2 (x {MM_H_TO_M3_S:.4f})")
    else:
        print(f"  Loaded obs as-is from column '{q_col}' (assumed m3/s)")
    return series


def nse(obs, sim):
    m = np.isfinite(obs) & np.isfinite(sim)
    if m.sum() < 2:
        return np.nan
    o, s = obs[m], sim[m]
    denom = ((o - o.mean()) ** 2).sum()
    if denom < 1e-10:
        return np.nan
    return 1.0 - float(((o - s) ** 2).sum() / denom)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--route-dir',  default=DEFAULT_ROUTE_DIR)
    parser.add_argument('--out-dir',    default=None)
    parser.add_argument('--usgs-csv',   default=DEFAULT_USGS_CSV)
    parser.add_argument('--plot-start', default=str(PLOT_START_DEFAULT))
    parser.add_argument('--plot-end',   default=str(PLOT_END_DEFAULT))
    parser.add_argument('--da-name',    default='routed_leadtime_da_full.parquet')
    parser.add_argument('--ol-name',    default='routed_leadtime_openloop_full.parquet')
    args = parser.parse_args()

    out_dir = args.out_dir or args.route_dir
    os.makedirs(out_dir, exist_ok=True)
    plot_start = pd.Timestamp(args.plot_start)
    plot_end   = pd.Timestamp(args.plot_end)

    print(f"Loading DA   parquet: {args.da_name}")
    da_long = load_parquet_long(os.path.join(args.route_dir, args.da_name))
    print(f"Loading OL   parquet: {args.ol_name}")
    ol_long = load_parquet_long(os.path.join(args.route_dir, args.ol_name))
    print(f"Loading USGS obs:     {args.usgs_csv}")
    obs = load_usgs(args.usgs_csv)

    print("Reconstructing DA time series...")
    da_rec = reconstruct(da_long)
    print(f"  {len(da_rec)} unique target hours")
    print("Reconstructing OL time series...")
    ol_rec = reconstruct(ol_long)

    da_rec = da_rec[(da_rec['valid_time'] >= plot_start) & (da_rec['valid_time'] <= plot_end)].copy()
    ol_rec = ol_rec[(ol_rec['valid_time'] >= plot_start) & (ol_rec['valid_time'] <= plot_end)].copy()
    da_rec['obs'] = da_rec['valid_time'].map(obs)
    ol_rec['obs'] = ol_rec['valid_time'].map(obs)

    nse_da = nse(da_rec['obs'].values, da_rec['median'].values)
    nse_ol = nse(ol_rec['obs'].values, ol_rec['median'].values)
    print(f"NSE (plot window) — DA: {nse_da:.3f}  |  OL: {nse_ol:.3f}")

    fig, ax = plt.subplots(figsize=(16, 5.5))
    fig.patch.set_facecolor("white")

    ax.axvspan(HELENE_START, HELENE_END, color="firebrick", alpha=0.08, zorder=0)
    ax.text(HELENE_START + pd.Timedelta(hours=12), 1.0, "Helene",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color="firebrick", ha="left", va="top", fontweight="bold")

    ax.fill_between(ol_rec['valid_time'], ol_rec['p05'], ol_rec['p95'],
                    color=OL_COLOR, alpha=0.18, zorder=2,
                    label="Open-loop [5th-95th pctile envelope]")
    ax.plot(ol_rec['valid_time'], ol_rec['median'],
            color=OL_COLOR, lw=1.2, linestyle='--', zorder=3,
            label=f"Open-loop median (NSE={nse_ol:.3f})")

    ax.fill_between(da_rec['valid_time'], da_rec['p05'], da_rec['p95'],
                    color=DA_COLOR, alpha=0.22, zorder=4,
                    label="DA [5th-95th pctile envelope]")
    ax.plot(da_rec['valid_time'], da_rec['median'],
            color=DA_COLOR, lw=1.4, zorder=5,
            label=f"DA median (NSE={nse_da:.3f})")

    obs_window = obs.loc[plot_start:plot_end]
    finite = np.isfinite(obs_window.values)
    ax.scatter(obs_window.index[finite], obs_window.values[finite],
               s=6, color=OBS_COLOR, linewidths=0, zorder=6, label="USGS obs")

    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("Discharge (m3/s)", fontsize=11)
    ax.set_title(
        f"Reconstructed forecast time series at USGS 03463300 — F5 (re-kriged variance)\n"
        f"Overlapping-leads pool from EnKF forecast ensemble "
        f"({plot_start.date()} - {plot_end.date()})",
        fontsize=12,
    )
    ax.set_xlim(plot_start, plot_end)
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    ax.grid(True, alpha=0.25, lw=0.4)
    ax.legend(loc='upper left', fontsize=9, frameon=True, framealpha=0.92)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "lead_time_reconstructed_timeseries.png")
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")

    helene_zoom_start = pd.Timestamp("2024-09-24 00:00:00")
    helene_zoom_end   = pd.Timestamp("2024-09-29 23:00:00")
    da_z = da_rec[(da_rec['valid_time'] >= helene_zoom_start)
                  & (da_rec['valid_time'] <= helene_zoom_end)]
    ol_z = ol_rec[(ol_rec['valid_time'] >= helene_zoom_start)
                  & (ol_rec['valid_time'] <= helene_zoom_end)]
    nse_da_z = nse(da_z['obs'].values, da_z['median'].values)
    nse_ol_z = nse(ol_z['obs'].values, ol_z['median'].values)
    print(f"NSE (Helene zoom) — DA: {nse_da_z:.3f}  |  OL: {nse_ol_z:.3f}")

    fig, ax = plt.subplots(figsize=(14, 5.5))
    fig.patch.set_facecolor("white")
    ax.fill_between(ol_z['valid_time'], ol_z['p05'], ol_z['p95'],
                    color=OL_COLOR, alpha=0.18, zorder=2, label="Open-loop [5th-95th]")
    ax.plot(ol_z['valid_time'], ol_z['median'], color=OL_COLOR,
            lw=1.4, linestyle='--', zorder=3,
            label=f"Open-loop median (NSE={nse_ol_z:.3f})")
    ax.fill_between(da_z['valid_time'], da_z['p05'], da_z['p95'],
                    color=DA_COLOR, alpha=0.22, zorder=4, label="DA [5th-95th]")
    ax.plot(da_z['valid_time'], da_z['median'], color=DA_COLOR,
            lw=1.6, zorder=5, label=f"DA median (NSE={nse_da_z:.3f})")
    obs_z = obs.loc[helene_zoom_start:helene_zoom_end]
    finite = np.isfinite(obs_z.values)
    ax.scatter(obs_z.index[finite], obs_z.values[finite],
               s=12, color=OBS_COLOR, linewidths=0, zorder=6, label="USGS obs")
    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("Discharge (m3/s)", fontsize=11)
    ax.set_title(
        "Reconstructed forecast time series at USGS 03463300 — Helene zoom\n"
        "F5 (re-kriged variance) | Sep 24-29, 2024 | overlapping-leads pool",
        fontsize=12,
    )
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    ax.grid(True, alpha=0.25, lw=0.4)
    ax.legend(loc='upper left', fontsize=9, frameon=True, framealpha=0.92)
    plt.tight_layout()
    out_zoom = os.path.join(out_dir, "lead_time_reconstructed_timeseries_helene.png")
    plt.savefig(out_zoom, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_zoom}")


if __name__ == "__main__":
    main()
