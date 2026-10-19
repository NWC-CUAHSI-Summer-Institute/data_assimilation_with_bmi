"""
Gauge-level lead-time forecast decay curve — F5 (re-kriged variance).

Routes 21-catchment forecast CSVs through T-route to USGS gauge 03463300 and
plots RMSE vs lead time (DA solid, open-loop dashed) against USGS observed Q.
Same 2-panel layout as the catchment-level plot, plus a regime split for the
Helene window.

The routing step is done separately by route_lead_time_forecasts.py;
this script reads the resulting parquets and produces the plot.

Inputs:
    <route>/routed_leadtime_da_full.parquet
    <route>/routed_leadtime_openloop_full.parquet
    Wide format: issue_time, lead_hour, member_00..member_19 (q_gauge_m3s)
Obs:
    /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv

Output:
    <route>/lead_time_decay_gauge_pooled.png
    <route>/lead_time_decay_gauge_by_regime.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DA_COLOR  = "tab:purple"
OL_COLOR  = "tab:gray"

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct_leadtime_routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 23:00:00")
STORM_OBS_THRESHOLD_M3S   = 20.0
LOWFLOW_OBS_THRESHOLD_M3S = 5.0

USGS_HELENE_PEAK_M3S = 1886.0


def load_parquet_long(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_parquet(path)

    if 'q_gauge_m3s' in df.columns and 'member' in df.columns:
        df['issue_time'] = pd.to_datetime(df['issue_time'])
        return df[['issue_time', 'lead_hour', 'member', 'q_gauge_m3s']]

    q_col_candidates = [c for c in df.columns if c.lower() in
                        ('q_gauge_m3s', 'q_m3s', 'q_gauge', 'q')]
    if 'member' in df.columns and q_col_candidates:
        qc = q_col_candidates[0]
        df['issue_time'] = pd.to_datetime(df['issue_time'])
        out = df[['issue_time', 'lead_hour', 'member', qc]].copy()
        return out.rename(columns={qc: 'q_gauge_m3s'})

    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    if not member_cols:
        raise ValueError(f"Couldn't identify member columns in {path}. "
                         f"Columns: {df.columns.tolist()}")
    df['issue_time'] = pd.to_datetime(df['issue_time'])
    keep = ['issue_time', 'lead_hour']
    if 'valid_time' in df.columns:
        keep.append('valid_time')
    return df[keep + member_cols].melt(
        id_vars=keep, value_vars=member_cols,
        var_name='member', value_name='q_gauge_m3s')


def load_usgs_obs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_candidates = [c for c in df.columns if c.lower() in
                       ('datetime', 'date', 'time', 'timestamp')]
    q_candidates    = [c for c in df.columns if 'q' in c.lower()
                       or 'flow' in c.lower() or 'discharge' in c.lower()]
    if not date_candidates or not q_candidates:
        raise ValueError(f"Could not identify date/Q columns in {usgs_csv}. "
                         f"Columns: {df.columns.tolist()}")
    dc, qc = date_candidates[0], q_candidates[0]
    df[dc] = pd.to_datetime(df[dc])
    series = df.set_index(dc)[qc].astype(float)
    if 'mm' in qc.lower():
        series = series * MM_H_TO_M3_S
        print(f"  Converted obs from mm/h (column '{qc}') to m3/s "
              f"using area = {WATERSHED_AREA_KM2} km2 (x {MM_H_TO_M3_S:.4f})")
    else:
        print(f"  Loaded obs as-is from column '{qc}' (assumed m3/s)")
    return series


def metrics_by_lead(df_long, obs_series, issue_mask=None):
    df = df_long.copy()
    if issue_mask is not None:
        df = df[df['issue_time'].isin(issue_mask)]
    if len(df) == 0:
        return None
    df['valid_time'] = df['issue_time'] + pd.to_timedelta(df['lead_hour'], unit='h')
    df['obs'] = df['valid_time'].map(obs_series)
    df = df.dropna(subset=['obs'])
    if len(df) == 0:
        return None

    grouped  = df.groupby(['issue_time', 'lead_hour'])
    ens_mean = grouped['q_gauge_m3s'].mean().reset_index(name='ens_mean')
    ens_std  = grouped['q_gauge_m3s'].std().reset_index(name='ens_std')
    obs_per  = grouped['obs'].first().reset_index(name='obs')
    panel    = ens_mean.merge(ens_std, on=['issue_time', 'lead_hour']).merge(
        obs_per, on=['issue_time', 'lead_hour'])

    leads = sorted(panel['lead_hour'].unique())
    rmse_mean = np.full(len(leads), np.nan)
    rmse_lo   = np.full(len(leads), np.nan)
    rmse_hi   = np.full(len(leads), np.nan)
    std_mean  = np.full(len(leads), np.nan)
    for i, L in enumerate(leads):
        psub = panel[panel['lead_hour'] == L]
        if len(psub) == 0:
            continue
        err = psub['ens_mean'].values - psub['obs'].values
        rmse_mean[i] = float(np.sqrt(np.mean(err ** 2)))
        dsub = df[df['lead_hour'] == L]
        member_rmses = []
        for m, g in dsub.groupby('member'):
            err_m = g['q_gauge_m3s'].values - g['obs'].values
            if len(err_m) >= 1:
                member_rmses.append(float(np.sqrt(np.mean(err_m ** 2))))
        if member_rmses:
            rmse_lo[i] = float(np.nanmin(member_rmses))
            rmse_hi[i] = float(np.nanmax(member_rmses))
        std_mean[i] = float(np.nanmean(psub['ens_std'].values))
    return np.asarray(leads), rmse_mean, rmse_lo, rmse_hi, std_mean


def plot_two_panel(out_path, da_metrics, ol_metrics, title):
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                          gridspec_kw={"height_ratios": [1.4, 1.0]})
    leads, rmse_da, rmse_da_lo, rmse_da_hi, std_da = da_metrics
    _,     rmse_ol, rmse_ol_lo, rmse_ol_hi, std_ol = ol_metrics

    ax_top.fill_between(leads, rmse_da_lo, rmse_da_hi,
                        color=DA_COLOR, alpha=0.18, zorder=2,
                        label="DA — member spread (min/max RMSE)")
    ax_top.fill_between(leads, rmse_ol_lo, rmse_ol_hi,
                        color=OL_COLOR, alpha=0.18, zorder=2,
                        label="Open-loop — member spread")
    ax_top.plot(leads, rmse_da, color=DA_COLOR, lw=2.4, marker='o',
                zorder=4, label="DA on (ensemble-mean RMSE)")
    ax_top.plot(leads, rmse_ol, color=OL_COLOR, lw=2.4, marker='s',
                linestyle='--', zorder=4, label="Open-loop (ensemble-mean RMSE)")
    ax_top.set_ylabel("RMSE vs USGS obs (m3/s)", fontsize=11)
    ax_top.set_title(title, fontsize=12)
    ax_top.grid(True, alpha=0.25)
    ax_top.legend(fontsize=9, loc='upper left', frameon=True, framealpha=0.92)

    ax_bot.plot(leads, std_da, color=DA_COLOR, lw=2.2, marker='o',
                label="DA — mean ensemble std")
    ax_bot.plot(leads, std_ol, color=OL_COLOR, lw=2.2, marker='s',
                linestyle='--', label="Open-loop — mean ensemble std")
    ax_bot.set_xlabel("Forecast lead time (hours)", fontsize=11)
    ax_bot.set_ylabel("Mean ensemble std (m3/s)", fontsize=11)
    ax_bot.set_xticks(np.arange(1, 19))
    ax_bot.grid(True, alpha=0.25)
    ax_bot.legend(fontsize=9, loc='upper left', frameon=True, framealpha=0.92)

    plt.tight_layout()
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def plot_three_panel_regime(out_path, regime_metrics):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10),
                              gridspec_kw={"height_ratios": [1.4, 1.0]})
    for col, (label, n_issue, da_m, ol_m) in enumerate(regime_metrics):
        ax_rmse = axes[0, col]
        ax_std  = axes[1, col]
        if da_m is None or ol_m is None:
            ax_rmse.set_title(f"{label}\n(no issue times match)", fontsize=11)
            ax_rmse.axis('off'); ax_std.axis('off')
            continue
        leads, rmse_da, rmse_da_lo, rmse_da_hi, std_da = da_m
        _,     rmse_ol, rmse_ol_lo, rmse_ol_hi, std_ol = ol_m
        ax_rmse.fill_between(leads, rmse_da_lo, rmse_da_hi,
                             color=DA_COLOR, alpha=0.18, zorder=2)
        ax_rmse.fill_between(leads, rmse_ol_lo, rmse_ol_hi,
                             color=OL_COLOR, alpha=0.18, zorder=2)
        ax_rmse.plot(leads, rmse_da, color=DA_COLOR, lw=2.2, marker='o',
                     zorder=4, label="DA on")
        ax_rmse.plot(leads, rmse_ol, color=OL_COLOR, lw=2.2, marker='s',
                     linestyle='--', zorder=4, label="Open-loop")
        ax_rmse.set_title(f"{label}\n({n_issue} issue times)", fontsize=11)
        ax_rmse.set_ylabel("RMSE (m3/s)", fontsize=10)
        ax_rmse.grid(True, alpha=0.25)
        ax_rmse.legend(fontsize=8, loc='upper left', frameon=True, framealpha=0.92)

        ax_std.plot(leads, std_da, color=DA_COLOR, lw=2.0, marker='o', label="DA")
        ax_std.plot(leads, std_ol, color=OL_COLOR, lw=2.0, marker='s',
                    linestyle='--', label="Open-loop")
        ax_std.set_xlabel("Forecast lead time (hours)", fontsize=10)
        ax_std.set_ylabel("Mean ensemble std (m3/s)", fontsize=10)
        ax_std.set_xticks(np.arange(1, 19))
        ax_std.grid(True, alpha=0.25)
        ax_std.legend(fontsize=8, loc='upper left', frameon=True, framealpha=0.92)

    fig.suptitle("Gauge-level lead-time error decay by regime — USGS 03463300\n"
                 "F5 (re-kriged variance direct)",
                 fontsize=13, y=0.995)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--route-dir', default=DEFAULT_ROUTE_DIR)
    parser.add_argument('--out-dir',   default=None)
    parser.add_argument('--usgs-csv',  default=DEFAULT_USGS_CSV)
    parser.add_argument('--da-name',   default='routed_leadtime_da_full.parquet')
    parser.add_argument('--ol-name',   default='routed_leadtime_openloop_full.parquet')
    args = parser.parse_args()
    out_dir = args.out_dir or args.route_dir
    os.makedirs(out_dir, exist_ok=True)

    da_path = os.path.join(args.route_dir, args.da_name)
    ol_path = os.path.join(args.route_dir, args.ol_name)
    df_da = load_parquet_long(da_path)
    df_ol = load_parquet_long(ol_path)
    obs_series = load_usgs_obs(args.usgs_csv)

    print(f"DA rows: {len(df_da):,} | unique issue_times: {df_da['issue_time'].nunique()}")
    print(f"OL rows: {len(df_ol):,} | unique issue_times: {df_ol['issue_time'].nunique()}")
    print(f"USGS Helene peak (reference): {USGS_HELENE_PEAK_M3S:.0f} m3/s")

    da_pooled = metrics_by_lead(df_da, obs_series)
    ol_pooled = metrics_by_lead(df_ol, obs_series)
    out_pooled = os.path.join(out_dir, "lead_time_decay_gauge_pooled.png")
    plot_two_panel(out_pooled, da_pooled, ol_pooled,
                   title="Gauge-level lead-time decay — USGS 03463300\n"
                         "F5 (re-kriged variance direct) — all issue times pooled")

    issue_times  = pd.to_datetime(sorted(df_da['issue_time'].unique()))
    obs_at_issue = pd.Series(issue_times, index=issue_times).map(obs_series)

    helene_mask  = (issue_times >= HELENE_START) & (issue_times <= HELENE_END)
    storm_mask   = obs_at_issue > STORM_OBS_THRESHOLD_M3S
    lowflow_mask = obs_at_issue < LOWFLOW_OBS_THRESHOLD_M3S

    regime_metrics = []
    for label, mask in [
        ("Helene window (Sep 24-28 2024)",                                         helene_mask),
        (f"Storm hours (obs(t0) > {STORM_OBS_THRESHOLD_M3S:.0f} m3/s)",           storm_mask),
        (f"Low flow (obs(t0) < {LOWFLOW_OBS_THRESHOLD_M3S:.0f} m3/s)",            lowflow_mask),
    ]:
        kept = set(pd.to_datetime(issue_times[mask]))
        da_m = metrics_by_lead(df_da, obs_series, kept)
        ol_m = metrics_by_lead(df_ol, obs_series, kept)
        regime_metrics.append((label, len(kept), da_m, ol_m))

    out_regime = os.path.join(out_dir, "lead_time_decay_gauge_by_regime.png")
    plot_three_panel_regime(out_regime, regime_metrics)


if __name__ == "__main__":
    main()
