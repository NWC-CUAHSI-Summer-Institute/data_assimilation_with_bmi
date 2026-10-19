"""
plot_lead_time_nse_gauge_f5.py

Gauge-level lead-time NSE decay for F5 (re-kriged variance, 20% holdout).

Reads routed lead-time forecast parquets produced by run_lead_time_forecast_sweep.py
(routed through run_route_leadtime_forecasts.py) and plots NSE vs lead time for
the DA-on (F5) and open-loop ensembles at USGS gauge 03463300.

NSE per lead time is computed by pooling all (issue_time, valid_time) pairs for
that lead hour across the full test period, then:
    NSE = 1 - sum((obs - ens_mean)^2) / sum((obs - mean(obs))^2)

A second panel shows the regime split (Helene window vs. storm hours vs. low flow).

Inputs:
    <route-dir>/routed_leadtime_da_full.parquet
    <route-dir>/routed_leadtime_openloop_full.parquet

Output:
    <out-dir>/lead_time_nse_gauge_f5_pooled.png
    <out-dir>/lead_time_nse_gauge_f5_by_regime.png

Run on server:
    python3 plot_lead_time_nse_gauge_f5.py
    python3 plot_lead_time_nse_gauge_f5.py --route-dir /path/to/routed --out-dir /path/to/out
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DA_COLOR  = "tab:purple"
OL_COLOR  = "tab:gray"

DEFAULT_ROUTE_DIR = "/mnt/disk2/suma_helen_poster/da_results/da_forecast_f5_rekrig_20pct/routed"
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"

WATERSHED_AREA_KM2 = 113.18
MM_H_TO_M3_S       = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 23:00:00")
STORM_THRESHOLD_M3S   = 20.0
LOWFLOW_THRESHOLD_M3S =  5.0


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_parquet(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_parquet(path)
    df['issue_time'] = pd.to_datetime(df['issue_time'])

    if 'q_gauge_m3s' in df.columns and 'member' in df.columns:
        return df[['issue_time', 'lead_hour', 'member', 'q_gauge_m3s']]

    member_cols = sorted(c for c in df.columns if c.startswith('member_'))
    if not member_cols:
        raise ValueError(f"Cannot identify member columns in {path}. "
                         f"Columns: {df.columns.tolist()}")
    keep = ['issue_time', 'lead_hour']
    long = df[keep + member_cols].melt(
        id_vars=keep, value_vars=member_cols,
        var_name='member', value_name='q_gauge_m3s')
    return long


def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    q_col    = next(c for c in df.columns if 'q' in c.lower() or 'flow' in c.lower()
                    or 'discharge' in c.lower())
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if 'mm' in q_col.lower():
        series = series * MM_H_TO_M3_S
        print(f"  Converted obs from mm/h ('{q_col}') → m³/s (×{MM_H_TO_M3_S:.4f})")
    return series


# ── NSE computation ───────────────────────────────────────────────────────────

def nse_by_lead(df_long, obs_series, issue_mask=None):
    """Return (leads, nse_values) pooling all issue times for each lead hour."""
    df = df_long.copy()
    if issue_mask is not None:
        df = df[df['issue_time'].isin(issue_mask)]
    if len(df) == 0:
        return None, None

    df['valid_time'] = df['issue_time'] + pd.to_timedelta(df['lead_hour'], unit='h')
    df['obs'] = df['valid_time'].map(obs_series)
    df = df.dropna(subset=['obs'])
    if len(df) == 0:
        return None, None

    ens_mean = (df.groupby(['issue_time', 'lead_hour'])['q_gauge_m3s']
                  .mean().reset_index(name='ens_mean'))
    obs_per  = (df.groupby(['issue_time', 'lead_hour'])['obs']
                  .first().reset_index(name='obs'))
    panel = ens_mean.merge(obs_per, on=['issue_time', 'lead_hour'])

    leads = sorted(panel['lead_hour'].unique())
    nse_vals = []
    for L in leads:
        sub = panel[panel['lead_hour'] == L]
        obs_v = sub['obs'].values
        sim_v = sub['ens_mean'].values
        denom = np.sum((obs_v - obs_v.mean()) ** 2)
        nse_L = 1.0 - np.sum((obs_v - sim_v) ** 2) / denom if denom > 0 else np.nan
        nse_vals.append(nse_L)

    return np.asarray(leads), np.asarray(nse_vals)


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_pooled(out_path, leads_da, nse_da, leads_ol, nse_ol):
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(leads_da, nse_da, color=DA_COLOR, lw=2.4, marker='o',
            zorder=4, label="F5 DA on (re-kriged σ²)")
    ax.plot(leads_ol, nse_ol, color=OL_COLOR, lw=2.4, marker='s',
            linestyle='--', zorder=4, label="Open loop (no DA)")

    ax.axhline(0, color='black', lw=0.8, linestyle=':', alpha=0.5)
    ax.set_xlabel("Forecast lead time (hours)", fontsize=11)
    ax.set_ylabel("NSE (ensemble-mean vs USGS obs)", fontsize=11)
    ax.set_title("Gauge-level lead-time NSE decay — F5 (re-kriged σ²), 20% holdout\n"
                 "USGS 03463300  |  all issue times pooled (Oct 2023 – Oct 2024)", fontsize=12)
    ax.set_xticks(np.arange(0, int(max(leads_da)) + 1, 1))
    ax.set_ylim(-0.1, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=10, loc='lower left', frameon=True, framealpha=0.92)

    plt.tight_layout()
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def plot_by_regime(out_path, regime_results):
    n = len(regime_results)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, (label, n_issue, res_da, res_ol) in zip(axes, regime_results):
        if res_da[0] is None:
            ax.set_title(f"{label}\n(no issue times)", fontsize=10)
            ax.axis('off')
            continue
        leads_da, nse_da = res_da
        leads_ol, nse_ol = res_ol
        ax.plot(leads_da, nse_da, color=DA_COLOR, lw=2.2, marker='o',
                label="F5 DA on")
        ax.plot(leads_ol, nse_ol, color=OL_COLOR, lw=2.2, marker='s',
                linestyle='--', label="Open loop")
        ax.axhline(0, color='black', lw=0.8, linestyle=':', alpha=0.5)
        ax.set_title(f"{label}\n({n_issue} issue times)", fontsize=10)
        ax.set_xlabel("Lead time (hours)", fontsize=10)
        ax.set_xticks(np.arange(0, int(max(leads_da)) + 1, 1))
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=9, loc='lower left', frameon=True, framealpha=0.9)

    axes[0].set_ylabel("NSE (ensemble-mean vs USGS obs)", fontsize=11)
    fig.suptitle("Lead-time NSE decay by regime — F5 (re-kriged σ²), USGS 03463300, 20% holdout",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

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

    print("Loading parquets...")
    df_da = load_parquet(os.path.join(args.route_dir, args.da_name))
    df_ol = load_parquet(os.path.join(args.route_dir, args.ol_name))
    obs   = load_usgs(args.usgs_csv)

    print(f"DA rows: {len(df_da):,}  |  OL rows: {len(df_ol):,}")

    # ---- Pooled NSE ----
    leads_da, nse_da = nse_by_lead(df_da, obs)
    leads_ol, nse_ol = nse_by_lead(df_ol, obs)
    plot_pooled(
        os.path.join(out_dir, "lead_time_nse_gauge_f5_pooled.png"),
        leads_da, nse_da, leads_ol, nse_ol)

    # ---- Regime split ----
    issue_times = pd.to_datetime(sorted(df_da['issue_time'].unique()))
    obs_at_issue = pd.Series(issue_times, index=issue_times).map(obs)

    helene_mask  = (issue_times >= HELENE_START) & (issue_times <= HELENE_END)
    storm_mask   = obs_at_issue > STORM_THRESHOLD_M3S
    lowflow_mask = obs_at_issue < LOWFLOW_THRESHOLD_M3S

    regime_results = []
    for label, mask in [
        (f"Helene window (Sep 24–28 2024)",                  helene_mask),
        (f"Storm hours (obs > {STORM_THRESHOLD_M3S:.0f} m³/s)",   storm_mask),
        (f"Low flow (obs < {LOWFLOW_THRESHOLD_M3S:.0f} m³/s)",    lowflow_mask),
    ]:
        kept = set(pd.to_datetime(issue_times[mask]))
        res_da = nse_by_lead(df_da, obs, kept)
        res_ol = nse_by_lead(df_ol, obs, kept)
        regime_results.append((label, len(kept), res_da, res_ol))

    plot_by_regime(
        os.path.join(out_dir, "lead_time_nse_gauge_f5_by_regime.png"),
        regime_results)

    print("Done.")


if __name__ == "__main__":
    main()
