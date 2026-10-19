"""
Lead-time decay curve, split by flow regime at issue time.

The pooled lead-time curve (plot_lead_time_decay.py) showed DA losing to
open-loop across most lead hours, but the test period is ~99% low-flow.
This script splits the same forecast CSVs by the flow regime at the issue
time t0, so we can see whether DA helps when it matters (storms / Helene)
and hurts when it doesn't (low flow).

Three regimes are partitioned on obs(t0) (kriging obs at the issue time):

  helene      issue_time ∈ [2024-09-24, 2024-09-28]   (the 5-day Helene window)
  storm       obs(t0) > 0.5 mm/h                       (high-flow issue times)
  low_flow    obs(t0) < 0.1 mm/h                       (low-flow issue times)

Same 2-panel layout per regime: RMSE vs lead (top), ensemble spread vs lead
(bottom). Three regime columns × 2 metric rows in one figure.

Inputs (from run_lead_time_forecast_sweep.py — no re-run needed):
    <leadtime>/<cat>/<cat>_lead_time_forecasts_da.csv
    <leadtime>/<cat>/<cat>_lead_time_forecasts_openloop.csv
Obs (from production EnKF):
    <da>/<cat>/<cat>_test_results.csv

Output:
    <leadtime>/<cat>/<cat>_lead_time_decay_by_regime.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DA_COLOR = "tab:purple"
OL_COLOR = "tab:gray"

DEFAULT_LEADTIME_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_lead_time_forecast"
DEFAULT_DA_DIR       = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"

# Module-level globals set by main() so helpers see them
CAT = "cat-1016300"
LEADTIME_DIR = DEFAULT_LEADTIME_DIR
DA_DIR       = DEFAULT_DA_DIR
OUT_PNG      = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_decay_by_regime.png")

# Regime definitions
HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 23:00:00")
STORM_OBS_THRESHOLD    = 0.5   # mm/h
LOWFLOW_OBS_THRESHOLD  = 0.1   # mm/h


def load_forecasts(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, parse_dates=['issue_time', 'valid_time'])
    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    return df, member_cols


def load_obs():
    obs_path = os.path.join(DA_DIR, CAT, f"{CAT}_test_results.csv")
    df = pd.read_csv(obs_path, parse_dates=['date'])
    return df.set_index('date')['obs_mm_h']


def regime_mask(issue_times, obs_at_issue, regime):
    """Boolean mask over issue_times for the named regime."""
    if regime == "helene":
        return (issue_times >= HELENE_START) & (issue_times <= HELENE_END)
    if regime == "storm":
        return obs_at_issue > STORM_OBS_THRESHOLD
    if regime == "low_flow":
        return obs_at_issue < LOWFLOW_OBS_THRESHOLD
    raise ValueError(regime)


def metrics_by_lead(df, member_cols, obs_series, issue_mask):
    """Compute (leads, rmse_mean, rmse_min, rmse_max, mean_std) restricted
    to issue times where issue_mask is True."""
    df = df.copy()
    df = df[df['issue_time'].isin(issue_mask)]
    if len(df) == 0:
        return None
    df['obs'] = df['valid_time'].map(obs_series)
    df = df.dropna(subset=['obs'])
    if len(df) == 0:
        return None
    df['ens_mean'] = df[member_cols].mean(axis=1)
    df['ens_std']  = df[member_cols].std(axis=1)

    leads = sorted(df['lead_hour'].unique())
    rmse_mean = np.full(len(leads), np.nan)
    rmse_lo   = np.full(len(leads), np.nan)
    rmse_hi   = np.full(len(leads), np.nan)
    std_mean  = np.full(len(leads), np.nan)
    for i, L in enumerate(leads):
        sub = df[df['lead_hour'] == L]
        if len(sub) == 0:
            continue
        err = sub['ens_mean'].values - sub['obs'].values
        rmse_mean[i] = float(np.sqrt(np.mean(err ** 2)))
        per_member = np.array([
            np.sqrt(np.mean((sub[c].values - sub['obs'].values) ** 2))
            for c in member_cols
        ])
        rmse_lo[i] = float(np.nanmin(per_member))
        rmse_hi[i] = float(np.nanmax(per_member))
        std_mean[i] = float(np.nanmean(sub['ens_std'].values))
    return np.asarray(leads), rmse_mean, rmse_lo, rmse_hi, std_mean


def plot_pair(ax_rmse, ax_std, da_metrics, ol_metrics,
              regime_label, n_issue_times):
    if da_metrics is None or ol_metrics is None:
        ax_rmse.set_title(f"{regime_label}\n(no issue times match)", fontsize=11)
        ax_rmse.axis('off'); ax_std.axis('off')
        return
    leads, rmse_da, rmse_da_lo, rmse_da_hi, std_da = da_metrics
    _,    rmse_ol, rmse_ol_lo, rmse_ol_hi, std_ol = ol_metrics

    # RMSE panel
    ax_rmse.fill_between(leads, rmse_da_lo, rmse_da_hi,
                         color=DA_COLOR, alpha=0.18, zorder=2)
    ax_rmse.fill_between(leads, rmse_ol_lo, rmse_ol_hi,
                         color=OL_COLOR, alpha=0.18, zorder=2)
    ax_rmse.plot(leads, rmse_da, color=DA_COLOR, lw=2.2, marker='o',
                 zorder=4, label="DA on")
    ax_rmse.plot(leads, rmse_ol, color=OL_COLOR, lw=2.2, marker='s',
                 linestyle='--', zorder=4, label="Open-loop")
    ax_rmse.set_title(f"{regime_label}\n({n_issue_times} issue times)", fontsize=11)
    ax_rmse.set_ylabel("RMSE (mm/h)", fontsize=10)
    ax_rmse.grid(True, alpha=0.25)
    ax_rmse.legend(fontsize=8, loc='upper left', frameon=True, framealpha=0.92)

    # Spread panel
    ax_std.plot(leads, std_da, color=DA_COLOR, lw=2.0, marker='o', label="DA")
    ax_std.plot(leads, std_ol, color=OL_COLOR, lw=2.0, marker='s',
                linestyle='--', label="Open-loop")
    ax_std.set_xlabel("Forecast lead time (hours)", fontsize=10)
    ax_std.set_ylabel("Mean ensemble std (mm/h)", fontsize=10)
    ax_std.set_xticks(np.arange(1, 19))
    ax_std.grid(True, alpha=0.25)
    ax_std.legend(fontsize=8, loc='upper left', frameon=True, framealpha=0.92)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',        default="cat-1016300")
    parser.add_argument('--leadtime-dir',  default=DEFAULT_LEADTIME_DIR,
                        help='Dir holding <cat>/<cat>_lead_time_forecasts_{da,openloop}.csv')
    parser.add_argument('--da-dir',        default=DEFAULT_DA_DIR,
                        help='Dir holding <cat>/<cat>_test_results.csv (for obs)')
    args = parser.parse_args()

    global CAT, LEADTIME_DIR, DA_DIR, OUT_PNG
    CAT = args.cat_id
    LEADTIME_DIR = args.leadtime_dir
    DA_DIR       = args.da_dir
    OUT_PNG      = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_decay_by_regime.png")

    da_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_da.csv")
    ol_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_openloop.csv")
    df_da, m_cols = load_forecasts(da_path)
    df_ol, _      = load_forecasts(ol_path)
    obs_series = load_obs()

    # Unique issue times across the run, with obs at each
    issue_times = pd.to_datetime(sorted(df_da['issue_time'].unique()))
    obs_at_issue = pd.Series(issue_times, index=issue_times).map(obs_series)

    regimes = [
        ("helene",    "Helene window (Sep 24–28 2024)"),
        ("storm",     f"Storm hours (obs(t0) > {STORM_OBS_THRESHOLD} mm/h)"),
        ("low_flow",  f"Low flow (obs(t0) < {LOWFLOW_OBS_THRESHOLD} mm/h)"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10),
                              gridspec_kw={"height_ratios": [1.4, 1.0]})
    for col, (key, label) in enumerate(regimes):
        mask = regime_mask(issue_times, obs_at_issue, key)
        kept = issue_times[mask]
        kept_set = set(pd.to_datetime(kept))
        da_m = metrics_by_lead(df_da, m_cols, obs_series, kept_set)
        ol_m = metrics_by_lead(df_ol, m_cols, obs_series, kept_set)
        plot_pair(axes[0, col], axes[1, col], da_m, ol_m, label, len(kept))

    fig.suptitle(
        f"Lead-time error decay by flow regime — {CAT}\n"
        f"Same forecast CSVs as the pooled view, partitioned on obs(t0).",
        fontsize=13, y=0.995,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
