"""
Catchment-level lead-time forecast error decay curve — F5 (re-kriged variance).

Two-panel figure for ONE catchment:
  TOP    — RMSE of ensemble-mean DA forecast vs Qkrig obs at each lead hour (1..18),
           with a shaded band for min/max per-member RMSE. Open-loop dashed.
  BOTTOM — Mean ensemble spread (std across 20 members) at each lead hour.

Inputs (from run_lead_time_forecast_sweep.py):
    <leadtime-dir>/<cat>/<cat>_lead_time_forecasts_da.csv
    <leadtime-dir>/<cat>/<cat>_lead_time_forecasts_openloop.csv
Obs:
    <da-dir>/<cat>/<cat>_test_results.csv  (obs_mm_h column)

Output:
    <leadtime-dir>/<cat>/<cat>_lead_time_decay.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DA_COLOR = "tab:purple"
OL_COLOR = "tab:gray"

F5_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct"
DEFAULT_KRIG_OBS_DIR = "/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance_rekrig"

DEFAULT_LEADTIME_DIR = F5_DIR
DEFAULT_DA_DIR       = F5_DIR

CAT          = "cat-1016300"
LEADTIME_DIR = DEFAULT_LEADTIME_DIR
DA_DIR       = DEFAULT_DA_DIR
KRIG_OBS_DIR = DEFAULT_KRIG_OBS_DIR
OUT_PNG      = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_decay.png")


def load_forecasts(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, parse_dates=['issue_time', 'valid_time'])
    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    return df, member_cols


def load_obs():
    # Prefer _test_results.csv if it exists; otherwise read from kriged obs file.
    test_path = os.path.join(DA_DIR, CAT, f"{CAT}_test_results.csv")
    if os.path.exists(test_path):
        df = pd.read_csv(test_path, parse_dates=['date'])
        return df.set_index('date')['obs_mm_h']
    krig_path = os.path.join(KRIG_OBS_DIR, f"{CAT}.csv")
    if not os.path.exists(krig_path):
        raise FileNotFoundError(
            f"No obs found: tried {test_path} and {krig_path}. "
            f"Pass --krig-obs-dir to the re-kriged obs directory.")
    df = pd.read_csv(krig_path)
    print(f"  kriged obs file columns: {list(df.columns)}")
    t_col = next(c for c in df.columns
                 if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    # Prefer columns that look like discharge obs (mm/h); exclude step/index columns
    obs_col = next(
        (c for c in df.columns
         if c != t_col and 'var' not in c.lower()
         and any(k in c.lower() for k in ('obs', 'q', 'krig', 'mm', 'flow', 'discharge'))
         and pd.api.types.is_numeric_dtype(df[c])),
        None,
    )
    if obs_col is None:
        # Fall back: any numeric column that isn't time-like or variance
        obs_col = next(c for c in df.columns
                       if c != t_col and 'var' not in c.lower()
                       and 'step' not in c.lower() and 'time' not in c.lower()
                       and 'index' not in c.lower()
                       and pd.api.types.is_numeric_dtype(df[c]))
    df[t_col] = pd.to_datetime(df[t_col])
    print(f"  obs from kriged file, column '{obs_col}'")
    return df.set_index(t_col)[obs_col].astype(float)


def metrics_by_lead(df, member_cols, obs_series):
    df = df.copy()
    df['obs'] = df['valid_time'].map(obs_series)
    df = df.dropna(subset=['obs'])
    df['ens_mean'] = df[member_cols].mean(axis=1)
    df['ens_std']  = df[member_cols].std(axis=1)

    leads = sorted(df['lead_hour'].unique())
    rmse_mean       = np.full(len(leads), np.nan)
    rmse_min_member = np.full(len(leads), np.nan)
    rmse_max_member = np.full(len(leads), np.nan)
    mean_ens_std    = np.full(len(leads), np.nan)
    for i, L in enumerate(leads):
        sub = df[df['lead_hour'] == L]
        if len(sub) == 0:
            continue
        err = sub['ens_mean'].values - sub['obs'].values
        rmse_mean[i] = float(np.sqrt(np.mean(err ** 2)))
        per_member_rmse = np.array([
            np.sqrt(np.mean((sub[c].values - sub['obs'].values) ** 2))
            for c in member_cols
        ])
        rmse_min_member[i] = float(np.nanmin(per_member_rmse))
        rmse_max_member[i] = float(np.nanmax(per_member_rmse))
        mean_ens_std[i]    = float(np.nanmean(sub['ens_std'].values))
    return (np.asarray(leads), rmse_mean, rmse_min_member,
            rmse_max_member, mean_ens_std)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',       default="cat-1016300")
    parser.add_argument('--leadtime-dir',  default=DEFAULT_LEADTIME_DIR,
                        help='Dir holding <cat>/<cat>_lead_time_forecasts_{da,openloop}.csv')
    parser.add_argument('--da-dir',        default=DEFAULT_DA_DIR,
                        help='Dir holding <cat>/<cat>_test_results.csv (for obs_mm_h)')
    parser.add_argument('--krig-obs-dir',  default=DEFAULT_KRIG_OBS_DIR,
                        help='Fallback: dir holding per-catchment kriged obs CSVs '
                             '(used when _test_results.csv is absent)')
    args = parser.parse_args()

    global CAT, LEADTIME_DIR, DA_DIR, KRIG_OBS_DIR, OUT_PNG
    CAT          = args.cat_id
    LEADTIME_DIR = args.leadtime_dir
    DA_DIR       = args.da_dir
    KRIG_OBS_DIR = args.krig_obs_dir
    OUT_PNG      = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_decay.png")

    da_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_da.csv")
    ol_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_openloop.csv")
    df_da, m_cols = load_forecasts(da_path)
    df_ol, _      = load_forecasts(ol_path)
    obs_series = load_obs()

    leads_da, rmse_da, rmse_da_lo, rmse_da_hi, std_da = metrics_by_lead(
        df_da, m_cols, obs_series)
    leads_ol, rmse_ol, rmse_ol_lo, rmse_ol_hi, std_ol = metrics_by_lead(
        df_ol, m_cols, obs_series)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                          gridspec_kw={"height_ratios": [1.4, 1.0]})

    ax_top.fill_between(leads_da, rmse_da_lo, rmse_da_hi,
                        color=DA_COLOR, alpha=0.18, zorder=2,
                        label="DA — member spread (min/max RMSE)")
    ax_top.fill_between(leads_ol, rmse_ol_lo, rmse_ol_hi,
                        color=OL_COLOR, alpha=0.18, zorder=2,
                        label="Open-loop — member spread")
    ax_top.plot(leads_da, rmse_da, color=DA_COLOR, lw=2.4, marker='o',
                zorder=4, label="DA on (ensemble-mean RMSE)")
    ax_top.plot(leads_ol, rmse_ol, color=OL_COLOR, lw=2.4, marker='s',
                linestyle='--', zorder=4, label="Open-loop (ensemble-mean RMSE)")
    ax_top.set_ylabel("RMSE vs Qkrig obs (mm/h)", fontsize=11)
    ax_top.set_title(
        f"Forecast lead-time error decay — {CAT}\n"
        f"F5 (re-kriged variance direct) — test period Oct 2023 – Oct 2024",
        fontsize=12,
    )
    ax_top.grid(True, alpha=0.25)
    ax_top.legend(fontsize=9, loc='upper left', frameon=True, framealpha=0.92)

    ax_bot.plot(leads_da, std_da, color=DA_COLOR, lw=2.2, marker='o',
                label="DA — mean ensemble std")
    ax_bot.plot(leads_ol, std_ol, color=OL_COLOR, lw=2.2, marker='s',
                linestyle='--', label="Open-loop — mean ensemble std")
    ax_bot.set_xlabel("Forecast lead time (hours)", fontsize=11)
    ax_bot.set_ylabel("Mean ensemble std (mm/h)", fontsize=11)
    ax_bot.set_xticks(np.arange(1, 19))
    ax_bot.grid(True, alpha=0.25)
    ax_bot.legend(fontsize=9, loc='upper left', frameon=True, framealpha=0.92)

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
