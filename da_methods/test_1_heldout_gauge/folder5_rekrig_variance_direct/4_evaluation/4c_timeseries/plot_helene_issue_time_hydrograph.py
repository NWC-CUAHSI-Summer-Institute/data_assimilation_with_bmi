"""
Per-issue-time forecast hydrograph diagnostic — F5 (re-kriged variance).

For two contrasting issue times — one storm-peak (Helene) and one low-flow —
show what the DA and open-loop ensembles actually produce over the 18-hour
forecast window, alongside the kriging observation.

Each panel shows:
  - All 20 DA members as faint purple lines + median (thick purple)
  - All 20 open-loop members as faint gray lines + median (thick gray, dashed)
  - Kriging obs as a thick black line (context: t0-6 through t0+18)
  - Vertical line at t0 + annotation with obs(t0) and lead-1 ensemble means

Inputs (from batch_run_f5_lead_time_sweep.sh):
    <leadtime>/<cat>/<cat>_lead_time_forecasts_da.csv
    <leadtime>/<cat>/<cat>_lead_time_forecasts_openloop.csv
Obs:
    <leadtime>/<cat>/<cat>_test_results.csv  if present, else
    <krig-obs-dir>/<cat>.csv  (qkrig_mm_hr column)

Output:
    <leadtime>/<cat>/<cat>_issue_time_hydrograph_helene_vs_lowflow.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

F5_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct"
DEFAULT_KRIG_OBS_DIR = "/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance_rekrig"

DEFAULT_LEADTIME_DIR = F5_DIR
DEFAULT_DA_DIR       = F5_DIR

CAT          = "cat-1016300"
LEADTIME_DIR = DEFAULT_LEADTIME_DIR
DA_DIR       = DEFAULT_DA_DIR
KRIG_OBS_DIR = DEFAULT_KRIG_OBS_DIR

DA_COLOR  = "tab:purple"
OL_COLOR  = "tab:gray"
OBS_COLOR = "black"

HELENE_T0  = pd.Timestamp("2024-09-26 12:00:00")
LOWFLOW_T0 = pd.Timestamp("2024-03-15 12:00:00")

CONTEXT_HOURS_BEFORE = 6
LEAD_HOURS_AFTER     = 18


def load_forecasts(path):
    df = pd.read_csv(path, parse_dates=['issue_time', 'valid_time'])
    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    return df, member_cols


def load_obs():
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
    t_col = next(c for c in df.columns
                 if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    obs_col = next(
        (c for c in df.columns
         if c != t_col and 'var' not in c.lower()
         and any(k in c.lower() for k in ('obs', 'q', 'krig', 'mm', 'flow', 'discharge'))
         and pd.api.types.is_numeric_dtype(df[c])),
        None,
    )
    if obs_col is None:
        obs_col = next(c for c in df.columns
                       if c != t_col and 'var' not in c.lower()
                       and 'step' not in c.lower() and 'time' not in c.lower()
                       and 'index' not in c.lower()
                       and pd.api.types.is_numeric_dtype(df[c]))
    df[t_col] = pd.to_datetime(df[t_col])
    return df.set_index(t_col)[obs_col].astype(float)


def slice_forecast(df, member_cols, t0):
    sub = df[df['issue_time'] == t0].sort_values('lead_hour')
    if len(sub) == 0:
        raise ValueError(f"No forecast rows for issue_time={t0}")
    valid_times = pd.to_datetime(sub['valid_time'].values)
    member_arr  = sub[member_cols].to_numpy(dtype=float)
    return valid_times, member_arr


def plot_panel(ax, t0, df_da, df_ol, member_cols, obs_series, title_extra=""):
    vt_da, q_da = slice_forecast(df_da, member_cols, t0)
    vt_ol, q_ol = slice_forecast(df_ol, member_cols, t0)

    ctx_start  = t0 - pd.Timedelta(hours=CONTEXT_HOURS_BEFORE)
    ctx_end    = t0 + pd.Timedelta(hours=LEAD_HOURS_AFTER)
    obs_window = obs_series.loc[ctx_start:ctx_end]

    for j in range(q_ol.shape[1]):
        ax.plot(vt_ol, q_ol[:, j], color=OL_COLOR, lw=0.6, alpha=0.45, zorder=2)
    for j in range(q_da.shape[1]):
        ax.plot(vt_da, q_da[:, j], color=DA_COLOR, lw=0.6, alpha=0.45, zorder=3)
    ax.plot(vt_ol, np.nanmedian(q_ol, axis=1), color=OL_COLOR, lw=2.6,
            linestyle="--", zorder=4, label="Open-loop median")
    ax.plot(vt_da, np.nanmedian(q_da, axis=1), color=DA_COLOR, lw=2.6,
            zorder=5, label="DA median")
    ax.plot(obs_window.index, obs_window.values, color=OBS_COLOR, lw=1.8,
            zorder=6, label="Qkrig (obs)")

    ax.axvline(t0, color="black", lw=1.0, linestyle=":", alpha=0.6, zorder=1)
    obs_t0        = float(obs_series.get(t0, np.nan))
    obs_lead1     = float(obs_series.get(t0 + pd.Timedelta(hours=1), np.nan))
    da_lead1_mean = float(np.nanmean(q_da[0, :]))
    ol_lead1_mean = float(np.nanmean(q_ol[0, :]))

    info_text = (
        f"t0 = {t0.strftime('%Y-%m-%d %H:%M')}\n"
        f"obs(t0)       = {obs_t0:.3f} mm/h\n"
        f"obs(t0+1)     = {obs_lead1:.3f} mm/h\n"
        f"DA mean(t0+1) = {da_lead1_mean:.3f} mm/h\n"
        f"OL mean(t0+1) = {ol_lead1_mean:.3f} mm/h"
    )
    ax.text(0.02, 0.97, info_text,
            transform=ax.transAxes,
            fontsize=9, va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      alpha=0.85, edgecolor='gray'))

    ax.set_title(title_extra, fontsize=12)
    ax.set_ylabel("Discharge (mm/h)", fontsize=10)
    ax.set_xlabel("Date", fontsize=10)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    ax.tick_params(axis='x', rotation=30, labelsize=8)
    ax.grid(True, alpha=0.2)
    ax.legend(loc='upper right', fontsize=9, frameon=True, framealpha=0.9)


def snap_to_nearest(t0, available_times):
    times = pd.DatetimeIndex(available_times)
    idx = np.abs((times - t0).total_seconds()).argmin()
    snapped = times[idx]
    if snapped != t0:
        print(f"  Snapped {t0} -> {snapped} (nearest available issue time)")
    return snapped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',        default="cat-1016300")
    parser.add_argument('--leadtime-dir',  default=DEFAULT_LEADTIME_DIR)
    parser.add_argument('--da-dir',        default=DEFAULT_DA_DIR)
    parser.add_argument('--krig-obs-dir',  default=DEFAULT_KRIG_OBS_DIR,
                        help='Dir holding per-catchment re-kriged obs CSVs '
                             '(used when _test_results.csv is absent)')
    parser.add_argument('--helene-t0',     default=None)
    parser.add_argument('--lowflow-t0',    default=None)
    args = parser.parse_args()

    global CAT, LEADTIME_DIR, DA_DIR, KRIG_OBS_DIR
    CAT          = args.cat_id
    LEADTIME_DIR = args.leadtime_dir
    DA_DIR       = args.da_dir
    KRIG_OBS_DIR = args.krig_obs_dir
    out_png = os.path.join(LEADTIME_DIR, CAT,
                           f"{CAT}_issue_time_hydrograph_helene_vs_lowflow.png")

    da_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_da.csv")
    ol_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_openloop.csv")
    df_da, m_cols = load_forecasts(da_path)
    df_ol, _      = load_forecasts(ol_path)
    obs_series = load_obs()

    available  = df_da['issue_time'].unique()
    helene_t0  = snap_to_nearest(
        pd.Timestamp(args.helene_t0)  if args.helene_t0  else HELENE_T0,  available)
    lowflow_t0 = snap_to_nearest(
        pd.Timestamp(args.lowflow_t0) if args.lowflow_t0 else pd.Timestamp("2024-03-15 00:00:00"),
        available)

    fig, (ax_h, ax_l) = plt.subplots(1, 2, figsize=(20, 7))

    plot_panel(ax_h, helene_t0, df_da, df_ol, m_cols, obs_series,
               title_extra=f"Helene peak issue time — {helene_t0.strftime('%Y-%m-%d %H:%M')}")
    plot_panel(ax_l, lowflow_t0, df_da, df_ol, m_cols, obs_series,
               title_extra=f"Low-flow issue time — {lowflow_t0.strftime('%Y-%m-%d %H:%M')}")

    fig.suptitle(
        f"Per-issue-time forecast hydrograph — {CAT} — F5 (re-kriged variance)\n"
        "Faint lines = 20 ensemble members.  Thick lines = ensemble medians.  "
        "Vertical dotted line = t0 (DA stops here for both scenarios).",
        fontsize=12, y=1.00,
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
