"""
Per-issue-time forecast hydrograph diagnostic.

For two contrasting issue times — one storm-peak (Helene) and one low-flow —
show what the DA and open-loop ensembles actually produce over the 18-hour
forecast window, alongside the kriging observation. Designed to answer in one
picture: is DA over-shooting, collapsing, or oscillating compared to open-loop?

Each panel shows:
  - All 20 DA members as faint purple lines + median (thick purple)
  - All 20 open-loop members as faint gray lines + median (thick gray, dashed)
  - Kriging obs as a thick black line (context: t0-6 through t0+18)
  - Vertical line at t0 + annotation with obs(t0) and lead-1 ensemble means

Inputs (from run_lead_time_forecast_sweep.py):
    <leadtime>/<cat>/<cat>_lead_time_forecasts_da.csv
    <leadtime>/<cat>/<cat>_lead_time_forecasts_openloop.csv
Obs (from production EnKF):
    <da>/<cat>/<cat>_test_results.csv

Output:
    <leadtime>/<cat>/<cat>_issue_time_hydrograph_helene_vs_lowflow.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

DEFAULT_LEADTIME_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_lead_time_forecast"
DEFAULT_DA_DIR       = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"

# Module-level globals set by main()
CAT = "cat-1016300"
LEADTIME_DIR = DEFAULT_LEADTIME_DIR
DA_DIR       = DEFAULT_DA_DIR

DA_COLOR = "tab:purple"
OL_COLOR = "tab:gray"
OBS_COLOR = "black"

# Two issue times — pick a Helene-peak one and a typical low-flow one
HELENE_T0  = pd.Timestamp("2024-09-26 12:00:00")
LOWFLOW_T0 = pd.Timestamp("2024-03-15 12:00:00")

CONTEXT_HOURS_BEFORE = 6   # hours of obs context shown before t0
LEAD_HOURS_AFTER    = 18   # forecast horizon

OUT_PNG = os.path.join(
    LEADTIME_DIR, CAT,
    f"{CAT}_issue_time_hydrograph_helene_vs_lowflow.png",
)


def load_forecasts(path):
    df = pd.read_csv(path, parse_dates=['issue_time', 'valid_time'])
    member_cols = sorted([c for c in df.columns if c.startswith('member_')])
    return df, member_cols


def load_obs():
    obs_path = os.path.join(DA_DIR, CAT, f"{CAT}_test_results.csv")
    df = pd.read_csv(obs_path, parse_dates=['date'])
    return df.set_index('date')['obs_mm_h']


def slice_forecast(df, member_cols, t0):
    """Return (valid_times, member_array shape (lead, N)) for a single issue time."""
    sub = df[df['issue_time'] == t0].sort_values('lead_hour')
    if len(sub) == 0:
        raise ValueError(f"No forecast rows for issue_time={t0}")
    valid_times = pd.to_datetime(sub['valid_time'].values)
    member_arr  = sub[member_cols].to_numpy(dtype=float)
    return valid_times, member_arr


def plot_panel(ax, t0, df_da, df_ol, member_cols, obs_series, title_extra=""):
    # Forecast trajectories
    vt_da, q_da = slice_forecast(df_da, member_cols, t0)
    vt_ol, q_ol = slice_forecast(df_ol, member_cols, t0)

    # Obs context window
    ctx_start = t0 - pd.Timedelta(hours=CONTEXT_HOURS_BEFORE)
    ctx_end   = t0 + pd.Timedelta(hours=LEAD_HOURS_AFTER)
    obs_window = obs_series.loc[ctx_start:ctx_end]

    # Open-loop members (drawn first so DA paints on top)
    for j in range(q_ol.shape[1]):
        ax.plot(vt_ol, q_ol[:, j], color=OL_COLOR, lw=0.6, alpha=0.45, zorder=2)
    # DA members
    for j in range(q_da.shape[1]):
        ax.plot(vt_da, q_da[:, j], color=DA_COLOR, lw=0.6, alpha=0.45, zorder=3)
    # Medians
    ax.plot(vt_ol, np.nanmedian(q_ol, axis=1), color=OL_COLOR, lw=2.6,
            linestyle="--", zorder=4, label="Open-loop median")
    ax.plot(vt_da, np.nanmedian(q_da, axis=1), color=DA_COLOR, lw=2.6,
            zorder=5, label="DA median")
    # Obs
    ax.plot(obs_window.index, obs_window.values, color=OBS_COLOR, lw=1.8,
            zorder=6, label="Qkrig (obs)")

    # t0 marker + annotation
    ax.axvline(t0, color="black", lw=1.0, linestyle=":", alpha=0.6, zorder=1)
    obs_t0 = float(obs_series.get(t0, np.nan))
    obs_lead1 = float(obs_series.get(t0 + pd.Timedelta(hours=1), np.nan))
    da_lead1_mean = float(np.nanmean(q_da[0, :]))
    ol_lead1_mean = float(np.nanmean(q_ol[0, :]))

    info_text = (
        f"t0 = {t0.strftime('%Y-%m-%d %H:%M')}\n"
        f"obs(t0)    = {obs_t0:.3f} mm/h\n"
        f"obs(t0+1)  = {obs_lead1:.3f} mm/h\n"
        f"DA mean(t0+1)  = {da_lead1_mean:.3f} mm/h\n"
        f"OL mean(t0+1)  = {ol_lead1_mean:.3f} mm/h"
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
    """Return the available issue_time closest to t0."""
    times = pd.DatetimeIndex(available_times)
    idx = np.abs((times - t0).total_seconds()).argmin()
    snapped = times[idx]
    if snapped != t0:
        print(f"  Snapped {t0} -> {snapped} (nearest available issue time)")
    return snapped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id',        default="cat-1016300")
    parser.add_argument('--leadtime-dir',  default=DEFAULT_LEADTIME_DIR,
                        help='Dir holding <cat>/<cat>_lead_time_forecasts_{da,openloop}.csv')
    parser.add_argument('--da-dir',        default=DEFAULT_DA_DIR,
                        help='Dir holding <cat>/<cat>_test_results.csv (for obs)')
    parser.add_argument('--helene-t0',     default=None,
                        help='Issue time for the Helene panel (default: 2024-09-26 12:00:00, '
                             'snapped to nearest available)')
    parser.add_argument('--lowflow-t0',    default=None,
                        help='Issue time for the low-flow panel (default: 2024-03-15 00:00:00, '
                             'snapped to nearest available)')
    args = parser.parse_args()

    global CAT, LEADTIME_DIR, DA_DIR, OUT_PNG
    CAT = args.cat_id
    LEADTIME_DIR = args.leadtime_dir
    DA_DIR       = args.da_dir
    OUT_PNG      = os.path.join(LEADTIME_DIR, CAT,
                                f"{CAT}_issue_time_hydrograph_helene_vs_lowflow.png")

    da_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_da.csv")
    ol_path = os.path.join(LEADTIME_DIR, CAT, f"{CAT}_lead_time_forecasts_openloop.csv")
    df_da, m_cols = load_forecasts(da_path)
    df_ol, _      = load_forecasts(ol_path)
    obs_series = load_obs()

    available = df_da['issue_time'].unique()
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
        f"Per-issue-time forecast hydrograph — {CAT}\n"
        "Faint lines = 20 ensemble members.  Thick lines = ensemble medians.  "
        "Vertical dotted line = t0 (DA stops here for both scenarios).",
        fontsize=12, y=1.00,
    )
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
