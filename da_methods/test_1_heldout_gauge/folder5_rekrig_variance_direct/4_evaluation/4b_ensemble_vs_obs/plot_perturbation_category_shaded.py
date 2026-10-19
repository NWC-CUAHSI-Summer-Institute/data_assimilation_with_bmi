"""
Per-catchment shaded ensemble-band plot, organized by perturbation category.
F5 re-kriged variance direct (1 gauge holdout).

Three categories:
    1. Initial states            (red)
    2. Meteorological forcings   (blue)
    3. Hydrological states       (green)

For each category, all 20 ensemble members are shown as a shaded band
(min-max envelope fill) plus a thicker median line in the same color.
Re-kriged Qkrig observation overlaid in black. Hurricane Helene peak window
shaded in pink.

Inputs (existing per-source sensitivity CSVs from run_perturbation_sensitivity.py):
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_init.csv      (20 members, init only)
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_forcing.csv   (20 members, forcing only)
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_process.csv   (20 members, process noise only)
    <OBS-DIR>/<CAT>.csv                             (re-kriged obs, qkrig_mm_hr column)

Outputs:
    <SEN-DIR>/<CAT>/<CAT>_perturbation_categories_linear.png
    <SEN-DIR>/<CAT>/<CAT>_perturbation_categories_log.png
"""
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

CAT = "cat-1016300"

SEN_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder5_rekrig_variance_direct"
OBS_DIR = "/mnt/disk2/1400_sites_helene/catchment_ts_no_03463300_dynamic_variance_rekrig"

OUT_LINEAR = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_linear.png")
OUT_LOG    = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_log.png")

PLOT_START = pd.Timestamp("2024-09-20")
PLOT_END   = pd.Timestamp("2024-10-05")

HELENE_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 00:00:00")

CATEGORIES = [
    ("init",    "Initial states only",          "tab:red"),
    ("forcing", "Meteorological forcings only", "tab:blue"),
    ("process", "Hydrological states only",     "tab:green"),
]


def load_members(source):
    path = os.path.join(SEN_DIR, CAT, f"{CAT}_sensitivity_{source}.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    return df["date"].values, df[member_cols].to_numpy(dtype=float)


def load_obs():
    p = os.path.join(OBS_DIR, f"{CAT}.csv")
    if not os.path.exists(p):
        return None, None
    df = pd.read_csv(p)
    time_col = next((c for c in df.columns
                     if c.lower() in ("datetime", "date", "time", "timestamp")),
                    df.columns[0])
    df[time_col] = pd.to_datetime(df[time_col])
    # rekrig files use qkrig_mm_hr column
    q_col = next((c for c in df.columns
                  if "qkrig" in c.lower() and "var" not in c.lower()), None)
    if q_col is None:
        q_col = next(c for c in df.columns
                     if c != time_col and "var" not in c.lower()
                     and "step" not in c.lower()
                     and pd.api.types.is_numeric_dtype(df[c]))
    return df[time_col].values, df[q_col].values


def plot_panel(ax, obs_dates, obs_vals, log_y=False):
    for source, label, color in CATEGORIES:
        dates, q = load_members(source)
        if dates is None:
            continue
        d = pd.to_datetime(dates)
        mask = (d >= PLOT_START) & (d <= PLOT_END)
        if mask.sum() == 0:
            continue
        q_window = q[mask, :]
        qmin   = np.nanmin(q_window,    axis=1)
        qmax   = np.nanmax(q_window,    axis=1)
        median = np.nanmedian(q_window, axis=1)

        ax.fill_between(d[mask], qmin, qmax,
                        color=color, alpha=0.30, zorder=2, edgecolor="none")
        ax.plot(d[mask], median,
                color=color, lw=1.7, alpha=0.95, zorder=3,
                label=f"{label} (N=20)")

    ax.axvspan(HELENE_START, HELENE_END, color="salmon", alpha=0.15, zorder=1)

    if obs_dates is not None:
        od = pd.to_datetime(obs_dates)
        om = (od >= PLOT_START) & (od <= PLOT_END)
        ax.plot(od[om], obs_vals[om],
                color="black", lw=1.4, label="Qkrig re-kriged (obs)", zorder=4)

    ax.set_xlabel("Date", fontsize=10)
    if log_y:
        ax.set_yscale("log")
        ax.set_ylabel(r"q (mm hr$^{-1}$) [log scale]", fontsize=10)
        ax.set_ylim(0.01, None)
    else:
        ax.set_ylabel("Discharge (mm/h)", fontsize=10)
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper left", fontsize=9, frameon=True, framealpha=0.9)


def main():
    global CAT, SEN_DIR, OUT_LINEAR, OUT_LOG
    parser = argparse.ArgumentParser()
    parser.add_argument("--cat-id", default=CAT)
    args = parser.parse_args()
    CAT = args.cat_id
    OUT_LINEAR = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_linear.png")
    OUT_LOG    = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_log.png")

    obs_dates, obs_vals = load_obs()

    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, obs_dates, obs_vals, log_y=False)
    fig.suptitle(
        f"Ensemble forecast by perturbation category — {CAT} — F5 (re-kriged variance)\n"
        "Sep 20 - Oct 5, 2024 (Hurricane Helene window) | "
        "Shaded bands = min-max envelope across 20 members. Lines = ensemble median.",
        fontsize=11, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(OUT_LINEAR, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LINEAR}")

    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, obs_dates, obs_vals, log_y=True)
    fig.suptitle(
        f"Ensemble forecast by perturbation category — {CAT} — F5 (re-kriged variance) — log scale\n"
        "Sep 20 - Oct 5, 2024 (Hurricane Helene window) | "
        "Shaded bands = min-max envelope across 20 members. Lines = ensemble median.",
        fontsize=11, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(OUT_LOG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LOG}")


if __name__ == "__main__":
    main()
