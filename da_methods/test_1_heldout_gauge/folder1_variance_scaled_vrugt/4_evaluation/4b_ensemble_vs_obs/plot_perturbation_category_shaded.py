"""
Per-catchment shaded ensemble-band plot, organized by perturbation category.

Three categories:
    1. Initial states            (red)
    2. Meteorological forcings   (blue)
    3. Hydrological states       (green)

For each category, all 20 ensemble members are shown as a shaded band
(min-max envelope fill) plus a thicker median line in the same color.
Qkrig observation overlaid in black. Hurricane Helene peak window shaded
in pink. Styled after a standard ensemble-forecast figure layout.

Inputs (existing per-source sensitivity CSVs from run_perturbation_sensitivity.py):
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_init.csv      (20 members, init only)
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_forcing.csv   (20 members, forcing only)
    <SEN-DIR>/<CAT>/<CAT>_sensitivity_process.csv   (20 members, process noise only)
    <DA-DIR>/<CAT>/<CAT>_test_results.csv           (Qkrig obs)

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
import matplotlib.patches as mpatches

CAT = "cat-1016300"

SEN_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
OBS_DIR = "/home/svyas/catchment_ts_no_03463300_gapfilled"

OUT_LINEAR = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_linear.png")
OUT_LOG    = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_log.png")

# Plot window — wider context, similar to the paper's Sep 10 - Oct 08
PLOT_START = pd.Timestamp("2024-09-20")
PLOT_END   = pd.Timestamp("2024-10-05")

# Helene peak band
HELENE_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 00:00:00")

# Category configuration: file suffix, display label, color
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
    time_col = 'datetime' if 'datetime' in df.columns else 'date'
    df[time_col] = pd.to_datetime(df[time_col])
    return df[time_col].values, df["qkrig"].values


def plot_panel(ax, obs_dates, obs_vals, log_y=False):
    handles_labels = []  # for the legend

    # Plot each category as a shaded band + median line
    for source, label, color in CATEGORIES:
        dates, q = load_members(source)
        if dates is None:
            continue
        d = pd.to_datetime(dates)
        mask = (d >= PLOT_START) & (d <= PLOT_END)
        if mask.sum() == 0:
            continue
        # Min-max envelope across all 20 members per timestep (widest possible band).
        # Bands are visually narrow even with min/max because perturbations are
        # tuned for production EnKF stability, not for max visible spread.
        q_window = q[mask, :]
        qmin   = np.nanmin(q_window,    axis=1)
        qmax   = np.nanmax(q_window,    axis=1)
        median = np.nanmedian(q_window, axis=1)

        ax.fill_between(d[mask], qmin, qmax,
                        color=color, alpha=0.30, zorder=2,
                        edgecolor="none")
        line, = ax.plot(d[mask], median,
                        color=color, lw=1.7, alpha=0.95, zorder=3,
                        label=f"{label} (N=20)")
        handles_labels.append((line, label))

    # Helene peak shaded band (vertical)
    ax.axvspan(HELENE_START, HELENE_END,
               color="salmon", alpha=0.15, zorder=1)
    ax.text((HELENE_START + (HELENE_END - HELENE_START) / 2),
            ax.get_ylim()[1] if not log_y else 1.0,
            "Helene peak",
            fontsize=9, color="salmon",
            ha="center", va="bottom", zorder=3)

    # Observation
    if obs_dates is not None:
        od = pd.to_datetime(obs_dates)
        om = (od >= PLOT_START) & (od <= PLOT_END)
        ax.plot(od[om], obs_vals[om],
                color="black", lw=1.4, label="Qkrig (obs)", zorder=4)

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
    global CAT, OUT_LINEAR, OUT_LOG
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id', default=CAT)
    args = parser.parse_args()
    CAT = args.cat_id
    OUT_LINEAR = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_linear.png")
    OUT_LOG    = os.path.join(SEN_DIR, CAT, f"{CAT}_perturbation_categories_log.png")

    obs_dates, obs_vals = load_obs()

    # ----- Linear-y -----
    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, obs_dates, obs_vals, log_y=False)
    fig.suptitle(
        f"Ensemble forecast by perturbation category - {CAT} - "
        f"Sep 20 - Oct 5, 2024 (Hurricane Helene window)\n"
        "Shaded bands = min-max envelope across 20 members.   "
        "Lines = ensemble median.   DA off in all three sub-experiments.",
        fontsize=11, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(OUT_LINEAR, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LINEAR}")

    # ----- Log-y (paper style) -----
    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, obs_dates, obs_vals, log_y=True)
    fig.suptitle(
        f"Ensemble forecast by perturbation category - {CAT} - log-scale q - "
        f"Sep 20 - Oct 5, 2024 (Hurricane Helene window)\n"
        "Shaded bands = min-max envelope across 20 members.   "
        "Lines = ensemble median.   DA off in all three sub-experiments.",
        fontsize=11, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(OUT_LOG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LOG}")


if __name__ == "__main__":
    main()
