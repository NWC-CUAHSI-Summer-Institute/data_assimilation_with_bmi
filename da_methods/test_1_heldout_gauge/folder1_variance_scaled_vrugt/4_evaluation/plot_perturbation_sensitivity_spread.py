"""
Plot ensemble-spread time series, one curve per perturbation source.

For each catchment and each of the three sensitivity sub-experiments, compute the
hourly std-dev of streamflow across the 20 members. Plot all three as colored
lines on the same axes so the relative contribution of each source is directly
comparable hour by hour.

This is the cleanest answer to "which perturbation source contributes most to
ensemble spread" — the spaghetti version is hard to read when the bundles
overlap. A single std-dev curve per source removes the visual clutter.

Color legend:
    red   = std-dev across initial-state-only ensemble
    blue  = std-dev across forcing-only ensemble
    green = std-dev across process-noise-only ensemble

Produces:
    helene_sensitivity_spread_main.png       (3-catchment subset for the main figure)
    helene_sensitivity_spread_appendix.png   (full 3 x 7 grid of all 21 catchments)

Consumes the same per-source CSVs produced by run_perturbation_sensitivity.py.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

ALL_CATS = [
    "cat-1016279", "cat-1016280", "cat-1016281", "cat-1016282", "cat-1016283",
    "cat-1016300", "cat-1016301", "cat-1016302", "cat-1016303", "cat-1016304",
    "cat-1016305", "cat-1016306", "cat-1016307", "cat-1016308", "cat-1016309",
    "cat-1016310", "cat-1016311", "cat-1016312", "cat-1016313", "cat-1016314",
    "cat-1016315",
]

MAIN_CATS = ["cat-1016311", "cat-1016300", "cat-1016302"]

DA_DIR  = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
SEN_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
OBS_DIR = "/home/svyas/catchment_ts_no_03463300_gapfilled"
MAIN_PNG     = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt/helene_sensitivity_spread_main.png"
APPENDIX_PNG = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt/helene_sensitivity_spread_appendix.png"

ZOOM_START = pd.Timestamp("2024-09-24")
ZOOM_END   = pd.Timestamp("2024-09-28")

SOURCE_COLORS = {
    "init":    "tab:red",
    "forcing": "tab:blue",
    "process": "tab:green",
}
SOURCE_LABELS = {
    "init":    "Initial state perturbation",
    "forcing": "Forcing perturbation (P, PET)",
    "process": "Process noise on states",
}


def load_member_spread(cat, source):
    """Returns (dates, std_per_hour) for the test period.

    std_per_hour is the std-dev of the 20 member streamflow values at each hour.
    """
    path = os.path.join(SEN_DIR, cat, f"{cat}_sensitivity_{source}.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    member_cols = [c for c in df.columns if c.startswith("member_")]
    member_arr = df[member_cols].to_numpy(dtype=float)
    # Hourly std-dev across members (sample std, ddof=1 to match Pyy convention)
    std = member_arr.std(axis=1, ddof=1)
    return df["date"].values, std


def load_qkrig_obs(cat):
    obs_p = os.path.join(OBS_DIR, f"{cat}.csv")
    if os.path.exists(obs_p):
        df = pd.read_csv(obs_p)
        time_col = 'datetime' if 'datetime' in df.columns else 'date'
        df[time_col] = pd.to_datetime(df[time_col])
        return df[time_col].values, df["qkrig"].values
    path = os.path.join(DA_DIR, cat, f"{cat}_test_results.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    return df["date"].values, df["obs_mm_h"].values


def plot_one_panel(ax, cat, show_obs=True):
    """Render std-dev time series for the three sources plus optional obs reference."""
    for source, color in SOURCE_COLORS.items():
        dates, std = load_member_spread(cat, source)
        if dates is None:
            continue
        d = pd.to_datetime(dates)
        mask = (d >= ZOOM_START) & (d <= ZOOM_END)
        if mask.sum() == 0:
            continue
        ax.plot(d[mask], std[mask], color=color, lw=1.6,
                label=SOURCE_LABELS[source], zorder=2)

    if show_obs:
        # Plot the observation on a secondary y-axis for context (storm timing)
        obs_dates, obs = load_qkrig_obs(cat)
        if obs_dates is not None:
            od = pd.to_datetime(obs_dates)
            mask = (od >= ZOOM_START) & (od <= ZOOM_END)
            ax2 = ax.twinx()
            ax2.plot(od[mask], obs[mask],
                     color="black", lw=1.0, alpha=0.35, zorder=1, label="Qkrig (obs, ref)")
            ax2.set_ylabel("Qkrig (mm/h)", fontsize=8, color="0.4")
            ax2.tick_params(axis="y", labelsize=7, colors="0.4")
            ax2.spines["right"].set_color("0.7")

    ax.set_title(cat, fontsize=10)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", labelsize=8)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.set_ylabel("Ensemble std-dev (mm/h)", fontsize=8)


def make_main_figure():
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
    for ax, cat in zip(axes, MAIN_CATS):
        plot_one_panel(ax, cat, show_obs=True)
    axes[-1].legend(fontsize=8, loc="upper right")
    fig.suptitle(
        "Ensemble spread (std-dev across 20 members) by perturbation source\n"
        "Hurricane Helene peak (Sep 24-28, 2024) | DA off | "
        "worst / median / best Run 3 catchment\n"
        "Grey line on right axis = Qkrig observation (for storm-timing context only)",
        fontsize=11, y=1.06,
    )
    plt.tight_layout()
    plt.savefig(MAIN_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {MAIN_PNG}")


def make_appendix_figure():
    fig, axes = plt.subplots(3, 7, figsize=(28, 12), sharex=True)
    axes = axes.flatten()
    for i, cat in enumerate(ALL_CATS):
        plot_one_panel(axes[i], cat, show_obs=True)
    axes[0].legend(fontsize=7, loc="upper right")
    fig.suptitle(
        "Ensemble spread (hourly std-dev across 20 members) by perturbation source - "
        "Hurricane Helene peak (Sep 24-28, 2024) - all 21 sub-catchments\n"
        "Red = init state | Blue = forcing | Green = process noise | "
        "Grey on right axis = Qkrig obs (storm-timing context only)",
        fontsize=13, y=1.00,
    )
    plt.tight_layout()
    plt.savefig(APPENDIX_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {APPENDIX_PNG}")


def main():
    make_main_figure()
    make_appendix_figure()


if __name__ == "__main__":
    main()
