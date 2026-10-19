"""
Plot perturbation-source sensitivity spaghetti plots.

Reads the three per-source sensitivity CSVs per catchment and plots all 20
member streamflow trajectories as colored bundles overlaid on one panel,
zoomed to the Hurricane Helene peak window (Sep 24 - Sep 28, 2024).

Color legend:
    red   = initial state perturbation only
    blue  = forcing perturbation only (precip + PET)
    green = process noise on hydrologic states only

Produces:
    helene_sensitivity_spaghetti_main.png       (3-catchment subset for the main figure)
    helene_sensitivity_spaghetti_appendix.png   (full 3 x 7 grid of all 21 catchments)

Requires the sensitivity runs to have completed first (see
run_perturbation_sensitivity.py).
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

# ---------- Configuration ----------
ALL_CATS = [
    "cat-1016279", "cat-1016280", "cat-1016281", "cat-1016282", "cat-1016283",
    "cat-1016300", "cat-1016301", "cat-1016302", "cat-1016303", "cat-1016304",
    "cat-1016305", "cat-1016306", "cat-1016307", "cat-1016308", "cat-1016309",
    "cat-1016310", "cat-1016311", "cat-1016312", "cat-1016313", "cat-1016314",
    "cat-1016315",
]

# Main figure: 3 catchments spanning the Run 3 KGE distribution.
# cat-1016311 = worst Run 3 (0.50); cat-1016300 = median (0.74); cat-1016302 = best (0.83)
MAIN_CATS = ["cat-1016311", "cat-1016300", "cat-1016302"]

DA_DIR  = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"   # for Qkrig obs
SEN_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_sensitivity"        # sensitivity CSVs
MAIN_PNG     = "/mnt/disk2/suma_helen_poster/da_results/v2_sensitivity/helene_sensitivity_spaghetti_main.png"
APPENDIX_PNG = "/mnt/disk2/suma_helen_poster/da_results/v2_sensitivity/helene_sensitivity_spaghetti_appendix.png"

# Helene 4-day peak zoom
ZOOM_START = pd.Timestamp("2024-09-24")
ZOOM_END   = pd.Timestamp("2024-09-28")

SOURCE_COLORS = {
    "init":    "tab:red",
    "forcing": "tab:blue",
    "process": "tab:green",
}
SOURCE_LABELS = {
    "init":    "Initial state only",
    "forcing": "Forcing only (P, PET)",
    "process": "Process noise on states only",
}


def load_sensitivity_csv(cat, source):
    """Returns (dates, q_matrix) for the test period (member columns)."""
    path = os.path.join(SEN_DIR, cat, f"{cat}_sensitivity_{source}.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    member_cols = [c for c in df.columns if c.startswith("member_")]
    return df["date"].values, df[member_cols].values


def load_qkrig_obs(cat):
    """Loads the Qkrig observation timeseries from the production-run CSV."""
    path = os.path.join(DA_DIR, cat, f"{cat}_test_results.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    return df["date"].values, df["obs_mm_h"].values


def plot_one_panel(ax, cat):
    """Render the three colored bundles plus the obs line for one catchment."""
    # Each source: 20 thin colored lines
    for source, color in SOURCE_COLORS.items():
        dates, q = load_sensitivity_csv(cat, source)
        if dates is None:
            continue
        df_dates = pd.to_datetime(dates)
        mask = (df_dates >= ZOOM_START) & (df_dates <= ZOOM_END)
        if mask.sum() == 0:
            continue
        # Plot all 20 members as thin transparent lines, plus a thicker mean line
        for i in range(q.shape[1]):
            ax.plot(df_dates[mask], q[mask, i],
                    color=color, lw=0.5, alpha=0.35, zorder=1)
        ax.plot(df_dates[mask], q[mask].mean(axis=1),
                color=color, lw=1.6, alpha=0.95, zorder=2,
                label=SOURCE_LABELS[source])

    # Black observation overlay
    obs_dates, obs = load_qkrig_obs(cat)
    if obs_dates is not None:
        od = pd.to_datetime(obs_dates)
        mask = (od >= ZOOM_START) & (od <= ZOOM_END)
        ax.plot(od[mask], obs[mask],
                color="black", lw=1.4, label="Qkrig (obs)", zorder=3)

    ax.set_title(cat, fontsize=10)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", labelsize=8)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))


def make_main_figure():
    """3-catchment side-by-side panel for the main paper / poster figure."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
    for ax, cat in zip(axes, MAIN_CATS):
        plot_one_panel(ax, cat)
    # Single legend on the rightmost panel (avoid clutter on the others)
    axes[-1].legend(fontsize=8, loc="upper right")
    axes[0].set_ylabel("Discharge (mm/h)", fontsize=11)
    fig.suptitle(
        "Ensemble spread by perturbation source — Hurricane Helene peak (Sep 24-28, 2024)\n"
        "20 members per source. DA off. Worst / median / best Run 3 catchment.",
        fontsize=12, y=1.02,
    )
    plt.tight_layout()
    plt.savefig(MAIN_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {MAIN_PNG}")


def make_appendix_figure():
    """Full 3 x 7 grid of all 21 catchments."""
    fig, axes = plt.subplots(3, 7, figsize=(28, 12), sharex=True)
    axes = axes.flatten()
    for i, cat in enumerate(ALL_CATS):
        plot_one_panel(axes[i], cat)
    # Single legend on the first panel
    axes[0].legend(fontsize=7, loc="upper right")
    fig.text(0.005, 0.5, "Discharge (mm/h)", va="center", rotation="vertical", fontsize=12)
    fig.suptitle(
        "Ensemble spread by perturbation source - Hurricane Helene peak (Sep 24-28, 2024) - all 21 sub-catchments\n"
        "Red = initial state only | Blue = forcing only | Green = process noise only | Black = Qkrig obs",
        fontsize=14, y=1.00,
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
