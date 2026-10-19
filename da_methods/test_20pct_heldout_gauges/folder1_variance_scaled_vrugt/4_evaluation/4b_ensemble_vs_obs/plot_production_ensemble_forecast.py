"""
Single-panel ensemble forecast plot for one catchment, paper-style.

Plots all 20 production-per-member streamflow trajectories overlaid on a
single time-series panel, with the Hurricane Helene window highlighted by a
pink shaded vertical band and the Qkrig observation drawn on top in black.

Ensemble forecast figure style (50/100/200 km
variogram-range ensemble panels): every member as a thin colored line, storm
window shaded, observation as a thick dark series, optional log-scaled y-axis.

Inputs:
  <PROD-DIR>/<CAT>/<CAT>_production_per_member.csv     (date + 20 members)
  <DA-DIR>/<CAT>/<CAT>_test_results.csv                (Qkrig obs)

Outputs:
  <PROD-DIR>/<CAT>/<CAT>_production_ensemble_forecast_linear.png
  <PROD-DIR>/<CAT>/<CAT>_production_ensemble_forecast_log.png
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from matplotlib import cm

# ----- Configuration -----
CAT = "cat-1016300"

PROD_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_production_per_member"
DA_DIR   = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"

# Plot window: a few weeks around Helene so the storm is in context
PLOT_START = pd.Timestamp("2024-09-20")
PLOT_END   = pd.Timestamp("2024-10-05")

# Helene peak window highlighted in pink
HELENE_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_END   = pd.Timestamp("2024-09-28 00:00:00")

OUT_LINEAR = os.path.join(PROD_DIR, CAT, f"{CAT}_production_ensemble_forecast_linear.png")
OUT_LOG    = os.path.join(PROD_DIR, CAT, f"{CAT}_production_ensemble_forecast_log.png")


def kge_score(obs, sim):
    m = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[m], sim[m]
    if len(o) < 2:
        return float("nan")
    denom = float(np.sqrt(((o - o.mean()) ** 2).sum() * ((s - s.mean()) ** 2).sum()))
    if denom == 0:
        return float("nan")
    r = float(((o - o.mean()) * (s - s.mean())).sum() / denom)
    alpha = float(s.std() / o.std()) if o.std() != 0 else float("nan")
    beta = float(s.mean() / o.mean()) if o.mean() != 0 else float("nan")
    return 1.0 - float(np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def main():
    prod_path = os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member.csv")
    obs_path  = os.path.join(DA_DIR,   CAT, f"{CAT}_test_results.csv")

    if not os.path.exists(prod_path):
        raise FileNotFoundError(
            f"Missing {prod_path}. Run run_production_per_member.py for {CAT} first.")
    if not os.path.exists(obs_path):
        print(f"Note: {obs_path} not found; will plot without observation overlay.")

    prod_df = pd.read_csv(prod_path, parse_dates=["date"])
    member_cols = sorted([c for c in prod_df.columns if c.startswith("member_")])
    member_arr  = prod_df[member_cols].to_numpy(dtype=float)
    dates       = pd.to_datetime(prod_df["date"].values)

    if os.path.exists(obs_path):
        obs_df = pd.read_csv(obs_path, parse_dates=["date"])
        obs_dates = pd.to_datetime(obs_df["date"].values)
        obs_vals  = obs_df["obs_mm_h"].values
    else:
        obs_dates = None
        obs_vals  = None

    # Compute per-member peak (within plot window) for the legend label
    pw_mask = (dates >= PLOT_START) & (dates <= PLOT_END)
    peaks = member_arr[pw_mask, :].max(axis=0)

    # Optional: compute per-member KGE vs obs across the plot window
    member_kges = np.full(len(member_cols), np.nan)
    if obs_dates is not None:
        obs_series = pd.Series(obs_vals, index=obs_dates)
        for i, _ in enumerate(member_cols):
            member_series = pd.Series(member_arr[:, i], index=dates)
            joined = pd.concat([obs_series, member_series], axis=1, join="inner").dropna()
            if len(joined) >= 2:
                member_kges[i] = kge_score(joined.iloc[:, 0].values, joined.iloc[:, 1].values)

    def plot_panel(ax, log_y=False):
        colormap = cm.get_cmap("turbo", len(member_cols))
        # Thin colored lines per member
        for i, col in enumerate(member_cols):
            label = f"{col} | peak={peaks[i]:.1f} mm/h"
            ax.plot(dates[pw_mask], member_arr[pw_mask, i],
                    color=colormap(i), lw=0.8, alpha=0.85, label=label, zorder=2)

        # Helene peak shaded band
        ax.axvspan(HELENE_START, HELENE_END, color="salmon", alpha=0.18, zorder=1)
        ax.text(HELENE_START + (HELENE_END - HELENE_START) / 2,
                ax.get_ylim()[1] if not log_y else 1.0,
                "Helene peak", fontsize=8, color="salmon",
                ha="center", va="bottom", zorder=3)

        # Observation as thick black series
        if obs_dates is not None:
            obs_mask = (obs_dates >= PLOT_START) & (obs_dates <= PLOT_END)
            ax.plot(obs_dates[obs_mask], obs_vals[obs_mask],
                    color="black", lw=1.6, label="Qkrig (obs)", zorder=4)

        ax.set_xlabel("Date", fontsize=10)
        if log_y:
            ax.set_yscale("log")
            ax.set_ylabel(r"q (mm hr$^{-1}$) [log scale]", fontsize=10)
            ax.set_ylim(0.01, max(member_arr.max(), (obs_vals.max() if obs_vals is not None else 0)) * 1.2)
        else:
            ax.set_ylabel("Discharge (mm/h)", fontsize=10)

        ax.tick_params(axis="x", rotation=30, labelsize=9)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.grid(True, alpha=0.2)

    # ---- LINEAR Y figure ----
    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, log_y=False)
    ax.legend(fontsize=6.5, loc="upper left", ncol=2,
              frameon=True, framealpha=0.85, markerfirst=False)
    fig.suptitle(
        f"Production ensemble forecast (N=20) at {CAT} - "
        f"Sep 20 - Oct 5, 2024 (Hurricane Helene window)",
        fontsize=12, y=0.99,
    )
    plt.tight_layout()
    plt.savefig(OUT_LINEAR, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LINEAR}")

    # ---- LOG Y figure (matches paper style) ----
    fig, ax = plt.subplots(1, 1, figsize=(20, 7))
    plot_panel(ax, log_y=True)
    ax.legend(fontsize=6.5, loc="upper left", ncol=2,
              frameon=True, framealpha=0.85, markerfirst=False)
    fig.suptitle(
        f"Production ensemble forecast (N=20) at {CAT} - log-scale q - "
        f"Sep 20 - Oct 5, 2024 (Hurricane Helene window)",
        fontsize=12, y=0.99,
    )
    plt.tight_layout()
    plt.savefig(OUT_LOG, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_LOG}")


if __name__ == "__main__":
    main()
