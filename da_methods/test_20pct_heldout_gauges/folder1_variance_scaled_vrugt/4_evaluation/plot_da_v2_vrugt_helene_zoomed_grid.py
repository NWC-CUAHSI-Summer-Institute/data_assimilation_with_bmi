"""
Hurricane Helene PEAK zoom: Run 3 (no DA) vs DA v2 (true EnKF + Vrugt R) vs Qkrig obs.
3 x 7 grid of all 21 catchments, zoomed to the 4-day peak window Sep 24 - Sep 28, 2024.

Same data as plot_da_v2_vrugt_helene_grid.py, just a tighter x-axis to see
the peak detail.

Produces: helene_da_v2_vrugt_vs_run3_grid_zoomed.png

Run inside the same conda env that has pandas + matplotlib (e.g. troute).
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

CATS = [
    "cat-1016279", "cat-1016280", "cat-1016281", "cat-1016282", "cat-1016283",
    "cat-1016300", "cat-1016301", "cat-1016302", "cat-1016303", "cat-1016304",
    "cat-1016305", "cat-1016306", "cat-1016307", "cat-1016308", "cat-1016309",
    "cat-1016310", "cat-1016311", "cat-1016312", "cat-1016313", "cat-1016314",
    "cat-1016315",
]

RUN3_DIR = "/mnt/disk2/suma_helen_poster/catchment_results_range100_run3"
DA_DIR   = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"
OUT_PNG  = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt/helene_da_v2_vrugt_vs_run3_grid_zoomed.png"
DA_LABEL = "DA v2 (true EnKF + Vrugt R, N=20)"

# Tighter 4-day peak window
ZOOM_START = pd.Timestamp("2024-09-24")
ZOOM_END   = pd.Timestamp("2024-09-28")


def kge(obs, sim):
    """Kling-Gupta Efficiency on full overlap of obs and sim."""
    m = ~(pd.isna(obs) | pd.isna(sim))
    o, s = obs[m], sim[m]
    if len(o) < 2:
        return float('nan')
    denom = (((o - o.mean()) ** 2).sum() * ((s - s.mean()) ** 2).sum()) ** 0.5
    if denom == 0:
        return float('nan')
    r = ((o - o.mean()) * (s - s.mean())).sum() / denom
    alpha = s.std() / o.std()
    beta = s.mean() / o.mean() if o.mean() != 0 else float('nan')
    return 1 - ((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2) ** 0.5


def main():
    fig, axes = plt.subplots(3, 7, figsize=(28, 12), sharex=True)
    axes = axes.flatten()

    for i, cat in enumerate(CATS):
        ax = axes[i]
        run3_csv = os.path.join(RUN3_DIR, cat, f"{cat}_test_results.csv")
        da_csv   = os.path.join(DA_DIR,   cat, f"{cat}_test_results.csv")

        if not (os.path.exists(run3_csv) and os.path.exists(da_csv)):
            ax.set_title(f"{cat} (missing)")
            continue

        df_run3 = pd.read_csv(run3_csv, parse_dates=["date"])
        df_da   = pd.read_csv(da_csv,   parse_dates=["date"])

        mask3 = (df_run3["date"] >= ZOOM_START) & (df_run3["date"] <= ZOOM_END)
        maskd = (df_da["date"]   >= ZOOM_START) & (df_da["date"]   <= ZOOM_END)

        # KGE on full test period (matches log-line value); displayed in legend
        kge_run3 = kge(df_run3["obs_mm_h"], df_run3["sim_mm_h"])
        kge_da   = kge(df_da["obs_mm_h"],   df_da["sim_mm_h"])

        ax.plot(df_da.loc[maskd,  "date"], df_da.loc[maskd,  "obs_mm_h"],
                color="black", lw=1.8, label="Qkrig (obs)", zorder=1)
        ax.plot(df_run3.loc[mask3, "date"], df_run3.loc[mask3, "sim_mm_h"],
                color="steelblue", lw=4.0, alpha=0.45, label=f"Run 3 KGE={kge_run3:.2f}", zorder=2)
        ax.plot(df_da.loc[maskd,  "date"], df_da.loc[maskd,  "sim_mm_h"],
                color="tomato", lw=1.4, label=f"DA KGE={kge_da:.2f}", zorder=3)

        ax.set_title(f"{cat}", fontsize=10)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6, 12, 18]))
        ax.legend(fontsize=6, loc="upper right")

    fig.suptitle(
        f"Hurricane Helene peak zoom (Sep 24 - Sep 28, 2024) - {DA_LABEL} vs Run 3 baseline",
        fontsize=14, y=1.00,
    )
    fig.text(0.005, 0.5, "Discharge (mm/h)", va="center", rotation="vertical", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
