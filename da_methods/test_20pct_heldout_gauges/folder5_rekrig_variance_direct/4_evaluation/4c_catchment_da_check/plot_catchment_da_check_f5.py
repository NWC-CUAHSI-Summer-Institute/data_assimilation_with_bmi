"""
plot_catchment_da_check_f5.py

Catchment-level sanity check for F5 (re-kriged variance, 20% holdout).

Reads <cat>_test_results.csv (DA ensemble mean output from Step 1) and plots:
  - sim_mm_h  : DA ensemble mean streamflow prediction
  - obs_mm_h  : Qkrig observation DA was trying to match
  - precip    : rainfall (inverted second axis)

Two panels per catchment:
  Panel 1 : full test period  (Oct 2023 - Oct 2024)
  Panel 2 : Helene zoom       (Sep 20 - Oct 5 2024)

NSE at catchment level (sim vs Qkrig) reported in title.
This figure checks DA in isolation — before routing — so any
remaining error at the gauge can be attributed to routing vs DA.

Run on server:
    python3 plot_catchment_da_check_f5.py
    python3 plot_catchment_da_check_f5.py --da-dir /path/to/da --out-dir /path/to/out
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_DA_DIR  = "/mnt/disk2/suma_helen_poster/da_results/folder5_rekrig_variance_direct"
DEFAULT_OUT_DIR = None   # falls back to DEFAULT_DA_DIR

CAT_ID = "cat-1016300"  # focal catchment draining to gauge 03463300

HELENE_ZOOM_START = pd.Timestamp("2024-09-20")
HELENE_ZOOM_END   = pd.Timestamp("2024-10-05")
HELENE_PEAK_START = pd.Timestamp("2024-09-24")
HELENE_PEAK_END   = pd.Timestamp("2024-09-30")

SIM_COLOR = "tomato"
OBS_COLOR = "black"


def nse(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2:
        return float("nan")
    denom = ((o - o.mean()) ** 2).sum()
    return float(1.0 - ((o - s) ** 2).sum() / denom) if denom > 0 else float("nan")


def plot_catchment(cat_id, df, out_dir):
    dates   = pd.to_datetime(df["date"])
    sim     = df["sim_mm_h"].values.astype(float)
    obs     = df["obs_mm_h"].values.astype(float)
    precip  = df["precip_mm_h"].values.astype(float)

    nse_full = nse(obs, sim)

    helene_mask = (dates >= HELENE_ZOOM_START) & (dates <= HELENE_ZOOM_END)
    nse_helene  = nse(obs[helene_mask], sim[helene_mask])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 9))
    fig.suptitle(
        f"{cat_id}  |  F5 catchment-level DA check (re-kriged σ², 20% holdout)\n"
        f"NSE = {nse_full:.3f} (full year)   NSE = {nse_helene:.3f} (Helene zoom)",
        fontsize=12,
    )

    # ── Panel 1: full period ──────────────────────────────────────────
    ax1.plot(dates, sim, color=SIM_COLOR, lw=1.2, label="DA ensemble mean  [mm/h]")
    ax1.plot(dates, obs, color=OBS_COLOR, lw=0.8, alpha=0.7, label="Qkrig obs  [mm/h]")
    ax1.set_ylabel("Streamflow (mm/h)", fontsize=10)
    ax1.set_title("Full test period: Oct 2023 – Oct 2024", fontsize=10)
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(True, alpha=0.2)
    ax1_twin = ax1.twinx()
    ax1_twin.bar(dates, precip, color="steelblue", alpha=0.35, width=0.04)
    ax1_twin.set_ylim([precip.max() * 4, 0])
    ax1_twin.set_ylabel("Precip (mm/h)", fontsize=9, color="steelblue")

    # ── Panel 2: Helene zoom ──────────────────────────────────────────
    d_h   = dates[helene_mask]
    sim_h = sim[helene_mask]
    obs_h = obs[helene_mask]
    pre_h = precip[helene_mask]

    ax2.axvspan(HELENE_PEAK_START, HELENE_PEAK_END, color="#ffcccc", alpha=0.35, zorder=0)
    ax2.plot(d_h, sim_h, color=SIM_COLOR, lw=2.0, label="DA ensemble mean  [mm/h]")
    ax2.plot(d_h, obs_h, color=OBS_COLOR, lw=1.5, label="Qkrig obs  [mm/h]")
    ax2.set_ylabel("Streamflow (mm/h)", fontsize=10)
    ax2.set_title("Helene zoom: Sep 20 – Oct 5 2024  (pink = Sep 24–30)", fontsize=10)
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.2)
    ax2_twin = ax2.twinx()
    ax2_twin.bar(d_h, pre_h, color="steelblue", alpha=0.40, width=0.04)
    ax2_twin.set_ylim([pre_h.max() * 4, 0])
    ax2_twin.set_ylabel("Precip (mm/h)", fontsize=9, color="steelblue")

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{cat_id}_catchment_da_check.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_path}  (NSE={nse_full:.3f})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--da-dir",  default=DEFAULT_DA_DIR)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    out_dir = args.out_dir or args.da_dir
    os.makedirs(out_dir, exist_ok=True)

    csv_path = os.path.join(args.da_dir, CAT_ID, f"{CAT_ID}_test_results.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(csv_path)
    df = pd.read_csv(csv_path)
    plot_catchment(CAT_ID, df, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
