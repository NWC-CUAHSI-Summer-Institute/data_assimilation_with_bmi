"""
Per-member input/output diagnostic for ONE catchment.

For each ensemble member, show a single combined panel (hydrograph-style) that
traces the inputs the member actually received and the output it produced:

    Each member panel uses an internal 2-row stack:
        TOP sub-panel:   perturbed precip (bars, left axis)
                         + perturbed PET (line, right twinx axis)
        BOTTOM sub-panel: simulated streamflow Q (purple)
                          + Qkrig observation (black dashed)
        Shared x-axis between the two sub-panels.

    Initial states for that member are shown in the title line.

Layout: 5 rows x 4 columns = 20 member panels.

Inputs (all from run_production_per_member.py):
    <PROD-DIR>/<CAT>/<CAT>_production_per_member.csv               (Q per member)
    <PROD-DIR>/<CAT>/<CAT>_production_per_member_precip.csv        (precip per member)
    <PROD-DIR>/<CAT>/<CAT>_production_per_member_pet.csv           (PET per member)
    <PROD-DIR>/<CAT>/<CAT>_production_per_member_initial_states.json
    <DA-DIR>/<CAT>/<CAT>_test_results.csv                          (Qkrig obs)

Output:
    <PROD-DIR>/<CAT>/<CAT>_per_member_io_diagnostic.png
"""
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec

# ----- Configuration -----
CAT = "cat-1016300"

PROD_DIR = "/mnt/disk2/suma_helen_poster/da_results/v2_production_per_member"
DA_DIR   = "/mnt/disk2/suma_helen_poster/da_results/v2_true_enkf_vrugt"
OUT_PNG  = os.path.join(PROD_DIR, CAT, f"{CAT}_per_member_io_diagnostic.png")

# Plot window (Helene 4-day zoom)
ZOOM_START = pd.Timestamp("2024-09-24")
ZOOM_END   = pd.Timestamp("2024-09-28")

PRECIP_COLOR = "tab:blue"
PET_COLOR    = "tab:orange"
Q_COLOR      = "tab:purple"
OBS_COLOR    = "black"

N_ROWS = 5      # grid rows of member panels
N_COLS = 4      # grid cols of member panels


def load_member_csv(path):
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    return df["date"].values, df[member_cols]


def main():
    q_dates,      q_df      = load_member_csv(os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member.csv"))
    precip_dates, precip_df = load_member_csv(os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member_precip.csv"))
    pet_dates,    pet_df    = load_member_csv(os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member_pet.csv"))

    init_states_path = os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member_initial_states.json")
    if os.path.exists(init_states_path):
        with open(init_states_path) as f:
            init_states = json.load(f)["members"]
    else:
        init_states = {}

    obs_path = os.path.join(DA_DIR, CAT, f"{CAT}_test_results.csv")
    if os.path.exists(obs_path):
        obs_df = pd.read_csv(obs_path, parse_dates=["date"])
        obs_dates = pd.to_datetime(obs_df["date"].values)
        obs_vals  = obs_df["obs_mm_h"].values
    else:
        obs_dates = None
        obs_vals  = None

    if q_df is None or precip_df is None or pet_df is None:
        raise FileNotFoundError(
            "Missing one of the per-member CSVs. Run run_production_per_member.py "
            "with the updated script that saves precip + PET + initial states.")

    member_cols = list(q_df.columns)
    N = len(member_cols)
    assert N <= N_ROWS * N_COLS, f"Grid too small for {N} members"

    qd = pd.to_datetime(q_dates)
    pd_ = pd.to_datetime(precip_dates)
    ed = pd.to_datetime(pet_dates)
    qmask = (qd >= ZOOM_START) & (qd <= ZOOM_END)
    pmask = (pd_ >= ZOOM_START) & (pd_ <= ZOOM_END)
    emask = (ed >= ZOOM_START) & (ed <= ZOOM_END)
    if obs_dates is not None:
        omask = (obs_dates >= ZOOM_START) & (obs_dates <= ZOOM_END)

    # Compute global y-axis ranges so all member panels are directly comparable
    p_ymax = float(np.nanmax(precip_df.values[pmask, :])) * 1.10
    e_ymax = float(np.nanmax(pet_df.values[emask, :])) * 1.10
    q_data_max = float(np.nanmax(q_df.values[qmask, :]))
    if obs_dates is not None:
        q_data_max = max(q_data_max, float(np.nanmax(obs_vals[omask])))
    q_ymax = q_data_max * 1.10

    fig = plt.figure(figsize=(22, 19))
    outer = gridspec.GridSpec(N_ROWS, N_COLS, hspace=0.75, wspace=0.30)

    for i in range(N):
        col = member_cols[i]
        row, c = divmod(i, N_COLS)
        # Each member cell is a 2-row inner grid (precip+PET on top, Q on bottom)
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row, c], hspace=0.20, height_ratios=[1.0, 1.6])

        # --- Top sub-panel: precip bars + PET line on twinx ---
        ax_top = fig.add_subplot(inner[0])
        ax_top.bar(pd_[pmask], precip_df[col].values[pmask],
                   width=0.04, color=PRECIP_COLOR, alpha=0.85, label="P")
        ax_top.set_ylim(0, p_ymax)
        ax_top.set_ylabel("P (mm/h)", fontsize=8, color=PRECIP_COLOR)
        ax_top.tick_params(axis="y", labelsize=7, labelcolor=PRECIP_COLOR)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_top.grid(True, alpha=0.15)
        # PET on twinx (different scale)
        ax_pet = ax_top.twinx()
        ax_pet.plot(ed[emask], pet_df[col].values[emask],
                    color=PET_COLOR, lw=1.0, label="PET")
        ax_pet.set_ylim(0, e_ymax)
        ax_pet.set_ylabel("PET (mm/h)", fontsize=8, color=PET_COLOR)
        ax_pet.tick_params(axis="y", labelsize=7, labelcolor=PET_COLOR)
        # Title: member name on line 1, initial states on line 2 (avoids overlap with neighbors)
        s = init_states.get(col, {})
        title = (f"{col}\n"
                 f"soil={s.get('soil_m', float('nan')):.3f} m  |  "
                 f"GW={s.get('gw_m', float('nan')):.4f} m  |  "
                 f"Nash[0]={s.get('nash0_m', float('nan')):.1e}  "
                 f"Nash[1]={s.get('nash1_m', float('nan')):.1e}")
        ax_top.set_title(title, fontsize=8, loc="left")

        # --- Bottom sub-panel: Q + obs (shares x-axis with top) ---
        ax_bot = fig.add_subplot(inner[1], sharex=ax_top)
        ax_bot.plot(qd[qmask], q_df[col].values[qmask],
                    color=Q_COLOR, lw=1.6, label=col)
        if obs_dates is not None:
            ax_bot.plot(obs_dates[omask], obs_vals[omask],
                        color=OBS_COLOR, lw=1.2, linestyle="--", label="Qkrig (obs)")
        ax_bot.set_ylim(0, q_ymax)
        ax_bot.set_ylabel("Q (mm/h)", fontsize=8)
        ax_bot.tick_params(labelsize=7)
        ax_bot.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax_bot.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax_bot.grid(True, alpha=0.15)
        if i == 0:
            ax_bot.legend(fontsize=7, loc="upper left", ncol=2,
                          frameon=True, framealpha=0.85)

    fig.suptitle(
        f"Per-member input/output diagnostic - {CAT} - "
        f"Hurricane Helene peak (Sep 24-28, 2024)\n"
        "Each member panel: TOP = perturbed precip (blue bars) + perturbed PET (orange line)   |   "
        "BOTTOM = simulated Q (purple) vs Qkrig obs (black dashed).   "
        "Title shows initial states at t=0.",
        fontsize=12, y=0.995,
    )
    plt.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
