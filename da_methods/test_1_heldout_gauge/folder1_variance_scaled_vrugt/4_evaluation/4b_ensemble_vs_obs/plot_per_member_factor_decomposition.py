"""
Per-member factor-decomposition plot for ONE catchment.

For each of the 20 ensemble members, show the member's actual forecast
(production: init + forcing + process + DA all active) alongside the member's
trajectory under each ISOLATED perturbation source (init only / forcing only /
process noise only). The viewer reads each panel as 'why did THIS member give
THIS forecast — which perturbation source pushed it where?'

Note: 'same member_i' across the four CSVs is a column-name correspondence
only; the underlying random draws are independent across the four runs. So
each panel is an illustrative comparison, not a rigorous matched-seed Shapley
decomposition. The story still reads correctly: each colored line shows what
ONE realization of that perturbation source produces, and the thick line
shows what the full production setup produces.

Inputs (per catchment):
  <prod-dir>/<cat-id>/<cat-id>_production_per_member.csv     (20 cols)
  <sen-dir>/<cat-id>/<cat-id>_sensitivity_init.csv           (20 cols)
  <sen-dir>/<cat-id>/<cat-id>_sensitivity_forcing.csv        (20 cols)
  <sen-dir>/<cat-id>/<cat-id>_sensitivity_process.csv        (20 cols)
  <da-dir>/<cat-id>/<cat-id>_test_results.csv                (for Qkrig obs)

Output:
  <prod-dir>/<cat-id>/<cat-id>_per_member_factor_decomp.png
"""
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ---------- Configuration ----------
CAT = "cat-1016300"           # change here to do a different catchment

PROD_DIR = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
SEN_DIR  = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
DA_DIR   = "/mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt"
OBS_DIR  = "/home/svyas/catchment_ts_no_03463300_gapfilled"
OUT_PNG  = os.path.join(PROD_DIR, CAT, f"{CAT}_per_member_factor_decomp.png")

ZOOM_START = pd.Timestamp("2024-09-24")
ZOOM_END   = pd.Timestamp("2024-09-28")

PROD_COLOR    = "tab:purple"
INIT_COLOR    = "tab:red"
FORCING_COLOR = "tab:blue"
PROC_COLOR    = "tab:green"
OBS_COLOR     = "black"


def load_member_csv(path):
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, parse_dates=["date"])
    member_cols = sorted([c for c in df.columns if c.startswith("member_")])
    return df["date"].values, df[member_cols]


def load_obs():
    obs_p = os.path.join(OBS_DIR, f"{CAT}.csv")
    if os.path.exists(obs_p):
        df = pd.read_csv(obs_p)
        time_col = 'datetime' if 'datetime' in df.columns else 'date'
        df[time_col] = pd.to_datetime(df[time_col])
        return df[time_col].values, df["qkrig"].values
    p = os.path.join(DA_DIR, CAT, f"{CAT}_test_results.csv")
    if not os.path.exists(p):
        return None, None
    df = pd.read_csv(p, parse_dates=["date"])
    return df["date"].values, df["obs_mm_h"].values


def main():
    global CAT, OUT_PNG
    parser = argparse.ArgumentParser()
    parser.add_argument('--cat-id', default=CAT)
    args = parser.parse_args()
    CAT = args.cat_id
    OUT_PNG = os.path.join(PROD_DIR, CAT, f"{CAT}_per_member_factor_decomp.png")

    prod_dates, prod = load_member_csv(os.path.join(PROD_DIR, CAT, f"{CAT}_production_per_member.csv"))
    init_dates, init = load_member_csv(os.path.join(SEN_DIR,  CAT, f"{CAT}_sensitivity_init.csv"))
    forc_dates, forc = load_member_csv(os.path.join(SEN_DIR,  CAT, f"{CAT}_sensitivity_forcing.csv"))
    proc_dates, proc = load_member_csv(os.path.join(SEN_DIR,  CAT, f"{CAT}_sensitivity_process.csv"))
    obs_dates, obs   = load_obs()

    if prod is None:
        raise FileNotFoundError(
            f"Missing {PROD_DIR}/{CAT}/{CAT}_production_per_member.csv — "
            "run run_production_per_member.py first.")
    if init is None or forc is None or proc is None:
        raise FileNotFoundError(
            f"Missing one of the sensitivity CSVs under {SEN_DIR}/{CAT}/. "
            "Run run_perturbation_sensitivity.py for all three sources first.")

    def to_mask(dates):
        d = pd.to_datetime(dates)
        return d, (d >= ZOOM_START) & (d <= ZOOM_END)

    pd_dates, pd_mask = to_mask(prod_dates)
    id_dates, id_mask = to_mask(init_dates)
    fd_dates, fd_mask = to_mask(forc_dates)
    cd_dates, cd_mask = to_mask(proc_dates)
    if obs_dates is not None:
        od, om = to_mask(obs_dates)
    else:
        od, om = None, None

    N = prod.shape[1]
    # Lay out 4 cols x 5 rows = 20 panels
    n_cols = 4
    n_rows = (N + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(22, 4 * n_rows), sharex=True)
    axes = axes.flatten()

    for i in range(N):
        ax = axes[i]
        col = f"member_{i:02d}"

        # Each member's INIT-only trajectory (thin red)
        if col in init.columns:
            ax.plot(id_dates[id_mask], init[col].values[id_mask],
                    color=INIT_COLOR,   lw=0.9, alpha=0.9, label="Init only")
        # FORCING-only (thin blue)
        if col in forc.columns:
            ax.plot(fd_dates[fd_mask], forc[col].values[fd_mask],
                    color=FORCING_COLOR,lw=0.9, alpha=0.9, label="Forcing only")
        # PROCESS-only (thin green)
        if col in proc.columns:
            ax.plot(cd_dates[cd_mask], proc[col].values[cd_mask],
                    color=PROC_COLOR,   lw=0.9, alpha=0.9, label="Process noise only")
        # PRODUCTION (thick purple) — actual forecast with DA on, all perturbations
        if col in prod.columns:
            ax.plot(pd_dates[pd_mask], prod[col].values[pd_mask],
                    color=PROD_COLOR, lw=2.2, alpha=0.95, label="Production (DA on)")
        # Qkrig observation (thick black dashed)
        if od is not None:
            ax.plot(od[om], obs[om],
                    color=OBS_COLOR, lw=1.4, linestyle="--", alpha=0.9, label="Qkrig (obs)")

        ax.set_title(f"{col}", fontsize=10)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")

    # Hide any extra blank axes
    for j in range(N, len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        f"Per-member factor decomposition - {CAT} - Hurricane Helene peak (Sep 24-28, 2024)\n"
        "Purple thick = production member (DA on, all perturbations) | "
        "Red = init only | Blue = forcing only | Green = process noise only | "
        "Black dashed = Qkrig obs",
        fontsize=12, y=1.00,
    )
    fig.text(0.005, 0.5, "Discharge (mm/h)", va="center", rotation="vertical", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
