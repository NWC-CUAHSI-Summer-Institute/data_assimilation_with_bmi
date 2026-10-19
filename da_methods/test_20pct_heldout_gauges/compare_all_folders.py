"""
compare_all_folders.py

Four-way comparison of DA analysis trajectories routed to gauge 03463300,
one curve per R-formula experiment:
    Folder 1 — Variance scaled using Vrugt  R(t)=(0.10*y)^2 + 0.001*sigma2_krig
    Folder 2 — Fixed R=0.07                 R=0.07 (constant)
    Folder 3 — Dynamic Vrugt seeded         same Vrugt formula, spliced obs, seed=42
    Folder 4 — Dynamic variance direct      R(t)=sigma2_krig directly

Reads routed_Q_test.csv (columns: date, Q_routed_m3s, Q_usgs_m3s) from each folder.
Folder 2 deterministic route not available — shown as absent with a note.

Outputs:
    <out-dir>/compare_all_folders_full.png
    <out-dir>/compare_all_folders_helene.png
    <out-dir>/compare_all_folders_kge_table.csv

Usage:
    python3 compare_all_folders.py --out-dir /mnt/disk2/suma_helen_poster/da_results/comparison_plots
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

HELENE_START = pd.Timestamp("2024-09-24")
HELENE_END   = pd.Timestamp("2024-09-29 23:00:00")

BASE = "/mnt/disk2/suma_helen_poster/da_results"

FOLDERS = [
    {
        "label": "F1: Variance scaled Vrugt\nR(t)=(0.10·y)²+0.001·σ²_krig",
        "short": "F1 Vrugt",
        "color": "#1f77b4",
        "lw": 2.0,
        "csv": f"{BASE}/vrugt_dynamic_routed/routed_Q_test.csv",
    },
    {
        "label": "F2: Fixed R=0.07",
        "short": "F2 R=0.07",
        "color": "#ff7f0e",
        "lw": 1.6,
        "csv": f"{BASE}/fixed_r007_routed/routed_Q_test.csv",
    },
    {
        "label": "F3: Dynamic Vrugt seeded\n(spliced obs, seed=42)",
        "short": "F3 Dyn Vrugt",
        "color": "#2ca02c",
        "lw": 1.6,
        "csv": f"{BASE}/dynamic_vrugt_seeded_routed/routed_Q_test.csv",
    },
    {
        "label": "F4: Dynamic variance direct\nR(t)=σ²_krig",
        "short": "F4 Direct σ²",
        "color": "#d62728",
        "lw": 1.4,
        "csv": f"{BASE}/dynamic_novrugt_seeded_routed/routed_Q_test.csv",
    },
    {
        "label": "F5: Re-kriged variance\nR(t)=σ²_krig (re-kriged network)",
        "short": "F5 Rekrig σ²",
        "color": "#9467bd",
        "lw": 2.2,
        "csv": f"{BASE}/folder5_rekrig_variance_direct/routed_Q_test.csv",
    },
]


def kge(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.std(o) == 0:
        return np.nan
    r = np.corrcoef(o, s)[0, 1]
    return 1.0 - np.sqrt((r-1)**2 + (np.std(s)/np.std(o)-1)**2
                         + (np.mean(s)/np.mean(o)-1)**2)


def nse(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    o, s = obs[mask], sim[mask]
    if len(o) < 2:
        return np.nan
    denom = np.sum((o - o.mean())**2)
    return 1.0 - np.sum((o - s)**2) / denom if denom > 0 else np.nan


def load_folder(cfg):
    p = cfg["csv"]
    if not os.path.exists(p):
        print(f"  MISSING: {p}")
        return None
    df = pd.read_csv(p, parse_dates=["date"]).set_index("date").sort_index()
    return df


def plot_comparison(dfs, obs, dates, helene_mask, out_path, helene_only=False):
    fig, ax = plt.subplots(figsize=(15, 6))

    if helene_only:
        plot_dates = dates[helene_mask]
        obs_plot   = obs[helene_mask]
    else:
        plot_dates = dates
        obs_plot   = obs

    ax.plot(plot_dates, obs_plot,
            color="black", lw=2.2, zorder=6, label="USGS obs")

    for cfg, df in zip(FOLDERS, dfs):
        if df is None:
            continue
        sim = df["Q_routed_m3s"].reindex(dates).values
        if helene_only:
            sim_plot = sim[helene_mask]
        else:
            sim_plot = sim

        kg = kge(obs[helene_mask] if helene_only else obs,
                 sim[helene_mask] if helene_only else sim)
        ns = nse(obs[helene_mask] if helene_only else obs,
                 sim[helene_mask] if helene_only else sim)

        ax.plot(plot_dates, sim_plot,
                color=cfg["color"], lw=cfg["lw"], zorder=4,
                label=f"{cfg['short']}  KGE={kg:+.3f}  NSE={ns:+.3f}")

    if not helene_only:
        ax.axvspan(HELENE_START, HELENE_END, color="gold", alpha=0.12, zorder=1)

    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    title = ("Helene window — " if helene_only else "Full test period — ")
    ax.set_title(title + "Four-folder R-formula comparison at gauge 03463300", fontsize=12)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.9)
    ax.grid(True, alpha=0.22)

    if helene_only:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        peak_usgs = np.nanmax(obs[helene_mask])
        ax.axhline(peak_usgs, color="black", lw=0.7, linestyle=":", alpha=0.5)
        ax.text(HELENE_END - pd.Timedelta(hours=6), peak_usgs * 1.02,
                f"USGS peak {peak_usgs:.0f} m³/s",
                fontsize=8.5, ha="right", color="black")
    else:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=f"{BASE}/comparison_plots")
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading routed CSVs...")
    dfs = [load_folder(cfg) for cfg in FOLDERS]

    # Use first available df for date index and USGS obs
    ref = next(df for df in dfs if df is not None)
    dates  = ref.index
    obs    = ref["Q_usgs_m3s"].values if "Q_usgs_m3s" in ref.columns else np.full(len(dates), np.nan)
    helene = ((dates >= HELENE_START) & (dates <= HELENE_END))

    # ── Plots ──────────────────────────────────────────────────────────────
    plot_comparison(dfs, obs, dates, helene,
                    os.path.join(args.out_dir, "compare_all_folders_full.png"),
                    helene_only=False)
    plot_comparison(dfs, obs, dates, helene,
                    os.path.join(args.out_dir, "compare_all_folders_helene.png"),
                    helene_only=True)

    # ── KGE table ──────────────────────────────────────────────────────────
    rows = []
    for cfg, df in zip(FOLDERS, dfs):
        if df is None:
            rows.append({"folder": cfg["short"], "full_kge": np.nan,
                         "full_nse": np.nan, "helene_kge": np.nan,
                         "helene_nse": np.nan, "helene_peak_m3s": np.nan})
            continue
        sim = df["Q_routed_m3s"].reindex(dates).values
        rows.append({
            "folder":          cfg["short"],
            "full_kge":        round(kge(obs, sim), 3),
            "full_nse":        round(nse(obs, sim), 3),
            "helene_kge":      round(kge(obs[helene], sim[helene]), 3),
            "helene_nse":      round(nse(obs[helene], sim[helene]), 3),
            "helene_peak_m3s": round(np.nanmax(sim[helene]), 1),
        })

    tbl = pd.DataFrame(rows)
    tbl_path = os.path.join(args.out_dir, "compare_all_folders_kge_table.csv")
    tbl.to_csv(tbl_path, index=False)
    print(f"Saved: {tbl_path}")
    print("\n" + tbl.to_string(index=False))
    print(f"\nUSGS Helene peak: {np.nanmax(obs[helene]):.1f} m³/s")


if __name__ == "__main__":
    main()
