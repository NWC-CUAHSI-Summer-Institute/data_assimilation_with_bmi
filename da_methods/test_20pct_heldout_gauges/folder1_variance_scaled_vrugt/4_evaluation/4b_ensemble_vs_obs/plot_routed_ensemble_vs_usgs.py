#!/usr/bin/env python3
"""
plot_routed_ensemble_vs_usgs.py
Plot T-route routed 600-member ensemble envelope vs USGS at gauge 03463300.

Input: routed_crossed_ensemble.parquet from run_route_crossed_ensemble.py
       Columns: issue_time, lead_hour, member_0000..member_0599  (Q in m3/s)

Usage:
    python3 plot_routed_ensemble_vs_usgs.py \
        --routed-pq  /mnt/disk2/.../folder1_vrugt_routed/routed_crossed_ensemble.parquet \
        --usgs-csv   /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv \
        --label      "F1 Vrugt — 1 gauge holdout" \
        --out-dir    ~/plots_f1_1gauge
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WATERSHED_AREA_KM2 = 113.18          # cat-1016300 drainage area
MM_H_TO_M3_S       = WATERSHED_AREA_KM2 * 1000.0 / 3600.0

PLOT_START = pd.Timestamp("2024-09-24 00:00:00")
PLOT_END   = pd.Timestamp("2024-09-28 18:00:00")
PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
PEAK_END   = pd.Timestamp("2024-09-28 06:00:00")
HELENE_PEAK = pd.Timestamp("2024-09-27 12:00:00")

COLOR_ENS = "#1f77b4"
COLOR_OBS = "black"


def build_envelope(pq_path, t_start, t_end):
    df = pd.read_parquet(pq_path)
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    member_cols = [c for c in df.columns if c.startswith("member_")]

    groups = defaultdict(list)
    for _, row in df.iterrows():
        vt = (row["issue_time"] + pd.Timedelta(hours=int(row["lead_hour"]))).floor("1h")
        if t_start <= vt <= t_end:
            groups[vt].extend(row[member_cols].values.astype(float).tolist())

    rows = []
    for vt in sorted(groups):
        v = np.array(groups[vt])
        v = v[~np.isnan(v)]
        if len(v) == 0:
            continue
        rows.append({"valid_time": vt,
                     "p05": np.percentile(v, 5),  "p25": np.percentile(v, 25),
                     "p50": np.percentile(v, 50), "p75": np.percentile(v, 75),
                     "p95": np.percentile(v, 95)})
    return pd.DataFrame(rows).set_index("valid_time")


def load_usgs(path):
    df = pd.read_csv(path)
    print(f"  USGS CSV columns: {list(df.columns)}")
    date_col = next(c for c in df.columns
                    if c.lower() in ("datetime", "date", "time", "timestamp"))
    q_col = next(c for c in df.columns
                 if "q" in c.lower() or "flow" in c.lower() or "discharge" in c.lower())
    print(f"  Using date_col={date_col!r}, q_col={q_col!r}")
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if "mm" in q_col.lower():
        series = series * MM_H_TO_M3_S
    return series.sort_index()


def compute_kge(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim))
    o, s = obs[mask], sim[mask]
    if len(o) < 2 or np.std(o) == 0:
        return np.nan
    r = np.corrcoef(o, s)[0, 1]
    return 1 - np.sqrt((r-1)**2 + (np.std(s)/np.std(o)-1)**2 + (np.mean(s)/np.mean(o)-1)**2)


def _shade(ax, env):
    ax.fill_between(env.index, env["p05"], env["p95"],
                    color=COLOR_ENS, alpha=0.18, label="5th–95th percentile")
    ax.fill_between(env.index, env["p25"], env["p75"],
                    color=COLOR_ENS, alpha=0.32, label="25th–75th percentile")
    ax.plot(env.index, env["p50"], color=COLOR_ENS, lw=1.5, label="Ensemble median (600 members)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--routed-pq", required=True)
    parser.add_argument("--usgs-csv",  required=True)
    parser.add_argument("--label",     default="F1 Vrugt — 1 gauge holdout")
    parser.add_argument("--out-dir",   required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Building full-window envelope...")
    env_full = build_envelope(args.routed_pq, PLOT_START, PLOT_END)
    print("Building peak-window envelope...")
    env_peak = build_envelope(args.routed_pq, PEAK_START, PEAK_END)

    print("Loading USGS obs...")
    try:
        obs_all = load_usgs(args.usgs_csv)
        obs_full = obs_all.loc[PLOT_START:PLOT_END]
        obs_peak = obs_all.loc[PEAK_START:PEAK_END]
        print(f"  Loaded {len(obs_all)} USGS obs rows; "
              f"full window: {len(obs_full)}, peak window: {len(obs_peak)}")
    except Exception as e:
        print(f"  Warning: USGS load failed ({e}) — plotting without obs")
        obs_full = obs_peak = None

    kge = np.nan
    if obs_full is not None and len(env_full) > 0:
        med = env_full["p50"]
        obs_hourly = obs_full.copy()
        obs_hourly.index = obs_hourly.index.floor("1h")
        obs_hourly = obs_hourly[~obs_hourly.index.duplicated(keep="first")]
        aligned = obs_hourly.reindex(med.index)
        print(f"  KGE alignment: {aligned.notna().sum()} matched out of {len(med)} timesteps")
        kge = compute_kge(aligned.values, med.values)

    # ── Figure 1: Full Helene window ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 4.5))
    _shade(ax, env_full)
    if obs_full is not None:
        ax.plot(obs_full.index, obs_full.values, color=COLOR_OBS, lw=1.6,
                label="USGS obs (gauge 03463300)")
    ax.axvspan(PEAK_START, PEAK_END, color="salmon", alpha=0.12, zorder=0, label="Helene peak window")
    ax.set_ylabel("Discharge (m³/s)")
    ax.set_xlabel("Date")
    ax.set_title(f"Routed 600-member ensemble vs USGS | {args.label}\n"
                 f"Ensemble median KGE = {kge:.3f}  |  Sep 24–28 2024")
    ax.legend(fontsize=9, loc="upper left")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    plt.tight_layout()
    out = out_dir / "routed_ensemble_vs_usgs_full.png"
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.close(fig)

    # ── Figure 2: Peak zoom ───────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4.5))
    _shade(ax, env_peak)
    if obs_peak is not None:
        ax.plot(obs_peak.index, obs_peak.values, color=COLOR_OBS, lw=1.6,
                label="USGS obs (gauge 03463300)")
        peak_val = float(obs_peak.max())
        ax.axhline(peak_val, color=COLOR_OBS, lw=0.7, ls=":", alpha=0.5)
        ax.annotate(f"USGS peak\n{peak_val:.0f} m³/s",
                    xy=(HELENE_PEAK, peak_val), xytext=(10, -40),
                    textcoords="offset points", fontsize=8, color=COLOR_OBS)
    ax.set_ylabel("Discharge (m³/s)")
    ax.set_xlabel("Date")
    ax.set_title(f"Helene peak zoom | {args.label}")
    ax.legend(fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %HUTC"))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    plt.tight_layout()
    out = out_dir / "routed_ensemble_vs_usgs_peak.png"
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
