"""
plot_da_arms_routed_f5.py

Routes both the forcing arm (ARM A, 30 members) and the hydro-state arm
(ARM B, 20 members) through T-Route Muskingum-Cunge for all 21 catchments
draining to gauge 03463300, then plots the comparison against USGS obs.

This is the correctly routed version of the 2ab comparison figure.
Each catchment's mm/h values are converted using that catchment's own
area from the GPKG (not the total watershed area approximation).

Reads:
  - <arms-dir>/<cat>/<cat>_da_forcing_arm.csv   (all 21 catchments, mm/h)
  - <arms-dir>/<cat>/<cat>_da_hydro_arm.csv     (all 21 catchments, mm/h)
  - USGS obs CSV

Output:
  <out-dir>/cat-1016300_2ab_arms_comparison_routed.png

Run on server (troute env):
    python3 plot_da_arms_routed_f5.py
"""

import os
import sys
import sqlite3
import types as _types
from functools import partial
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates

# ── Coverage stub required for numba/troute on Python 3.10 ───────────────────
_stub = _types.ModuleType('coverage.types')
for _cls in ['Tracer', 'TTraceData', 'TShouldTraceFn', 'TFileDisposition',
             'TShouldStartContextFn', 'TWarnFn', 'TTraceFn']:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules['coverage.types'] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

# ── Constants ─────────────────────────────────────────────────────────────────
DEFAULT_ARMS_DIR = "/mnt/disk2/suma_helen_poster/da_results/da_arms_f5_rekrig_20pct"
DEFAULT_USGS_CSV = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
DEFAULT_GPKG     = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                    "gage-03463300_subset.gpkg")
DEFAULT_OUT_DIR  = None  # falls back to --arms-dir

CATS = [
    'cat-1016279', 'cat-1016280', 'cat-1016281', 'cat-1016282', 'cat-1016283',
    'cat-1016300', 'cat-1016301', 'cat-1016302', 'cat-1016303', 'cat-1016304',
    'cat-1016305', 'cat-1016306', 'cat-1016307', 'cat-1016308', 'cat-1016309',
    'cat-1016310', 'cat-1016311', 'cat-1016312', 'cat-1016313', 'cat-1016314',
    'cat-1016315',
]

TERMINAL_INT     = 1016283   # wb-1016283 = gauge 03463300
DT               = 3600.0
QTS_SUBDIVISIONS = 1
WATERSHED_AREA_KM2 = 113.18  # for USGS unit conversion only

HELENE_START      = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END        = pd.Timestamp("2024-09-30 23:00:00")
HELENE_PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
HELENE_PEAK_END   = pd.Timestamp("2024-09-28 00:00:00")

FORCING_COLOR = "#4878CF"   # blue
HYDRO_COLOR   = "#2ca02c"   # green


# ── T-Route helpers ───────────────────────────────────────────────────────────

def read_network(gpkg):
    con = sqlite3.connect(gpkg)
    fp_attr = pd.read_sql(
        'SELECT link, "to", BtmWdth, TopWdth, TopWdthCC, n, nCC, ChSlp, So, Length_m '
        'FROM "flowpath-attributes"', con)
    fp = pd.read_sql('SELECT divide_id, areasqkm FROM flowpaths', con)
    con.close()
    return fp_attr, fp


def build_connections(fp_attr):
    link_set = {int(r[3:]) for r in fp_attr['link']}
    connections = {}
    for _, row in fp_attr.iterrows():
        us = int(row['link'][3:])
        ds = int(row['to'][4:])
        connections[us] = [ds] if ds in link_set else []
    return connections


def build_param_df(fp_attr):
    rows = [{'seg_id': int(r['link'][3:]),
             'dt':    float(DT),
             'bw':    float(r['BtmWdth']),   'tw':   float(r['TopWdth']),
             'twcc':  float(r['TopWdthCC']),  'dx':   float(r['Length_m']),
             'n':     float(r['n']),           'ncc':  float(r['nCC']),
             'cs':    float(r['ChSlp']),       's0':   float(r['So']),
             'alt':   0.0}
            for _, r in fp_attr.iterrows()]
    return pd.DataFrame(rows).set_index('seg_id').sort_index().astype('float32')


def route_one(reaches_wTypes, upstreams, param_df, q0_df,
              qlat_arr, nts, terminal_pos):
    e1i   = np.zeros(0, dtype='int32')
    e1f   = np.zeros(0, dtype='float32')
    e2f   = np.zeros((0, nts), dtype='float32')
    e00f32 = np.zeros((0, 0), dtype='float32')
    e00f64 = np.zeros((0, 0), dtype='float64')
    e00i32 = np.zeros((0, 0), dtype='int32')

    results = compute_network_structured(
        nts, DT, QTS_SUBDIVISIONS,
        reaches_wTypes, upstreams,
        param_df.index.values.astype('int64'),
        param_df.columns.values,
        param_df.values,
        q0_df.values.astype('float32'),
        qlat_arr.astype('float32'),
        [], e00f64, {}, e00i32, False,
        '2024-09-24_00:00:00',
        e2f, e1i, e1i, e1i, e1f, e1f, 0.0,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1i, [], e1i, e1i, e1f, e1i, e1i,
        e1i, e1i, e1f, e1i, e1f, e1i, e1i, e00f32,
    )
    seg_ids = np.asarray(results[0])
    fvd     = np.asarray(results[1])
    Q_var   = fvd[terminal_pos, :nts]
    Q_time  = fvd[terminal_pos, 0::3]
    return Q_var if Q_var.max() > Q_time.max() else Q_time


def route_issue_time(t0, arm_dfs, member_cols, lead_hours,
                     reaches_wTypes, upstreams, param_df, q0_df,
                     seg_ids_sorted, area_map, terminal_pos):
    """Route all members for one issue time. Returns DataFrame index=valid_time."""
    n_segs    = len(param_df)
    n_leads   = len(lead_hours)
    n_members = len(member_cols)

    qlat_cube = np.zeros((n_members, n_segs, n_leads), dtype='float32')
    for cat, df in arm_dfs.items():
        sid     = int(cat[4:])
        area_m2 = area_map.get(sid, 0.0)
        if sid not in seg_ids_sorted:
            continue
        seg_pos = int(np.where(seg_ids_sorted == sid)[0][0])
        sub = df[df['issue_time'] == t0].sort_values('lead_hour')
        if sub.empty:
            continue
        vals = np.maximum(sub[member_cols].to_numpy(dtype='float32'), 0.0)
        vals = vals / 1000.0 / 3600.0 * area_m2   # mm/h → m³/s
        qlat_cube[:, seg_pos, :vals.shape[0]] = vals.T

    Q_members = np.full((n_leads, n_members), np.nan, dtype='float32')
    for m in range(n_members):
        Q_m = route_one(reaches_wTypes, upstreams, param_df, q0_df,
                        qlat_cube[m], n_leads, terminal_pos)
        Q_members[:, m] = Q_m[:n_leads]

    valid_times = pd.DatetimeIndex([t0 + pd.Timedelta(hours=int(h)) for h in lead_hours])
    return pd.DataFrame(np.maximum(Q_members, 0.0),
                        columns=member_cols, index=valid_times)


# ── Route all issue times for one arm; collect all Q values per valid_time ────

def route_arm_all_inits(arm_dfs, member_cols, all_issues, lead_hours,
                        reaches_wTypes, upstreams, param_df, q0_df,
                        seg_ids_sorted, area_map, terminal_pos, label):
    """Route every issue time for one arm.

    Returns a dict: valid_time (Timestamp) → 1-D array of all member Q values
    across all issue times that cover that valid_time.
    """
    q_by_vt = defaultdict(list)
    day_range = pd.date_range(HELENE_START.normalize(), HELENE_END.normalize(), freq='D')

    for day in day_range:
        diffs = np.abs((all_issues - day).total_seconds())
        t0    = all_issues[diffs.argmin()]
        print(f"  [{label}] init {day.date()} → matched {t0} ...")

        piv = route_issue_time(
            t0, arm_dfs, member_cols, lead_hours,
            reaches_wTypes, upstreams, param_df, q0_df,
            seg_ids_sorted, area_map, terminal_pos)

        for vt, row in piv.iterrows():
            q_by_vt[vt].extend(row.dropna().tolist())

    return q_by_vt


def summarise(q_by_vt):
    """From {valid_time: [Q values]}, return (times, median, q5, q95)."""
    times  = sorted(q_by_vt.keys())
    median = np.array([np.median(q_by_vt[t]) for t in times])
    q5     = np.array([np.percentile(q_by_vt[t],  5) for t in times])
    q95    = np.array([np.percentile(q_by_vt[t], 95) for t in times])
    return times, median, q5, q95


# ── USGS loader ───────────────────────────────────────────────────────────────

def load_usgs(usgs_csv):
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns
                    if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    q_col    = next(c for c in df.columns
                    if any(k in c.lower() for k in ('q', 'flow', 'discharge')))
    df[date_col] = pd.to_datetime(df[date_col])
    obs = df.set_index(date_col)[q_col].astype(float)
    if 'mm' in q_col.lower():
        obs = obs * WATERSHED_AREA_KM2 * 1000.0 / 3600.0
    return obs


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--arms-dir', default=DEFAULT_ARMS_DIR)
    parser.add_argument('--usgs-csv', default=DEFAULT_USGS_CSV)
    parser.add_argument('--gpkg',     default=DEFAULT_GPKG)
    parser.add_argument('--out-dir',  default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    out_dir = args.out_dir or args.arms_dir
    os.makedirs(out_dir, exist_ok=True)

    # ── Load forcing arm CSVs ─────────────────────────────────────────────────
    print("Loading forcing arm CSVs...")
    fa_dfs = {}
    for cat in CATS:
        p = os.path.join(args.arms_dir, cat, f"{cat}_da_forcing_arm.csv")
        if os.path.exists(p):
            df = pd.read_csv(p)
            df['issue_time'] = pd.to_datetime(df['issue_time'])
            fa_dfs[cat] = df
    print(f"  Loaded {len(fa_dfs)}/21 catchments (forcing arm)")

    fa_ref      = fa_dfs['cat-1016300']
    fa_members  = sorted(c for c in fa_ref.columns if c.startswith('member_'))
    all_issues  = pd.to_datetime(sorted(fa_ref['issue_time'].unique()))
    lead_hours  = sorted(fa_ref['lead_hour'].unique())
    print(f"  {len(all_issues)} issue times | {len(lead_hours)} leads | "
          f"{len(fa_members)} forcing members")

    # ── Load hydro arm CSVs ───────────────────────────────────────────────────
    print("Loading hydro arm CSVs...")
    ha_dfs = {}
    for cat in CATS:
        p = os.path.join(args.arms_dir, cat, f"{cat}_da_hydro_arm.csv")
        if os.path.exists(p):
            df = pd.read_csv(p)
            df['issue_time'] = pd.to_datetime(df['issue_time'])
            ha_dfs[cat] = df
    ha_ref     = ha_dfs['cat-1016300']
    ha_members = sorted(c for c in ha_ref.columns if c.startswith('member_'))
    print(f"  Loaded {len(ha_dfs)}/21 catchments (hydro arm) | "
          f"{len(ha_members)} hydro members")

    # ── Build T-Route network (once, shared for both arms) ────────────────────
    print("Building T-Route network from GPKG...")
    fp_attr, fp = read_network(args.gpkg)
    connections = build_connections(fp_attr)
    rconn       = nhd_network.reverse_network(connections)
    path_func   = partial(nhd_network.split_at_junction, rconn)
    reach_list  = nhd_network.dfs_decomposition(rconn, path_func)
    reaches_wTypes = [(r, 0) for r in reach_list]
    upstreams      = dict(rconn)
    param_df       = build_param_df(fp_attr)
    seg_ids_sorted = param_df.index.values
    n_segs         = len(param_df)

    area_map = {int(r['divide_id'][4:]): r['areasqkm'] * 1e6
                for _, r in fp.iterrows()
                if r['divide_id'] and str(r['divide_id']).startswith('cat-')}

    terminal_pos = int(np.where(seg_ids_sorted == TERMINAL_INT)[0][0])
    q0_df = pd.DataFrame(np.zeros((n_segs, 3), dtype='float32'),
                         index=param_df.index, columns=['qu0', 'qd0', 'h0'])
    print(f"  {n_segs} segments | terminal wb-{TERMINAL_INT} at index {terminal_pos}")

    # ── Route forcing arm (all issue times) ───────────────────────────────────
    print("\nRouting forcing arm (7 init times × 30 members)...")
    fa_q = route_arm_all_inits(
        fa_dfs, fa_members, all_issues, lead_hours,
        reaches_wTypes, upstreams, param_df, q0_df,
        seg_ids_sorted, area_map, terminal_pos, label="forcing")
    fa_times, fa_med, fa_q5, fa_q95 = summarise(fa_q)

    # ── Route hydro arm (all issue times) ─────────────────────────────────────
    print("\nRouting hydro arm (7 init times × 20 members)...")
    ha_q = route_arm_all_inits(
        ha_dfs, ha_members, all_issues, lead_hours,
        reaches_wTypes, upstreams, param_df, q0_df,
        seg_ids_sorted, area_map, terminal_pos, label="hydro")
    ha_times, ha_med, ha_q5, ha_q95 = summarise(ha_q)

    # ── USGS obs ──────────────────────────────────────────────────────────────
    print("\nLoading USGS obs...")
    usgs     = load_usgs(args.usgs_csv)
    obs_mask = ((usgs.index >= HELENE_START) &
                (usgs.index <= HELENE_END + pd.Timedelta(hours=24)))

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(17, 6))

    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color='salmon', alpha=0.12, zorder=0, label='_nolegend_')
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color='firebrick', ha='left', va='top')

    # Forcing arm band + median
    ax.fill_between(fa_times, fa_q5, fa_q95,
                    color=FORCING_COLOR, alpha=0.25, linewidth=0, zorder=2,
                    label=f"Forcing arm spread [5th–95th] (N={len(fa_members)} members)")
    ax.plot(fa_times, fa_med,
            color=FORCING_COLOR, lw=2.2, zorder=3,
            label="Forcing arm — grand median")

    # Hydro arm band + median
    ax.fill_between(ha_times, ha_q5, ha_q95,
                    color=HYDRO_COLOR, alpha=0.25, linewidth=0, zorder=2,
                    label=f"Hydro-state arm spread [5th–95th] (N={len(ha_members)} members)")
    ax.plot(ha_times, ha_med,
            color=HYDRO_COLOR, lw=2.2, zorder=3,
            label="Hydro-state arm — grand median")

    # USGS obs
    ax.plot(usgs[obs_mask].index, usgs[obs_mask].values,
            color='black', lw=2.6, zorder=6, label="USGS obs")

    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_title(
        "2a vs 2b — Forcing arm (blue) vs Hydro-state arm (green)  |  DA on  |  "
        "All 21 catchments T-Route routed to gauge 03463300\n"
        "Shaded = 5th–95th percentile across all issue times  |  "
        "Line = grand median  |  F5 re-kriged σ²  |  20% gauge holdout",
        fontsize=10)
    ax.set_xlim(HELENE_START, HELENE_END + pd.Timedelta(hours=24))
    ax.legend(fontsize=9, loc='upper left', framealpha=0.92)
    ax.grid(True, alpha=0.22)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha='right', fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "cat-1016300_2ab_arms_comparison_routed.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
