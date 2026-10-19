"""
plot_hydro_arm_with_openloop_f5.py

Plots the 2b hydro-state arm (DA on, 20 members: initial state uncertainty)
for cat-1016300 during Hurricane Helene (Sep 24-30 2024), with open-loop
grand median overlaid as a black dashed line.

Reads:
  - <arms-dir>/<cat>/<cat>_da_hydro_arm.csv  (all 21 catchments, mm/h)
  - <route-dir>/routed_leadtime_openloop_full.parquet  (open loop m³/s at gauge)
  - USGS obs CSV

Routes arm CSVs through T-Route Muskingum-Cunge to gauge 03463300,
one init time per Helene day (Sep 24-30), 20 members each.

Run on server (troute env):
    python3 plot_hydro_arm_with_openloop_f5.py
"""

import argparse
import os
import sys
import sqlite3
import types as _types
from functools import partial

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
DEFAULT_ARMS_DIR  = "/mnt/disk2/suma_helen_poster/da_results/da_arms_f5_rekrig_20pct"
DEFAULT_ROUTE_DIR = ("/mnt/disk2/suma_helen_poster/da_results/"
                     "da_forecast_f5_rekrig_20pct/routed")
DEFAULT_USGS_CSV  = "/mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv"
DEFAULT_GPKG      = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
                     "gage-03463300_subset.gpkg")
DEFAULT_OUT_DIR   = None   # falls back to --arms-dir

CATS = [
    'cat-1016279', 'cat-1016280', 'cat-1016281', 'cat-1016282', 'cat-1016283',
    'cat-1016300', 'cat-1016301', 'cat-1016302', 'cat-1016303', 'cat-1016304',
    'cat-1016305', 'cat-1016306', 'cat-1016307', 'cat-1016308', 'cat-1016309',
    'cat-1016310', 'cat-1016311', 'cat-1016312', 'cat-1016313', 'cat-1016314',
    'cat-1016315',
]

TERMINAL_INT       = 1016283   # wb-1016283 = gauge 03463300
DT                 = 3600.0
QTS_SUBDIVISIONS   = 1
WATERSHED_AREA_KM2 = 113.18   # total drainage area to gauge 03463300

HELENE_START = pd.Timestamp("2024-09-24 00:00:00")
HELENE_END   = pd.Timestamp("2024-09-30 23:00:00")

INIT_COLORS = [
    "#1f77b4",  # Sep-24
    "#ff7f0e",  # Sep-25
    "#2ca02c",  # Sep-26
    "#d62728",  # Sep-27
    "#9467bd",  # Sep-28
    "#8c564b",  # Sep-29
    "#e377c2",  # Sep-30
]


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
    link_set = {int(r[3:]) for r in fp_attr['link']}  # "wb-XXXX"
    connections = {}
    for _, row in fp_attr.iterrows():
        us = int(row['link'][3:])
        ds = int(row['to'][4:])            # "nex-XXXX"
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
    """Route one member's qlat_arr (n_segs × nts) through MC; return Q at terminal."""
    e1i  = np.zeros(0, dtype='int32')
    e1f  = np.zeros(0, dtype='float32')
    e2f  = np.zeros((0, nts), dtype='float32')
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
    """Route all 20 members for one issue time. Returns DataFrame index=valid_time."""
    n_segs   = len(param_df)
    n_leads  = len(lead_hours)
    n_members = len(member_cols)

    qlat_cube = np.zeros((n_members, n_segs, n_leads), dtype='float32')
    for cat, df in arm_dfs.items():
        sid     = int(cat[4:])            # "cat-XXXX"
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
        obs = obs * WATERSHED_AREA_KM2 * 1000.0 / 3600.0   # mm/h → m³/s
    return obs


# ── Open loop helper ──────────────────────────────────────────────────────────

def load_ol_grand_median(route_dir):
    """Return open-loop grand median (m³/s) indexed by valid_time for Helene."""
    path  = os.path.join(route_dir, 'routed_leadtime_openloop_full.parquet')
    df_ol = pd.read_parquet(path)
    df_ol['issue_time'] = pd.to_datetime(df_ol['issue_time'])

    helene = df_ol[(df_ol['issue_time'] >= HELENE_START) &
                   (df_ol['issue_time'] <= HELENE_END)]

    if 'q_gauge_m3s' in helene.columns:
        helene = helene.copy()
        helene['valid_time'] = (helene['issue_time'] +
                                pd.to_timedelta(helene['lead_hour'], unit='h'))
        grand = helene.groupby('valid_time')['q_gauge_m3s'].median()
    else:
        mem_cols = sorted(c for c in helene.columns if c.startswith('member_'))
        helene = helene.copy()
        helene['ens_mean'] = helene[mem_cols].mean(axis=1)
        helene['valid_time'] = (helene['issue_time'] +
                                pd.to_timedelta(helene['lead_hour'], unit='h'))
        grand = helene.groupby('valid_time')['ens_mean'].median()

    return grand.clip(lower=0).sort_index()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--arms-dir',  default=DEFAULT_ARMS_DIR)
    parser.add_argument('--route-dir', default=DEFAULT_ROUTE_DIR)
    parser.add_argument('--usgs-csv',  default=DEFAULT_USGS_CSV)
    parser.add_argument('--gpkg',      default=DEFAULT_GPKG)
    parser.add_argument('--out-dir',   default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    out_dir = args.out_dir or args.arms_dir
    os.makedirs(out_dir, exist_ok=True)

    # ── Load hydro arm CSVs ───────────────────────────────────────────────────
    print("Loading hydro arm CSVs...")
    arm_dfs = {}
    for cat in CATS:
        p = os.path.join(args.arms_dir, cat, f"{cat}_da_hydro_arm.csv")
        if os.path.exists(p):
            df = pd.read_csv(p)
            df['issue_time'] = pd.to_datetime(df['issue_time'])
            arm_dfs[cat] = df
    if 'cat-1016300' not in arm_dfs:
        raise FileNotFoundError(f"cat-1016300_da_hydro_arm.csv not found in {args.arms_dir}")
    print(f"  Loaded {len(arm_dfs)}/21 catchments")

    ref_df      = arm_dfs['cat-1016300']
    member_cols = sorted(c for c in ref_df.columns if c.startswith('member_'))
    all_issues  = pd.to_datetime(sorted(ref_df['issue_time'].unique()))
    lead_hours  = sorted(ref_df['lead_hour'].unique())
    print(f"  {len(all_issues)} issue times | {len(lead_hours)} leads | {len(member_cols)} members")

    # ── Build T-Route network (once) ──────────────────────────────────────────
    print("Building T-Route network from GPKG...")
    fp_attr, fp = read_network(args.gpkg)
    connections = build_connections(fp_attr)
    rconn       = nhd_network.reverse_network(connections)
    path_func   = partial(nhd_network.split_at_junction, rconn)
    reach_list  = nhd_network.dfs_decomposition(rconn, path_func)
    reaches_wTypes = [(r, 0) for r in reach_list]
    upstreams   = dict(rconn)
    param_df    = build_param_df(fp_attr)
    n_segs      = len(param_df)
    seg_ids_sorted = param_df.index.values

    area_map = {int(r['divide_id'][4:]): r['areasqkm'] * 1e6
                for _, r in fp.iterrows()
                if r['divide_id'] and str(r['divide_id']).startswith('cat-')}

    terminal_pos = int(np.where(seg_ids_sorted == TERMINAL_INT)[0][0])
    q0_df = pd.DataFrame(np.zeros((n_segs, 3), dtype='float32'),
                         index=param_df.index, columns=['qu0', 'qd0', 'h0'])
    print(f"  {n_segs} segments | terminal wb-{TERMINAL_INT} at index {terminal_pos}")

    # ── Load open loop + USGS ─────────────────────────────────────────────────
    print("Loading open loop parquet...")
    ol_grand = load_ol_grand_median(args.route_dir)
    print("Loading USGS obs...")
    usgs = load_usgs(args.usgs_csv)

    # ── Route one init time per Helene day + plot ─────────────────────────────
    print("Routing arms through T-Route (7 init times × 20 members)...")
    day_range = pd.date_range(HELENE_START.normalize(), HELENE_END.normalize(), freq='D')

    HELENE_PEAK_START = pd.Timestamp("2024-09-26 12:00:00")
    HELENE_PEAK_END   = pd.Timestamp("2024-09-28 00:00:00")

    fig, ax = plt.subplots(figsize=(17, 6))
    ax.axvspan(HELENE_PEAK_START, HELENE_PEAK_END,
               color='salmon', alpha=0.12, zorder=0, label='_nolegend_')
    ax.text(HELENE_PEAK_START + pd.Timedelta(hours=6), 1.0, "Helene peak",
            transform=ax.get_xaxis_transform(),
            fontsize=9, color='firebrick', ha='left', va='top')

    legend_patches = []

    for i, day in enumerate(day_range):
        color = INIT_COLORS[i % len(INIT_COLORS)]
        diffs = np.abs((all_issues - day).total_seconds())
        t0    = all_issues[diffs.argmin()]
        print(f"  Init {day.date()} → matched {t0} ...")

        piv = route_issue_time(
            t0, arm_dfs, member_cols, lead_hours,
            reaches_wTypes, upstreams, param_df, q0_df,
            seg_ids_sorted, area_map, terminal_pos)

        vt  = piv.index
        mn  = piv.min(axis=1)
        mx  = piv.max(axis=1)
        med = piv.median(axis=1)

        for col in piv.columns:
            ax.plot(vt, piv[col].values, color=color, lw=0.4, alpha=0.18, zorder=2)
        ax.fill_between(vt, mn, mx, color=color, alpha=0.22, linewidth=0, zorder=2)
        ax.plot(vt, med.values, color=color, lw=1.8, alpha=0.92, zorder=3)

        legend_patches.append(mpatches.Patch(color=color, alpha=0.85,
                                             label=f"Init {day.strftime('%Y-%m-%d')}"))

    # ── Open loop grand median ────────────────────────────────────────────────
    ol_mask = ((ol_grand.index >= HELENE_START) &
               (ol_grand.index <= HELENE_END + pd.Timedelta(hours=24)))
    ax.plot(ol_grand[ol_mask].index, ol_grand[ol_mask].values,
            color='black', lw=2.0, linestyle='--', zorder=5,
            label="Open loop (no DA) — grand median")

    # ── USGS obs ──────────────────────────────────────────────────────────────
    obs_mask = ((usgs.index >= HELENE_START) &
                (usgs.index <= HELENE_END + pd.Timedelta(hours=24)))
    ax.plot(usgs[obs_mask].index, usgs[obs_mask].values,
            color='black', lw=2.6, zorder=6, label="USGS obs")

    # ── Legend ────────────────────────────────────────────────────────────────
    extra_lines = [
        plt.Line2D([0], [0], color='black', lw=2.0, linestyle='--',
                   label="Open loop (no DA) — grand median"),
        plt.Line2D([0], [0], color='black', lw=2.6, label="USGS obs"),
    ]
    ax.legend(handles=legend_patches + extra_lines,
              fontsize=8, loc='upper left', framealpha=0.92, ncol=2)

    ax.set_ylabel("Discharge (m³/s)", fontsize=11)
    ax.set_xlabel("Date (UTC)", fontsize=11)
    ax.set_title(
        "2b — Hydro-state arm (DA on, 20 members): initial state uncertainty  |  "
        "Sep 24–30 2024\n"
        "cat-1016300  |  Shaded = member min–max  |  Line = median per init time",
        fontsize=11)
    ax.set_xlim(HELENE_START, HELENE_END + pd.Timedelta(hours=24))
    ax.grid(True, alpha=0.22)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha='right', fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "cat-1016300_2b_hydro_arm_helene.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
