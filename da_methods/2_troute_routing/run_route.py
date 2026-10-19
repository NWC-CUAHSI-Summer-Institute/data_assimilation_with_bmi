#!/usr/bin/env python3
"""
run_route.py
Route a single deterministic DA trajectory through T-route Muskingum-Cunge
to gauge 03463300 (South Toe River Near Celo, NC).

Reads per-catchment *_test_results.csv files (sim_mm_h column) from a DA
results directory, routes the full timeseries through the channel network,
and writes routed_Q_test.csv at the gauge outlet.

Works for any DA run that produces the standard _test_results.csv layout:
    columns: date, sim_mm_h, obs_mm_h, precip_mm_h

Usage:
    /home/svyas/miniconda3/envs/troute/bin/python run_route.py \\
        --gpkg    /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg \\
        --da-dir  /mnt/disk2/1400_sites_helene/da_results_dynamic_vrugt_seeded \\
        --out-dir /mnt/disk2/suma_helen_poster/da_results/dynamic_vrugt_seeded_routed \\
        --usgs-csv /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv \\
        --kv-dir  /mnt/disk2/1400_sites_helene/catchment_ts_03463300_spliced_dyn_helene

    # For the no-Vrugt (raw variance) seeded run:
    /home/svyas/miniconda3/envs/troute/bin/python run_route.py \\
        --gpkg    /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg \\
        --da-dir  /mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded \\
        --out-dir /mnt/disk2/suma_helen_poster/da_results/dynamic_novrugt_seeded_routed \\
        --usgs-csv /mnt/disk2/suma_helen_poster/03463300_usgs_hourly_2018_2024.csv \\
        --kv-dir  /mnt/disk2/1400_sites_helene/catchment_ts_03463300_spliced_dyn_helene
"""

import argparse
import os
import sys
import sqlite3
import types as _types
from functools import partial

import numpy as np
import pandas as pd

# Coverage stub — required for numba/troute on Python 3.10
_stub = _types.ModuleType('coverage.types')
for _cls in ['Tracer', 'TTraceData', 'TShouldTraceFn', 'TFileDisposition',
             'TShouldStartContextFn', 'TWarnFn', 'TTraceFn']:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules['coverage.types'] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

TERMINAL_INT     = 1016283   # wb-1016283 → gauge 03463300
DT               = 3600.0    # 1-hour timestep (seconds)
QTS_SUBDIVISIONS = 1
WATERSHED_AREA_KM2 = 113.18

CATS = [
    'cat-1016279', 'cat-1016280', 'cat-1016281', 'cat-1016282', 'cat-1016283',
    'cat-1016300', 'cat-1016301', 'cat-1016302', 'cat-1016303', 'cat-1016304',
    'cat-1016305', 'cat-1016306', 'cat-1016307', 'cat-1016308', 'cat-1016309',
    'cat-1016310', 'cat-1016311', 'cat-1016312', 'cat-1016313', 'cat-1016314',
    'cat-1016315',
]


def seg_int(wb_id):
    return int(wb_id[3:])


def read_network(gpkg):
    con = sqlite3.connect(gpkg)
    fp_attr = pd.read_sql(
        'SELECT link, "to", BtmWdth, TopWdth, TopWdthCC, n, nCC, ChSlp, So, Length_m '
        'FROM "flowpath-attributes"', con)
    fp = pd.read_sql('SELECT divide_id, areasqkm FROM flowpaths', con)
    con.close()
    return fp_attr, fp


def build_connections(fp_attr):
    link_set = {seg_int(r) for r in fp_attr['link']}
    connections = {}
    for _, row in fp_attr.iterrows():
        us = seg_int(row['link'])
        ds = int(row['to'][4:])
        connections[us] = [ds] if ds in link_set else []
    return connections


def build_param_df(fp_attr):
    rows = [{
        'seg_id': seg_int(r['link']),
        'dt':  float(DT),
        'bw':  float(r['BtmWdth']),   'tw':   float(r['TopWdth']),
        'twcc':float(r['TopWdthCC']), 'dx':   float(r['Length_m']),
        'n':   float(r['n']),          'ncc':  float(r['nCC']),
        'cs':  float(r['ChSlp']),      's0':   float(r['So']),
        'alt': 0.0,
    } for _, r in fp_attr.iterrows()]
    df = pd.DataFrame(rows).set_index('seg_id').sort_index()
    return df.astype('float32')


def route_timeseries(reaches_wTypes, upstreams, param_df, q0_df,
                     qlat_arr, nts, terminal_pos):
    """Route qlat_arr (n_segs × nts, m³/s) and return Q at terminal (nts,)."""
    e1i    = np.zeros(0, dtype='int32')
    e1f    = np.zeros(0, dtype='float32')
    e2f    = np.zeros((0, nts), dtype='float32')
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
        '2023-10-01_00:00:00',
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


def load_usgs(usgs_csv):
    if not usgs_csv or not os.path.exists(usgs_csv):
        return None
    df = pd.read_csv(usgs_csv)
    date_col = next(c for c in df.columns
                    if c.lower() in ('datetime', 'date', 'time', 'timestamp'))
    q_col = next(c for c in df.columns
                 if any(k in c.lower() for k in ('q', 'flow', 'discharge')))
    df[date_col] = pd.to_datetime(df[date_col])
    series = df.set_index(date_col)[q_col].astype(float)
    if 'mm' in q_col.lower():
        series = series * (WATERSHED_AREA_KM2 * 1000.0 / 3600.0)
    return series


def load_krig(kv_dir, dates):
    """Load Qkrig (mm/h) from the outlet catchment obs CSV and convert to m³/s."""
    for fname in ('cat-1016300.csv',):
        p = os.path.join(kv_dir, fname)
        if not os.path.exists(p):
            print(f'  WARNING: Qkrig file not found: {p}')
            return None
        df = pd.read_csv(p)
        df.columns = [c.strip() for c in df.columns]
        # Accept any date-like first column
        date_col = next((c for c in df.columns
                         if c.lower() in ('datetime', 'date', 'time', 'timestamp')),
                        df.columns[0])
        df[date_col] = pd.to_datetime(df[date_col])
        # Accept qkrig, qkrig_mm_hr, or any qkrig column without 'var'
        q_col = next((c for c in df.columns
                      if 'qkrig' in c.lower() and 'var' not in c.lower()), None)
        if q_col is None:
            print(f'  WARNING: no qkrig column in {p}. Columns: {list(df.columns)}')
            return None
        df = df.set_index(date_col).sort_index()
        area_m2 = 113.18 * 1e6   # full watershed — Qkrig at outlet catchment
        krig_m3s = df[q_col].astype(float) / 1000.0 / 3600.0 * area_m2
        return krig_m3s.reindex(dates)
    return None


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpkg',     required=True,
                        help='GeoPackage with flowpath-attributes and flowpaths tables')
    parser.add_argument('--da-dir',   required=True,
                        help='Dir containing <cat>/<cat>_test_results.csv (sim_mm_h column)')
    parser.add_argument('--out-dir',  required=True,
                        help='Output directory for routed_Q_test.csv')
    parser.add_argument('--usgs-csv', default=None,
                        help='USGS obs CSV for KGE/NSE summary (optional)')
    parser.add_argument('--kv-dir',   default=None,
                        help='Obs dir with per-catchment Qkrig CSVs (optional)')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Network ────────────────────────────────────────────────────────────
    print('Reading network from GPKG...')
    fp_attr, fp = read_network(args.gpkg)
    connections  = build_connections(fp_attr)
    rconn        = nhd_network.reverse_network(connections)
    path_func    = partial(nhd_network.split_at_junction, rconn)
    reach_list   = nhd_network.dfs_decomposition(rconn, path_func)
    reaches_wTypes = [(r, 0) for r in reach_list]
    upstreams    = dict(rconn)
    param_df     = build_param_df(fp_attr)
    n_segs       = len(param_df)

    area_map = {int(r['divide_id'][4:]): r['areasqkm'] * 1e6
                for _, r in fp.iterrows()
                if r['divide_id'] and str(r['divide_id']).startswith('cat-')}
    print(f'  {n_segs} segments | {len(area_map)} catchment areas loaded')

    seg_ids_sorted = param_df.index.values
    terminal_pos   = int(np.where(seg_ids_sorted == TERMINAL_INT)[0][0])
    print(f'  Terminal wb-{TERMINAL_INT} at index {terminal_pos}')

    q0_df = pd.DataFrame(
        np.zeros((n_segs, 3), dtype='float32'),
        index=param_df.index,
        columns=['qu0', 'qd0', 'h0'])

    # ── Load per-catchment sim_mm_h ────────────────────────────────────────
    print('Loading DA outputs...')
    cat_series = {}
    missing    = []
    for cat in CATS:
        sid = int(cat[4:])
        p   = os.path.join(args.da_dir, cat, f'{cat}_test_results.csv')
        if os.path.exists(p):
            df = pd.read_csv(p, parse_dates=['date'])
            df = df.set_index('date').sort_index()
            cat_series[sid] = df['sim_mm_h'].astype(float)
            continue
        # Fallback 1: use ensemble median from production per-member CSV
        p_prod = os.path.join(args.da_dir, cat, f'{cat}_production_per_member.csv')
        if os.path.exists(p_prod):
            df = pd.read_csv(p_prod, parse_dates=['date'])
            df = df.set_index('date').sort_index()
            mem_cols = [c for c in df.columns if c.startswith('member_')]
            cat_series[sid] = df[mem_cols].median(axis=1).astype(float)
            continue
        # Fallback 2: reconstruct nowcast from forcing arm lead_hour=1
        # valid_time = issue_time + 1h; covers only the Helene assimilation window
        p_arm = os.path.join(args.da_dir, cat, f'{cat}_da_forcing_arm.csv')
        if os.path.exists(p_arm):
            df = pd.read_csv(p_arm, parse_dates=['issue_time'])
            now = df[df['lead_hour'] == 1].copy()
            now['valid_time'] = now['issue_time'] + pd.Timedelta(hours=1)
            now = now.set_index('valid_time').sort_index()
            mem_cols = [c for c in now.columns if c.startswith('member_')]
            cat_series[sid] = now[mem_cols].median(axis=1).astype(float)
            continue
        missing.append(cat)

    if missing:
        print(f'  WARNING: missing outputs for {missing} — zero inflow assumed')
    print(f'  Loaded {len(cat_series)}/21 catchments')

    # Common time index from first available catchment
    ref_dates = next(iter(cat_series.values())).index
    nts       = len(ref_dates)
    print(f'  Timesteps: {nts}  ({ref_dates[0]} → {ref_dates[-1]})')

    # ── Build lateral inflow array (n_segs × nts, m³/s) ───────────────────
    print('Building lateral inflow array...')
    qlat = np.zeros((n_segs, nts), dtype='float32')
    seg_pos = {sid: i for i, sid in enumerate(seg_ids_sorted)}

    for sid, q_mm_h in cat_series.items():
        if sid not in seg_pos:
            continue
        area_m2 = area_map.get(sid, 0.0)
        if area_m2 == 0.0:
            print(f'  WARNING: no area for cat-{sid}')
            continue
        q_m3s = q_mm_h.reindex(ref_dates).fillna(0.0).values / 1000.0 / 3600.0 * area_m2
        qlat[seg_pos[sid], :] = q_m3s.astype('float32')

    # ── Route ──────────────────────────────────────────────────────────────
    print(f'Routing {nts} timesteps through T-route...')
    Q_routed = route_timeseries(
        reaches_wTypes, upstreams, param_df, q0_df, qlat, nts, terminal_pos)
    print(f'  Done. Peak routed Q = {Q_routed.max():.2f} m³/s')

    # ── Build output DataFrame ─────────────────────────────────────────────
    out_df = pd.DataFrame({'date': ref_dates, 'Q_routed_m3s': Q_routed})
    out_df = out_df.set_index('date')

    # Attach USGS obs
    usgs = load_usgs(args.usgs_csv)
    if usgs is not None:
        out_df['Q_usgs_m3s'] = usgs.reindex(ref_dates)

    # Attach Qkrig
    if args.kv_dir:
        krig = load_krig(args.kv_dir, ref_dates)
        if krig is not None:
            out_df['Q_krig_m3s'] = krig.values

    out_path = os.path.join(args.out_dir, 'routed_Q_test.csv')
    out_df.to_csv(out_path)
    print(f'Saved: {out_path}')

    # ── KGE / NSE summary ─────────────────────────────────────────────────
    if 'Q_usgs_m3s' in out_df.columns:
        obs = out_df['Q_usgs_m3s'].values
        sim = out_df['Q_routed_m3s'].values

        helene = ((out_df.index >= '2024-09-24') &
                  (out_df.index <= '2024-09-29 23:00:00'))

        print('\n── Routed vs USGS ─────────────────────────────────────────')
        print(f'  Full period :  KGE={kge(obs, sim):+.3f}  NSE={nse(obs, sim):+.3f}'
              f'  peak_sim={sim.max():.1f}  peak_obs={np.nanmax(obs):.1f} m³/s')
        if helene.sum() > 0:
            oh, sh = obs[helene], sim[helene]
            print(f'  Helene window:  KGE={kge(oh, sh):+.3f}  NSE={nse(oh, sh):+.3f}'
                  f'  peak_sim={sh.max():.1f}  peak_obs={np.nanmax(oh):.1f} m³/s'
                  f'  ({sh.max()/np.nanmax(oh)*100:.0f}% of USGS)')

        if 'Q_krig_m3s' in out_df.columns:
            krig_v = out_df['Q_krig_m3s'].values
            print(f'\n── Qkrig-routed vs USGS ───────────────────────────────────')
            print(f'  Full period :  KGE={kge(obs, krig_v):+.3f}  '
                  f'NSE={nse(obs, krig_v):+.3f}')
            if helene.sum() > 0:
                ok, sk = obs[helene], krig_v[helene]
                print(f'  Helene window:  KGE={kge(ok, sk):+.3f}  '
                      f'NSE={nse(ok, sk):+.3f}')
        print()


if __name__ == '__main__':
    main()
