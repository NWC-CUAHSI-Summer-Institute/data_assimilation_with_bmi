#!/usr/bin/env python3
"""
run_route_leadtime_forecasts.py
Route per-catchment lead-time forecast CSVs through T-route to gauge 03463300.

Reads the DA and open-loop lead-time forecast CSVs produced by
run_lead_time_forecast_sweep.py (one per catchment), routes all 20 members
through Muskingum-Cunge, and writes two output parquets at the gauge outlet.

Input (per catchment, from run_lead_time_forecast_sweep.py):
    <forecast-dir>/<cat-id>/<cat-id>_lead_time_forecasts_da.csv
    <forecast-dir>/<cat-id>/<cat-id>_lead_time_forecasts_openloop.csv
    Columns: issue_time, lead_hour, valid_time, member_00..member_19  (q in mm/h)

Output:
    <out-dir>/routed_leadtime_da_full.parquet
    <out-dir>/routed_leadtime_openloop_full.parquet
    Columns: issue_time, lead_hour, member_00..member_19  (Q in m³/s at gauge)

Usage:
    python3 run_route_leadtime_forecasts.py \\
        --gpkg /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg \\
        --forecast-dir /mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt \\
        --out-dir /mnt/disk2/suma_helen_poster/da_results_1gauge_heldout/folder1_vrugt_leadtime_routed
"""

import argparse
import sqlite3
import os
import sys
import types as _types
from functools import partial

import numpy as np
import pandas as pd

# Patch coverage.types for numba/troute on Python 3.10
_stub = _types.ModuleType('coverage.types')
for _cls in ['Tracer','TTraceData','TShouldTraceFn','TFileDisposition',
             'TShouldStartContextFn','TWarnFn','TTraceFn']:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules['coverage.types'] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

TERMINAL_INT   = 1016283   # wb-1016283 = gauge 03463300
DT             = 3600.0    # 1-hour timestep (seconds)
QTS_SUBDIVISIONS = 1

CATS = [
    'cat-1016279','cat-1016280','cat-1016281','cat-1016282','cat-1016283',
    'cat-1016300','cat-1016301','cat-1016302','cat-1016303','cat-1016304',
    'cat-1016305','cat-1016306','cat-1016307','cat-1016308','cat-1016309',
    'cat-1016310','cat-1016311','cat-1016312','cat-1016313','cat-1016314',
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
        'dt': float(DT),
        'bw': float(r['BtmWdth']),  'tw': float(r['TopWdth']),
        'twcc': float(r['TopWdthCC']), 'dx': float(r['Length_m']),
        'n':  float(r['n']),         'ncc': float(r['nCC']),
        'cs': float(r['ChSlp']),     's0':  float(r['So']),
        'alt': 0.0,
    } for _, r in fp_attr.iterrows()]
    df = pd.DataFrame(rows).set_index('seg_id').sort_index()
    return df.astype('float32')


def route_one(reaches_wTypes, upstreams, param_df, q0_df,
              qlat_arr, nts, n_segs, terminal_pos):
    """Route a single member's qlat through t-route; return Q at terminal."""
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
    Q_var = fvd[terminal_pos, :nts]
    Q_time = fvd[terminal_pos, 0::3]
    return Q_var if Q_var.max() > Q_time.max() else Q_time


def route_scenario(label, cat_data, seg_ids_sorted, area_map,
                   reaches_wTypes, upstreams, param_df, q0_df,
                   member_cols, issue_times, n_leads, n_members, n_segs,
                   terminal_pos, out_path):
    """Route one scenario (da or openloop) and save output parquet."""
    all_records = []

    for t_idx, t0 in enumerate(issue_times):
        t0_str = t0.strftime('%Y-%m-%d %H:%M:%S')

        qlat_cube = np.zeros((n_members, n_segs, n_leads), dtype='float32')

        for sid, df in cat_data.items():
            seg_pos = int(np.where(seg_ids_sorted == sid)[0][0]) if sid in seg_ids_sorted else -1
            if seg_pos < 0:
                continue
            area = area_map.get(sid, 0.0)

            sub = df[df['issue_time'] == t0].sort_values('lead_hour')
            if sub.empty:
                continue

            mem_vals = sub[member_cols].to_numpy(dtype='float32')  # (n_leads, n_members)
            mem_vals = mem_vals / 1000.0 / 3600.0 * area           # mm/h -> m³/s
            qlat_cube[:, seg_pos, :mem_vals.shape[0]] = mem_vals.T

        Q_members = np.full((n_leads, n_members), np.nan, dtype='float32')
        for m in range(n_members):
            qlat_m = qlat_cube[m]  # (n_segs, n_leads)
            qlat_df_m = pd.DataFrame(
                qlat_m, index=param_df.index,
                columns=range(n_leads))
            qlat_df_m = qlat_df_m.reindex(param_df.index, fill_value=0.0)
            Q_t = route_one(reaches_wTypes, upstreams, param_df, q0_df,
                            qlat_df_m.values.astype('float32'),
                            n_leads, n_segs, terminal_pos)
            Q_members[:, m] = Q_t[:n_leads]

        for lead in range(1, n_leads + 1):
            row = {'issue_time': t0_str, 'lead_hour': lead}
            for m, col in enumerate(member_cols):
                row[col] = float(Q_members[lead - 1, m])
            all_records.append(row)

        if (t_idx + 1) % 20 == 0 or (t_idx + 1) == len(issue_times):
            print(f'  [{label}] Routed {t_idx + 1}/{len(issue_times)} issue times')

    df_out = pd.DataFrame(all_records)
    df_out.to_parquet(out_path, index=False)
    print(f'Saved: {out_path}  shape={df_out.shape}')
    return df_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpkg',         required=True,
                        help='GPKG hydrofabric (gage-03463300_subset.gpkg)')
    parser.add_argument('--forecast-dir', required=True,
                        help='Dir with <cat>/<cat>_lead_time_forecasts_da.csv etc.')
    parser.add_argument('--out-dir',      required=True,
                        help='Output dir for routed parquets')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print('Reading network from GPKG...')
    fp_attr, fp = read_network(args.gpkg)
    connections = build_connections(fp_attr)
    rconn = nhd_network.reverse_network(connections)
    path_func = partial(nhd_network.split_at_junction, rconn)
    reach_list = nhd_network.dfs_decomposition(rconn, path_func)
    reaches_wTypes = [(r, 0) for r in reach_list]
    upstreams = dict(rconn)
    param_df = build_param_df(fp_attr)
    n_segs = len(param_df)

    area_map = {int(r['divide_id'][4:]): r['areasqkm'] * 1e6
                for _, r in fp.iterrows()
                if r['divide_id'] and str(r['divide_id']).startswith('cat-')}
    print(f'  {len(fp_attr)} segments | total area {sum(area_map.values())/1e6:.2f} km2')

    seg_ids_sorted = param_df.index.values
    terminal_pos = int(np.where(seg_ids_sorted == TERMINAL_INT)[0][0])
    print(f'  Terminal: wb-{TERMINAL_INT} at param_df index {terminal_pos}')

    q0_df = pd.DataFrame(
        np.zeros((n_segs, 3), dtype='float32'),
        index=param_df.index,
        columns=['qu0', 'qd0', 'h0'])

    for scenario in ('da', 'openloop'):
        print(f'\nLoading per-catchment lead-time CSVs ({scenario})...')
        cat_data = {}
        missing = []
        for cat in CATS:
            sid = int(cat[4:])
            csv_path = os.path.join(args.forecast_dir, cat,
                                    f'{cat}_lead_time_forecasts_{scenario}.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df['issue_time'] = pd.to_datetime(df['issue_time'])
                cat_data[sid] = df
            else:
                missing.append(cat)
        if missing:
            print(f'  WARNING: missing CSVs for: {missing}')
            print('  These catchments will have zero lateral inflow.')
        print(f'  Loaded {len(cat_data)}/21 catchments')

        ref_df = next(iter(cat_data.values()))
        issue_times = sorted(ref_df['issue_time'].unique())
        member_cols = sorted([c for c in ref_df.columns if c.startswith('member_')])
        n_members = len(member_cols)
        n_leads   = int(ref_df['lead_hour'].max())
        print(f'  {len(issue_times)} issue times | {n_leads} leads | {n_members} members')

        out_path = os.path.join(args.out_dir, f'routed_leadtime_{scenario}_full.parquet')
        route_scenario(
            scenario, cat_data, seg_ids_sorted, area_map,
            reaches_wTypes, upstreams, param_df, q0_df,
            member_cols, issue_times, n_leads, n_members, n_segs,
            terminal_pos, out_path)

    print('\nDone.')


if __name__ == '__main__':
    main()
