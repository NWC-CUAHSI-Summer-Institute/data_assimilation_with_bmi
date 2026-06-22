#!/usr/bin/env python3
import argparse
import os
import sqlite3
import sys
import time
import types as _types
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd

# coverage.types stub -- required so numba / t-route import under Python 3.10.
# MUST run before troute is imported.
_stub = _types.ModuleType("coverage.types")
for _cls in ["Tracer", "TTraceData", "TShouldTraceFn", "TFileDisposition",
             "TShouldStartContextFn", "TWarnFn", "TTraceFn"]:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules["coverage.types"] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

# ---------------------------------------------------------------------------
# Gauge / network constants (GAUGE-SPECIFIC -- change these to swap the gauge)
# ---------------------------------------------------------------------------
TERMINAL_INT       = 1016283     # outlet segment wb-<id> = gauge 03463300
WATERSHED_AREA_KM2 = 113.18      # total basin area (sum of the divides)
CATS = [
    "cat-1016279", "cat-1016280", "cat-1016281", "cat-1016282", "cat-1016283",
    "cat-1016300", "cat-1016301", "cat-1016302", "cat-1016303", "cat-1016304",
    "cat-1016305", "cat-1016306", "cat-1016307", "cat-1016308", "cat-1016309",
    "cat-1016310", "cat-1016311", "cat-1016312", "cat-1016313", "cat-1016314",
    "cat-1016315",
]

# Routing constants (not gauge-specific)
DT               = 3600.0        # 1-hour timestep (seconds)
QTS_SUBDIVISIONS = 1

# Default I/O (all overridable on the CLI)
FC_DIR      = "/mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_wide_spread"
DA_TRAJ_DIR = "/mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_da"
OL_TRAJ_DIR = "/mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_ol"
OUT_DIR     = "/mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_wide_spread/routed"
GPKG        = ("/mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/"
               "gage-03463300_subset.gpkg")

# globals shared with fork()ed workers (copy-on-write, read-only)
G = {}


# ===========================================================================
# 1. CHANNEL NETWORK  (build the Muskingum-Cunge network from the GeoPackage)
# ===========================================================================
def seg_int(wb_id):
    """'wb-1016283' -> 1016283."""
    return int(wb_id[3:])


def read_network(gpkg):
    """Read the flowpath-attributes (channel params + topology) and flowpaths
    (divide_id + area) tables from the GeoPackage."""
    con = sqlite3.connect(gpkg)
    fp_attr = pd.read_sql(
        'SELECT link, "to", BtmWdth, TopWdth, TopWdthCC, n, nCC, ChSlp, So, Length_m '
        'FROM "flowpath-attributes"', con)
    fp = pd.read_sql("SELECT divide_id, areasqkm FROM flowpaths", con)
    con.close()
    return fp_attr, fp


def build_connections(fp_attr):
    """upstream segment -> [downstream segment] (empty list at the outlet)."""
    link_set = {seg_int(r) for r in fp_attr["link"]}
    connections = {}
    for _, row in fp_attr.iterrows():
        us = seg_int(row["link"])
        ds = int(row["to"][4:])
        connections[us] = [ds] if ds in link_set else []
    return connections


def build_param_df(fp_attr):
    """Per-segment Muskingum-Cunge channel parameters, indexed by seg_id."""
    rows = [{
        "seg_id": seg_int(r["link"]),
        "dt":   float(DT),
        "bw":   float(r["BtmWdth"]),   "tw":   float(r["TopWdth"]),
        "twcc": float(r["TopWdthCC"]), "dx":   float(r["Length_m"]),
        "n":    float(r["n"]),         "ncc":  float(r["nCC"]),
        "cs":   float(r["ChSlp"]),     "s0":   float(r["So"]),
        "alt":  0.0,
    } for _, r in fp_attr.iterrows()]
    return pd.DataFrame(rows).set_index("seg_id").sort_index().astype("float32")


def build_network(gpkg):
    """Assemble everything route_flow() needs: reach decomposition, upstream map,
    param table, terminal position, per-catchment segment position, and the
    mm/h -> m^3/s factor per catchment (= divide_area_m2 / 1000 / 3600)."""
    fp_attr, fp = read_network(gpkg)
    connections = build_connections(fp_attr)
    rconn = nhd_network.reverse_network(connections)
    reach_list = nhd_network.dfs_decomposition(
        rconn, partial(nhd_network.split_at_junction, rconn))
    param_df = build_param_df(fp_attr)
    seg_ids = param_df.index.values

    area_map = {int(r["divide_id"][4:]): r["areasqkm"] * 1e6
                for _, r in fp.iterrows()
                if r["divide_id"] and str(r["divide_id"]).startswith("cat-")}
    net = {
        "reaches_wTypes": [(r, 0) for r in reach_list],
        "upstreams": dict(rconn),
        "param_df": param_df,
        "n_segs": len(param_df),
        "terminal_pos": int(np.where(seg_ids == TERMINAL_INT)[0][0]),
        "seg_pos": {int(c[4:]): int(np.where(seg_ids == int(c[4:]))[0][0]) for c in CATS},
        "factor": {int(c[4:]): area_map[int(c[4:])] / 1000.0 / 3600.0 for c in CATS},
        "q0_df": pd.DataFrame(np.zeros((len(param_df), 3), dtype="float32"),
                              index=param_df.index, columns=["qu0", "qd0", "h0"]),
    }
    return net


# ===========================================================================
# 2. ROUTING  (one member through the whole network)
# ===========================================================================
def route_flow(net, qlat_arr, nts):
    """Route one member's lateral inflow (n_segs x nts, m^3/s) through the network
    and return the terminal FLOW (m^3/s) for all nts steps.
    compute_network_structured returns (n_segs, nts*3) interleaved
    (flow, velocity, depth) per timestep; flow is the 0::3 stride."""
    e1i = np.zeros(0, dtype="int32")
    e1f = np.zeros(0, dtype="float32")
    e2f = np.zeros((0, nts), dtype="float32")
    e00f32 = np.zeros((0, 0), dtype="float32")
    e00f64 = np.zeros((0, 0), dtype="float64")
    e00i32 = np.zeros((0, 0), dtype="int32")
    res = compute_network_structured(
        nts, DT, QTS_SUBDIVISIONS,
        net["reaches_wTypes"], net["upstreams"],
        net["param_df"].index.values.astype("int64"),
        net["param_df"].columns.values,
        net["param_df"].values,
        net["q0_df"].values.astype("float32"),
        qlat_arr.astype("float32"),
        [], e00f64, {}, e00i32, False,
        "2024-09-20_00:00:00",
        e2f, e1i, e1i, e1i, e1f, e1f, 0.0,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1i, [], e1i, e1i, e1f, e1i, e1i,
        e1i, e1i, e1f, e1i, e1f, e1i, e1i, e00f32,
    )
    fvd = np.asarray(res[1])
    return fvd[net["terminal_pos"], 0::3]


def route_issue(k):
    """Route all members for issue index k. Returns (k, (n_leads, n_members)).

    Builds qlat (n_members, n_segs, warmup+leads): for each catchment segment, the
    warm-up block (real antecedent runoff in 'trajectory' mode, or a constant
    lead-1 proxy) followed by the 18 forecast leads; routes each member; slices
    the warm-up off, keeping the leads on a channel that's already full."""
    net = G["net"]
    n_segs, n_members = net["n_segs"], G["n_members"]
    n_leads, warmup, mode = G["n_leads"], G["warmup_h"], G["warmup_mode"]
    n_route = warmup + n_leads
    qlat = np.zeros((n_members, n_segs, n_route), dtype="float32")

    for sid in G["sids"]:
        pos = net["seg_pos"][sid]
        f = net["factor"][sid]
        leads = G["fc"][sid][k].T * f                    # (n_members, n_leads) m^3/s
        if warmup > 0:
            if mode == "trajectory":
                t_idx = G["traj_t0_idx"][k]              # row of t0 in the trajectory
                qlat[:, pos, :warmup] = G["traj"][sid][t_idx - warmup + 1:t_idx + 1].T * f
            else:  # constant lead-1 proxy
                qlat[:, pos, :warmup] = leads[:, 0:1]
        qlat[:, pos, warmup:] = leads

    out = np.empty((n_leads, n_members), dtype="float32")
    for m in range(n_members):
        Q = route_flow(net, qlat[m], n_route)
        out[:, m] = Q[warmup:warmup + n_leads]
    return k, out


# ===========================================================================
# 3. SCENARIO DRIVER  (load forecasts/trajectories, route, save)
# ===========================================================================
def run_scenario(tag, fc_name, traj_dir, net, args):
    t0 = time.time()
    print(f"[{tag}] loading forecasts ({args.warmup_mode} warm-up) for {len(CATS)} catchments...")
    member_cols, issues, n_leads = None, None, None
    fc, traj, traj_times = {}, {}, None
    for cat in CATS:
        sid = int(cat[4:])
        f = pd.read_parquet(os.path.join(args.fc_dir, cat,
                                         f"{cat}_forecast_leadtime_{fc_name}.parquet"))
        if member_cols is None:
            member_cols = [c for c in f.columns if c.startswith("member_")]
            issues = sorted(f["issue_time"].unique().tolist())
            n_leads = int(f["lead_hour"].max())
            if args.max_issues > 0:
                issues = issues[:args.max_issues]
        sub = f[f["issue_time"].isin(issues)].sort_values(["issue_time", "lead_hour"])
        # (n_issues, n_leads, n_members) -- member columns are the SAME order in every catchment
        fc[sid] = (sub[member_cols].values.astype("float32")
                   .reshape(len(issues), n_leads, len(member_cols)))
        if args.warmup_mode == "trajectory":
            tr = pd.read_parquet(os.path.join(traj_dir, cat, f"{cat}_member_q.parquet"))
            if traj_times is None:
                traj_times = tr["time"].astype(str).tolist()
            traj[sid] = tr[member_cols].values.astype("float32")   # (n_hours, n_members)

    G.update({
        "net": net, "fc": fc, "sids": [int(c[4:]) for c in CATS],
        "n_members": len(member_cols), "n_leads": n_leads,
        "warmup_h": args.warmup_h, "warmup_mode": args.warmup_mode,
    })
    if args.warmup_mode == "trajectory":
        t_row = {t: i for i, t in enumerate(traj_times)}
        missing = [t for t in issues if t not in t_row]
        if missing:
            raise SystemExit(f"[{tag}] {len(missing)} issue times absent from {traj_dir} "
                             f"member_q (e.g. {missing[:2]}); cannot trajectory-warm-up.")
        G["traj"] = traj
        G["traj_t0_idx"] = [t_row[t] for t in issues]

    print(f"[{tag}] routing {len(issues)} issues x {len(member_cols)} members "
          f"({args.warmup_h}h {args.warmup_mode} warm-up + {n_leads} leads each), "
          f"{args.procs} procs...")
    results = {}
    with Pool(args.procs) as pool:
        for i, (k, out) in enumerate(pool.imap_unordered(route_issue,
                                                         range(len(issues)), chunksize=2)):
            results[k] = out
            if (i + 1) % 24 == 0:
                print(f"  [{tag}] {i + 1}/{len(issues)} issues ({time.time() - t0:.0f}s)")

    meta, rows = [], []
    for k, t in enumerate(issues):
        for lead in range(1, n_leads + 1):
            meta.append((t, lead))
            rows.append(results[k][lead - 1])
    out_df = pd.concat([pd.DataFrame(meta, columns=["issue_time", "lead_hour"]),
                        pd.DataFrame(np.vstack(rows), columns=member_cols)], axis=1)
    os.makedirs(args.out_dir, exist_ok=True)
    pq = os.path.join(args.out_dir, f"routed_forecast_leadtime_{fc_name}.parquet")
    out_df.to_parquet(pq, index=False)
    print(f"[{tag}] saved {pq} ({len(out_df)} rows, {len(member_cols)} members) "
          f"in {time.time() - t0:.0f}s")
    G.pop("traj", None)
    G.pop("traj_t0_idx", None)


# ===========================================================================
# 4. CLI
# ===========================================================================
def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--scenario", default="both", choices=["da", "openloop", "both"])
    p.add_argument("--warmup-mode", default="trajectory", choices=["trajectory", "constant"],
                   help="trajectory = actual antecedent member_q (correct); "
                        "constant = hold lead-1 inflow (no trajectory files needed)")
    p.add_argument("--warmup-h", type=int, default=48)
    p.add_argument("--procs", type=int, default=20)
    p.add_argument("--max-issues", type=int, default=0,
                   help="route only the first N issue times (0 = all; for smoke tests)")
    p.add_argument("--gpkg", default=GPKG)
    p.add_argument("--fc-dir", default=FC_DIR)
    p.add_argument("--da-traj-dir", default=DA_TRAJ_DIR)
    p.add_argument("--ol-traj-dir", default=OL_TRAJ_DIR)
    p.add_argument("--out-dir", default=OUT_DIR)
    args = p.parse_args()

    print("Building channel network from GPKG...")
    net = build_network(args.gpkg)
    print(f"  {net['n_segs']} segments | terminal wb-{TERMINAL_INT} at pos {net['terminal_pos']}")

    if args.scenario in ("da", "both"):
        run_scenario("DA", "da", args.da_traj_dir, net, args)
    if args.scenario in ("openloop", "both"):
        run_scenario("OL", "openloop", args.ol_traj_dir, net, args)


if __name__ == "__main__":
    main()
