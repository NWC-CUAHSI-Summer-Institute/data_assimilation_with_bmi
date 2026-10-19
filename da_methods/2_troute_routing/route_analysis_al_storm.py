#!/usr/bin/env python3
"""
route_analysis_al_storm.py
Routes the full DA + OL analysis trajectory (Oct 1 - Oct 9) through T-Route
5-min Muskingum-Cunge for all 12 Alabama storm gauges.

Reads:  <da-traj-dir>/<cat>/<cat>_member_q.parquet  (mm/h, hourly per member)
        <ol-traj-dir>/<cat>/<cat>_member_q.parquet
Writes: <fc-dir>/routed/<gauge>/routed_analysis_da.parquet
        <fc-dir>/routed/<gauge>/routed_analysis_ol.parquet
Columns: time (str), member_0000..member_NNNN (m3/s)

These are prepended to the spaghetti fan plot to show DA correction Oct 1-Oct 9.
"""
import argparse
import sys
import time
import types as _types
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

# coverage.types stub — must run before troute import
_stub = _types.ModuleType("coverage.types")
for _cls in ["Tracer", "TTraceData", "TShouldTraceFn", "TFileDisposition",
             "TShouldStartContextFn", "TWarnFn", "TTraceFn"]:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules["coverage.types"] = _stub

sys.path.insert(0, str(Path(__file__).parent))
from route_forecast import build_network, route_flow, DT, SUBSTEPS

AL_GAUGES = [
    "02361000", "02369800", "02371500", "02372250", "02374500", "02408540",
    "02422500", "02450250", "02464000", "02465493", "02469800", "03574500",
]

G = {}


def _route_one_member(m):
    net = G["net"]
    qlat = G["qlat"]          # (n_members, n_segs, n_hours) m^3/s
    n_hours = qlat.shape[2]
    qf = np.repeat(qlat[m], SUBSTEPS, axis=1)   # (n_segs, n_hours*SUBSTEPS)
    Q  = route_flow(net, qf, n_hours * SUBSTEPS, DT / SUBSTEPS)
    return Q[SUBSTEPS - 1::SUBSTEPS][:n_hours]   # hourly terminal flow (m^3/s)


def route_scenario(tag, traj_dir, net, out_path, procs):
    cats   = net["cats"]
    n_segs = net["n_segs"]

    print(f"  [{tag}] loading trajectories for {len(cats)} catchments...")
    member_cols = None
    times       = None
    n_hours     = None
    qlat_parts  = {}

    for cat in cats:
        sid = int(cat[4:])
        pq  = Path(traj_dir) / cat / f"{cat}_member_q.parquet"
        if not pq.exists():
            continue
        tr = pd.read_parquet(pq)
        if member_cols is None:
            member_cols = [c for c in tr.columns if c.startswith("member_")]
            times   = tr["time"].astype(str).tolist()
            n_hours = len(tr)
        if sid not in net["seg_pos"]:
            continue
        f = net["factor"][sid]
        qlat_parts[sid] = tr[member_cols].values.astype("float32") * f   # (n_hours, n_members) m^3/s

    if not qlat_parts or member_cols is None:
        print(f"  [{tag}] no trajectory files found — skip")
        return

    n_members = len(member_cols)
    qlat = np.zeros((n_members, n_segs, n_hours), dtype="float32")
    for sid, arr in qlat_parts.items():
        qlat[:, net["seg_pos"][sid], :] = arr.T    # (n_members, n_hours)

    G["net"]  = net
    G["qlat"] = qlat

    print(f"  [{tag}] routing {n_members} members x {n_hours} hours, {procs} procs...")
    t0 = time.time()
    with Pool(procs) as pool:
        results = pool.map(_route_one_member, range(n_members))

    out = np.vstack(results).T                     # (n_hours, n_members) m^3/s
    df  = pd.DataFrame(out, columns=member_cols)
    df.insert(0, "time", times)
    df.to_parquet(out_path, index=False)
    print(f"  [{tag}] saved {out_path.name} ({len(df)} h x {n_members} members) "
          f"in {time.time()-t0:.0f}s")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Route full DA + OL analysis trajectories to gauge for all 12 AL gauges.")
    parser.add_argument("--da-traj-dir", required=True)
    parser.add_argument("--ol-traj-dir", required=True)
    parser.add_argument("--fc-dir",      required=True,
                        help="forecast dir; output goes to <fc-dir>/routed/<gauge>/")
    parser.add_argument("--gpkg-dir",    default="/mnt/disk2/al_storm_2026/gpkg")
    parser.add_argument("--procs",       type=int, default=10)
    parser.add_argument("--gauge",       default=None, help="single gauge (default: all 12)")
    parser.add_argument("--scenario",    default="both", choices=["da", "openloop", "both"])
    args = parser.parse_args()

    gauges = [args.gauge] if args.gauge else AL_GAUGES

    for gauge in gauges:
        gpkg_matches = sorted(Path(args.gpkg_dir).glob(f"*{gauge}*.gpkg"))
        if not gpkg_matches:
            print(f"\n[{gauge}] SKIP — no GPKG in {args.gpkg_dir}")
            continue
        gpkg = str(gpkg_matches[0])
        print(f"\n{'='*60}\nGauge {gauge}  |  GPKG: {Path(gpkg).name}")

        net = build_network(gpkg)
        print(f"  {len(net['cats'])} catchments | {net['area_km2']:.1f} km2 | "
              f"outlet cat-{net['terminal_int']}")

        out_dir = Path(args.fc_dir) / "routed" / gauge
        out_dir.mkdir(parents=True, exist_ok=True)

        if args.scenario in ("da", "both"):
            route_scenario("DA", args.da_traj_dir, net,
                           out_dir / "routed_analysis_da.parquet", args.procs)
        if args.scenario in ("openloop", "both"):
            route_scenario("OL", args.ol_traj_dir, net,
                           out_dir / "routed_analysis_ol.parquet", args.procs)

    print("\nDone. Re-run plot_forecast_all_gauges.py to include Oct 1–Oct 9 analysis.")


if __name__ == "__main__":
    main()
