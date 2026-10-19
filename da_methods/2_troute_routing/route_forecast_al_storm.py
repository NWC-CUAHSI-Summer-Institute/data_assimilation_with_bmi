#!/usr/bin/env python3
"""
route_forecast_al_storm.py -- Route 18h ensemble forecast through T-Route
for the 12 Alabama storm gauges.

Reads:   /mnt/disk2/svyas_clean_run/forecast_al_storm_2026/<cat>/<cat>_forecast18h.parquet
         (18 timesteps x 600 members, mm/h)
         /mnt/disk2/al_storm_2026/gpkg/gage-{gauge}_subset.gpkg

Writes:  /home/svyas/qkrig_DA/results/Alabama_gauges/{gauge}/gauge_{gauge}_forecast18h.parquet
         columns: valid_time, lead_h, member_0000 .. member_0NNN (m³/s)
         /home/svyas/qkrig_DA/results/Alabama_gauges/{gauge}/gauge_{gauge}_forecast18h.png
         ensemble fan plot (m³/s vs valid time)

USAGE
    python3 route_forecast_al_storm.py [--gauge XXXXXXXX]
"""

import argparse
import sqlite3
import sys
import types as _types
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── T-Route stub (same as route_and_plot_al_storm.py) ────────────────────────
_stub = _types.ModuleType("coverage.types")
for _cls in ["Tracer", "TTraceData", "TShouldTraceFn", "TFileDisposition",
             "TShouldStartContextFn", "TWarnFn", "TTraceFn"]:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules["coverage.types"] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

# ── Paths ─────────────────────────────────────────────────────────────────────
GPKG_DIR    = Path("/mnt/disk2/al_storm_2026/gpkg")
FC_DIR      = Path("/mnt/disk2/svyas_clean_run/forecast_al_storm_2026")
OUT_BASE    = Path("/home/svyas/qkrig_DA/results/Alabama_gauges")

# ── Constants ─────────────────────────────────────────────────────────────────
DT = 3600.0; SUBSTEPS = 12; QTS_SUBDIVISIONS = 1
LEAD_HOURS = 18

AL_GAUGES = [
    "02361000", "02369800", "02371500", "02372250", "02374500", "02408540",
    "02422500", "02450250", "02464000", "02465493", "02469800", "03574500",
]


# ── Network helpers (identical to route_and_plot_al_storm.py) ─────────────────
def seg_int(s): return int(s[3:])

def derive_basin(gpkg):
    con = sqlite3.connect(str(gpkg))
    div = pd.read_sql("SELECT divide_id, areasqkm FROM divides", con)
    cats = sorted(str(c) for c in div["divide_id"] if str(c).startswith("cat-"))
    area_km2 = float(div["areasqkm"].sum())
    try:
        nt = pd.read_sql(
            "SELECT divide_id, tot_drainage_areasqkm FROM network WHERE divide_id IS NOT NULL", con)
        outlet = str(nt.loc[nt["tot_drainage_areasqkm"].idxmax(), "divide_id"])
    except Exception:
        outlet = cats[-1]
    con.close()
    return cats, area_km2, outlet

def build_network(gpkg):
    cats, area_km2, outlet_cat = derive_basin(gpkg)
    terminal_int = int(outlet_cat[4:])
    con = sqlite3.connect(str(gpkg))
    fp_attr = pd.read_sql(
        'SELECT link,"to",BtmWdth,TopWdth,TopWdthCC,n,nCC,ChSlp,So,Length_m FROM "flowpath-attributes"', con)
    fp = pd.read_sql("SELECT divide_id, areasqkm FROM flowpaths", con)
    con.close()
    link_set = {seg_int(r) for r in fp_attr["link"]}
    connections = {
        seg_int(row["link"]): (
            [int(row["to"][4:])] if int(row["to"][4:]) in link_set else []
        )
        for _, row in fp_attr.iterrows()
    }
    rconn = nhd_network.reverse_network(connections)
    reach_list = nhd_network.dfs_decomposition(
        rconn, partial(nhd_network.split_at_junction, rconn))
    param_df = pd.DataFrame([{
        "seg_id": seg_int(r["link"]), "dt": DT,
        "bw": float(r["BtmWdth"]), "tw": float(r["TopWdth"]),
        "twcc": float(r["TopWdthCC"]), "dx": float(r["Length_m"]),
        "n": float(r["n"]), "ncc": float(r["nCC"]),
        "cs": float(r["ChSlp"]), "s0": float(r["So"]), "alt": 0.0,
    } for _, r in fp_attr.iterrows()]).set_index("seg_id").sort_index().astype("float32")
    seg_ids = param_df.index.values
    area_map = {
        int(r["divide_id"][4:]): r["areasqkm"] * 1e6
        for _, r in fp.iterrows()
        if r["divide_id"] and str(r["divide_id"]).startswith("cat-")
    }
    valid = [c for c in cats if int(c[4:]) in link_set]
    return {
        "reaches_wTypes": [(r, 0) for r in reach_list],
        "upstreams": dict(rconn),
        "param_df": param_df,
        "n_segs": len(param_df),
        "cats": valid,
        "area_km2": area_km2,
        "terminal_int": terminal_int,
        "terminal_pos": int(np.where(seg_ids == terminal_int)[0][0]),
        "seg_pos": {
            int(c[4:]): int(np.where(seg_ids == int(c[4:]))[0][0]) for c in valid
        },
        "factor": {int(c[4:]): area_map.get(int(c[4:]), 0) / 1000.0 / 3600.0 for c in valid},
        "q0_df": pd.DataFrame(
            np.zeros((len(param_df), 3), dtype="float32"),
            index=param_df.index, columns=["qu0", "qd0", "h0"],
        ),
    }

def route_flow(net, qlat, nts, dt):
    e1i = np.zeros(0, dtype="int32"); e1f = np.zeros(0, dtype="float32")
    e2f = np.zeros((0, nts), dtype="float32"); e00f32 = np.zeros((0, 0), dtype="float32")
    e00f64 = np.zeros((0, 0), dtype="float64"); e00i32 = np.zeros((0, 0), dtype="int32")
    res = compute_network_structured(
        nts, dt, QTS_SUBDIVISIONS, net["reaches_wTypes"], net["upstreams"],
        net["param_df"].index.values.astype("int64"), net["param_df"].columns.values,
        net["param_df"].values, net["q0_df"].values.astype("float32"),
        qlat.astype("float32"),
        [], e00f64, {}, e00i32, False, "2026-10-08_00:00:00",
        e2f, e1i, e1i, e1i, e1f, e1f, 0.0, e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f, e2f, e1i, e1i, [], e1i, e1i, e1f,
        e1i, e1i, e1i, e1i, e1f, e1i, e1f, e1i, e1i, e00f32,
    )
    return np.asarray(res[1])[net["terminal_pos"], 0::3]

def route_one_member(net, runoff_dict, nts):
    """Route one member's runoff (mm/h) through T-Route, return outlet q (m³/s)."""
    qlat = np.zeros((net["n_segs"], nts), dtype="float32")
    for cat, arr in runoff_dict.items():
        sid = int(cat[4:])
        if sid not in net["seg_pos"]: continue
        qlat[net["seg_pos"][sid], :] = np.asarray(arr, dtype="float32") * net["factor"][sid]
    qf = np.repeat(qlat, SUBSTEPS, axis=1)
    Q  = route_flow(net, qf, nts * SUBSTEPS, DT / SUBSTEPS)
    return Q[SUBSTEPS - 1::SUBSTEPS][:nts]


# ── Load forecast parquets for one gauge's catchments ─────────────────────────
def load_gauge_forecasts(gauge, cats, fc_dir, fc_tag=""):
    """Return dict: cat_id -> array shape (LEAD_HOURS, n_members) in mm/h."""
    data = {}
    member_cols = None
    missing = 0
    tag_sfx = f"_{fc_tag}" if fc_tag else ""
    for cat in cats:
        pq = fc_dir / cat / f"{cat}_forecast18h{tag_sfx}.parquet"
        if not pq.exists():
            missing += 1
            continue
        df = pd.read_parquet(pq)
        if member_cols is None:
            member_cols = [c for c in df.columns if c.startswith("member_")]
        arr = df[member_cols].values.astype("float32")   # (18, n_members)
        data[cat] = arr
    if missing:
        print(f"  {gauge}: {missing}/{len(cats)} catchments missing forecast parquet")
    return data, member_cols


def route_gauge(gauge, args):
    """Route 600-member 18h forecast for one gauge, return summary dict."""
    gpkg_matches = sorted(GPKG_DIR.glob(f"*{gauge}*.gpkg"))
    if not gpkg_matches:
        print(f"{gauge}  SKIP — no GPKG in {GPKG_DIR}")
        return None

    print(f"\n── {gauge} ── building network...")
    try:
        net = build_network(gpkg_matches[0])
    except Exception as e:
        print(f"{gauge}  SKIP — build_network: {e}")
        return None

    cats = net["cats"]
    print(f"  {len(cats)} catchments | loading forecast parquets...")
    fc_data, member_cols = load_gauge_forecasts(gauge, cats, FC_DIR, getattr(args, "fc_tag", ""))
    if not fc_data:
        print(f"{gauge}  SKIP — no forecast data")
        return None

    n_members = len(member_cols)
    n_loaded  = len(fc_data)
    print(f"  {n_loaded}/{len(cats)} catchments with data | {n_members} members | routing...")

    # Route each member: shape (n_members, LEAD_HOURS)
    q_ensemble = np.zeros((n_members, LEAD_HOURS), dtype="float32")
    for mi in range(n_members):
        # build per-member runoff dict: cat_id -> (LEAD_HOURS,) array
        runoff_dict = {cat: arr[:, mi] for cat, arr in fc_data.items()}
        q_ensemble[mi, :] = route_one_member(net, runoff_dict, LEAD_HOURS)
        if (mi + 1) % 100 == 0:
            print(f"    routed {mi+1}/{n_members} members")

    # Get valid_times from one parquet
    sample_pq = next(FC_DIR / cat / f"{cat}_forecast18h.parquet"
                     for cat in cats if cat in fc_data)
    sample_df  = pd.read_parquet(sample_pq)
    valid_times = sample_df["valid_time"].tolist()
    lead_hs     = sample_df["lead_h"].tolist()

    # Save ensemble parquet: rows = 18 leads, cols = member_0000..
    out_dir = OUT_BASE / gauge
    out_dir.mkdir(parents=True, exist_ok=True)

    out_cols = [f"member_{i:04d}" for i in range(n_members)]
    fc_df = pd.DataFrame(q_ensemble.T, columns=out_cols, dtype=np.float32)
    fc_df.insert(0, "valid_time", valid_times)
    fc_df.insert(1, "lead_h", lead_hs)
    tag_sfx = f"_{args.out_tag}" if getattr(args, "out_tag", "") else ""
    fc_df.to_parquet(out_dir / f"gauge_{gauge}_forecast18h{tag_sfx}.parquet", index=False)

    # Mean / std per lead
    mean = q_ensemble.mean(axis=0)
    std  = q_ensemble.std(axis=0, ddof=1)
    peak_mean = float(mean.max())
    peak_lead = int(np.argmax(mean)) + 1

    # Plot
    dates = pd.to_datetime(valid_times)
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.fill_between(dates, np.clip(mean - 2*std, 0, None), mean + 2*std,
                    color="steelblue", alpha=0.25, label="±2 std")
    for pct, alpha in ((5, 0.1), (25, 0.15)):
        lo = np.percentile(q_ensemble, pct, axis=0)
        hi = np.percentile(q_ensemble, 100 - pct, axis=0)
        ax.fill_between(dates, lo, hi, color="steelblue", alpha=alpha)
    ax.plot(dates, mean, "steelblue", lw=2, label="ensemble mean")
    ax.plot(dates, np.percentile(q_ensemble, 50, axis=0),
            "navy", lw=1, ls="--", label="median")
    ax.set_ylabel("Streamflow (m³/s)")
    ax.set_title(f"Gauge {gauge} | 18h forecast | {n_members} members | "
                 f"peak={peak_mean:.2f} m³/s at lead {peak_lead}h")
    ax.legend(fontsize=9)
    ax.set_xlim(dates[0], dates[-1])
    plt.tight_layout()
    plt.savefig(out_dir / f"gauge_{gauge}_forecast18h{tag_sfx}.png", dpi=130, bbox_inches="tight")
    plt.close()
    print(f"  peak mean={peak_mean:.3f} m³/s at lead {peak_lead}h  |  saved -> {out_dir.name}/")

    return {"gauge": gauge, "n_cats_loaded": n_loaded, "n_members": n_members,
            "peak_mean_m3s": peak_mean, "peak_lead_h": peak_lead,
            "p10_peak_m3s": float(np.percentile(q_ensemble.max(axis=1), 10)),
            "p90_peak_m3s": float(np.percentile(q_ensemble.max(axis=1), 90))}


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gauge", default=None,
                        help="single gauge ID (omit for all 12)")
    parser.add_argument("--fc-dir",   default=str(FC_DIR))
    parser.add_argument("--out-base", default=str(OUT_BASE))
    parser.add_argument("--fc-tag",   default="",
                        help="suffix on input forecast parquets, e.g. '2026100912' reads _forecast18h_2026100912.parquet")
    parser.add_argument("--out-tag",  default="",
                        help="suffix on output parquet/png, same value as --fc-tag")
    args = parser.parse_args()

    FC_DIR   = Path(args.fc_dir)
    OUT_BASE = Path(args.out_base)
    gauges   = [args.gauge] if args.gauge else AL_GAUGES

    results = []
    for g in gauges:
        r = route_gauge(g, args)
        if r:
            results.append(r)

    if results:
        df = pd.DataFrame(results)
        out_csv = OUT_BASE / "forecast18h_routed_summary.csv"
        df.to_csv(out_csv, index=False)
        print(f"\n── Routed forecast summary ──")
        print(df[["gauge","peak_mean_m3s","peak_lead_h","p10_peak_m3s","p90_peak_m3s"]]
              .to_string(index=False))
        print(f"\nSaved: {out_csv}")
    print("\nDone.")
