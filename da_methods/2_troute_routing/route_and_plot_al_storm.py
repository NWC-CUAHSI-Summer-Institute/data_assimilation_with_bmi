"""Route calibrated CFE output and plot vs USGS for Alabama storm gauges.

Routing:  /mnt/disk2/al_storm_2026/gpkg/gage-{gauge}_subset.gpkg
Cal:      /mnt/disk2/svyas_clean_run/calibration_al_storm_2026/{cat}/{cat}_cal_results.csv
USGS obs: /mnt/disk2/al_storm_2026/usgs_obs_12/{gauge}.csv  (q_m3s column)
Output:   results/Alabama_gauges/{gauge}/gauge_{gauge}_routed.png

Run on dualearth1:
  python3 route_and_plot_al_storm.py [--gauge XXXXXXXX]
"""
import argparse, sqlite3, sys, types as _types
from functools import partial
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── T-Route stub (mirrors route_and_plot_heldin_v2.py) ──────────────────────
_stub = _types.ModuleType("coverage.types")
for _cls in ["Tracer","TTraceData","TShouldTraceFn","TFileDisposition",
             "TShouldStartContextFn","TWarnFn","TTraceFn"]:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules["coverage.types"] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

# ── Paths ────────────────────────────────────────────────────────────────────
GPKG_DIR  = Path("/mnt/disk2/al_storm_2026/gpkg")
CAL_DIR   = Path("/mnt/disk2/svyas_clean_run/calibration_al_storm_2026")
USGS_DIR  = Path("/mnt/disk2/al_storm_2026/usgs_obs_12")
OUT_BASE  = Path("/home/svyas/qkrig_DA/results/Alabama_gauges")

# ── Constants ────────────────────────────────────────────────────────────────
DT = 3600.0; SUBSTEPS = 12; QTS_SUBDIVISIONS = 1
CAL_START = "2020-01-01 00:00:00"
CAL_END   = "2022-12-31 23:00:00"

AL_GAUGES = [
    "02361000","02369800","02371500","02372250","02374500","02408540",
    "02422500","02450250","02464000","02465493","02469800","03574500",
]

# ── Routing helpers (identical to route_and_plot_heldin_v2.py) ───────────────
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
    except:
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
        "area_m2": area_km2 * 1e6,
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
        [], e00f64, {}, e00i32, False, "2020-01-01_00:00:00",
        e2f, e1i, e1i, e1i, e1f, e1f, 0.0, e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f, e2f, e1i, e1i, [], e1i, e1i, e1f,
        e1i, e1i, e1i, e1i, e1f, e1i, e1f, e1i, e1i, e00f32,
    )
    return np.asarray(res[1])[net["terminal_pos"], 0::3]

def route_ts(net, runoff_dict, nts):
    qlat = np.zeros((net["n_segs"], nts), dtype="float32")
    for cat, arr in runoff_dict.items():
        sid = int(cat[4:])
        if sid not in net["seg_pos"]: continue
        qlat[net["seg_pos"][sid], :] = np.asarray(arr, dtype="float32") * net["factor"][sid]
    qf = np.repeat(qlat, SUBSTEPS, axis=1)
    Q  = route_flow(net, qf, nts * SUBSTEPS, DT / SUBSTEPS)
    return Q[SUBSTEPS - 1::SUBSTEPS][:nts]

def nse(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim)) & (obs >= 0)
    if mask.sum() < 2: return np.nan
    o, s = obs[mask], sim[mask]
    return 1 - np.sum((o - s) ** 2) / np.sum((o - np.mean(o)) ** 2)

def kge_score(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim)) & (obs > 0)
    if mask.sum() < 2: return np.nan
    o, s = obs[mask], sim[mask]
    r = np.corrcoef(o, s)[0, 1]
    return 1 - np.sqrt((r - 1) ** 2 + (s.std() / o.std() - 1) ** 2 + (s.mean() / o.mean() - 1) ** 2)

# ── USGS obs (local CSV) ─────────────────────────────────────────────────────
def load_usgs_obs(gauge):
    """Read USGS hourly obs from local CSV; returns a Series of q_m3s indexed by time."""
    csv = USGS_DIR / f"{gauge}.csv"
    if not csv.exists():
        print(f"  {gauge}: no USGS CSV at {csv}")
        return None
    df = pd.read_csv(csv, parse_dates=["time"]).set_index("time")
    df.index = df.index.tz_localize(None)   # strip tz if present
    series = df["q_m3s"].astype("float64")
    series[series < 0] = np.nan
    return series

# ── Main ─────────────────────────────────────────────────────────────────────
def run_gauge(gauge, times):
    nts = len(times)
    gpkg_matches = sorted(GPKG_DIR.glob(f"*{gauge}*.gpkg"))
    if not gpkg_matches:
        print(f"{gauge}  SKIP — no GPKG in {GPKG_DIR}"); return None

    print(f"\n── {gauge} ── building network...")
    try:
        net = build_network(gpkg_matches[0])
    except Exception as e:
        print(f"{gauge}  SKIP — build_network: {e}"); return None

    cats = net["cats"]
    print(f"  {len(cats)} catchments, area={net['area_km2']:.1f} km²")

    sim_dict = {}
    missing  = 0
    for cat in cats:
        f = CAL_DIR / cat / f"{cat}_cal_results.csv"
        if f.exists():
            try:
                df = pd.read_csv(f, parse_dates=["date"]).set_index("date").reindex(times)
                sim_dict[cat] = df["sim_mm_h"].fillna(0).values
            except Exception:
                missing += 1
        else:
            missing += 1
    print(f"  sim loaded: {len(sim_dict)}/{len(cats)}  (missing: {missing})")

    if not sim_dict:
        print(f"{gauge}  SKIP — no sim data"); return None

    print("  routing CFE...")
    q_sim = route_ts(net, sim_dict, nts)

    print("  loading USGS obs...")
    usgs_series = load_usgs_obs(gauge)
    if usgs_series is None:
        print(f"{gauge}  SKIP — no USGS obs"); return None
    q_usgs = usgs_series.reindex(times).values

    kge = kge_score(q_usgs, q_sim)
    nse_val = nse(q_usgs, q_sim)
    print(f"  KGE_USGS={kge:.3f}  NSE_USGS={nse_val:.3f}")

    # Plot
    fig, ax = plt.subplots(figsize=(18, 5))
    ax.plot(times, q_usgs, color="black", linestyle=":", linewidth=1.0,
            alpha=0.85, label="USGS obs", zorder=3)
    ax.plot(times, q_sim,  color="#FFD700", linewidth=0.9, alpha=0.9,
            label="CFE + Troute", zorder=2)
    ax.set_ylabel("Streamflow (m³/s)")
    ax.set_title(
        f"Gauge {gauge}  |  Alabama storm 2026  |  Cal 2020–2022  |  "
        f"KGE={kge:.3f}  NSE={nse_val:.3f}"
    )
    ax.legend(fontsize=9)
    valid_usgs = q_usgs[~np.isnan(q_usgs)]
    y_max = max(float(np.nanpercentile(valid_usgs, 99.5)) * 1.3, 5.0) if len(valid_usgs) else 20.0
    ax.set_ylim(0, y_max)
    ax.set_xlim(times[0], times[-1])
    plt.tight_layout()

    out_dir = OUT_BASE / gauge
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"gauge_{gauge}_routed.png"
    plt.savefig(str(out_png), dpi=130, bbox_inches="tight")
    plt.close()
    print(f"  saved → {out_png}")
    return {"gauge": gauge, "KGE_USGS": kge, "NSE_USGS": nse_val}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gauge", default=None,
                        help="Single gauge ID to process (omit for all 12)")
    args = parser.parse_args()

    gauges = [args.gauge] if args.gauge else AL_GAUGES
    times  = pd.date_range(CAL_START, CAL_END, freq="h")

    results = []
    for g in gauges:
        res = run_gauge(g, times)
        if res: results.append(res)

    if results:
        df = pd.DataFrame(results)
        print("\n── Summary ──")
        print(df.to_string(index=False))
        print(f"Median KGE_USGS: {df['KGE_USGS'].median():.3f}")
        print(f"Median NSE_USGS: {df['NSE_USGS'].median():.3f}")
    print("\nAll done.")
