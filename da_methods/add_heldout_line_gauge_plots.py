#!/usr/bin/env python3
"""
add_heldout_line_gauge_plots.py

Gauge-level calibration plots with 5 lines:
  - USGS obs              black dots
  - Krig obs + T-Route    blue
  - Held-in CFE + T-Route yellow
  - Held-out CFE + T-Route red
  - Lumped CFE            green

Usage (dualearth1, troute env):
  ~/miniconda3/envs/troute/bin/python3 ~/add_heldout_line_gauge_plots.py --all
  ~/miniconda3/envs/troute/bin/python3 ~/add_heldout_line_gauge_plots.py --gauge 02018000
"""

import argparse
import sqlite3
import sys
import types as _types
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_stub = _types.ModuleType("coverage.types")
for _cls in ["Tracer","TTraceData","TShouldTraceFn","TFileDisposition",
             "TShouldStartContextFn","TWarnFn","TTraceFn"]:
    setattr(_stub, _cls, type(_cls, (), {}))
sys.modules["coverage.types"] = _stub

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import compute_network_structured

# ── Paths ─────────────────────────────────────────────────────────────────────
GPKG_DIR        = Path("/mnt/disk2/1400_sites_helene/gpkg")
HELDIN_V2_DIR   = Path("/mnt/disk2/svyas_clean_run/calibration_heldin_v2")
HELDOUT_V2_DIR  = Path("/mnt/disk2/svyas_clean_run/calibration_heldout_v2")
KRIG_OBS_DIR    = Path("/mnt/disk2/1400_sites_helene/krig_obs_fold0_heldout_v2")
GAUGE_CAL_DIR   = Path("/mnt/disk2/svyas_clean_run")
OUT_DIR         = Path("/home/svyas/qkrig_DA/results/Gauge Level calibration")

# Special gauges: override cal and krig obs paths
SPECIAL_HELDIN_DIR  = {
    "03456500": Path("/mnt/disk2/svyas_clean_run/calibration_03456500_heldin"),
    "03463300": Path("/mnt/disk2/suma_helen_poster/catchment_results_range100"),
}
SPECIAL_HELDOUT_DIR = {
    "03456500": Path("/mnt/disk2/svyas_clean_run/calibration_03456500"),
    "03463300": Path("/mnt/disk2/suma_helen_poster/catchment_results_1gauge_heldout"),
}
SPECIAL_KRIG_DIR = {
    "03456500": Path("/mnt/disk2/1400_sites_helene/catchment_ts_03456500_heldin_dynamic_variance"),
    "03463300": Path("/mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig"),
}
SPECIAL_KRIG_COL = {
    "03456500": "qkrig_mm_hr",
    "03463300": "qkrig_mm_hr",
}

# Area override for gauges where special dirs cover only a subset of catchments.
# 03456500: only 20 East Fork cats (cat-1016565..1016584) have cal results;
# the full GPKG may include the entire upstream network at higher area.
SPECIAL_AREA_KM2 = {
    "03456500": 133.38,   # USGS drain_area_va; matches the KV obs mm/h conversion
}

# Terminal segment override: wb-1016569 is the true East Fork outlet.
# The 43-cat GPKG never links EF segs to the main-network terminal.
SPECIAL_TERMINAL_INT = {
    "03456500": 1016569,
}

DT               = 3600.0
SUBSTEPS         = 12
QTS_SUBDIVISIONS = 1
CAL_START        = "2020-01-01 00:00:00"
CAL_END          = "2022-12-31 23:00:00"

GAUGES = [
    "01606500","01667500","02013000","02014000","02018000",
    "02038850","02046000","02055100","02082950","02096846",
    "02149000","02193340","02198100","02221525","02408540",
    "02450250","02465493","03164000","03285000","03488000","03604000",
    "03456500","03463300",
]

# ── Routing ───────────────────────────────────────────────────────────────────
def seg_int(s): return int(s[3:])

def build_network(gpkg):
    con = sqlite3.connect(str(gpkg))
    div = pd.read_sql("SELECT divide_id, areasqkm FROM divides", con)
    cats = sorted(str(c) for c in div["divide_id"] if str(c).startswith("cat-"))
    area_km2 = float(div["areasqkm"].sum())
    try:
        nt = pd.read_sql("SELECT divide_id, tot_drainage_areasqkm FROM network WHERE divide_id IS NOT NULL", con)
        outlet = str(nt.loc[nt["tot_drainage_areasqkm"].idxmax(), "divide_id"])
    except:
        outlet = cats[-1]
    fp_attr = pd.read_sql('SELECT link,"to",BtmWdth,TopWdth,TopWdthCC,n,nCC,ChSlp,So,Length_m FROM "flowpath-attributes"', con)
    fp      = pd.read_sql("SELECT divide_id, areasqkm FROM flowpaths", con)
    con.close()
    terminal_int = int(outlet[4:])
    link_set = {seg_int(r) for r in fp_attr["link"]}
    connections = {seg_int(row["link"]): ([int(row["to"][4:])] if int(row["to"][4:]) in link_set else [])
                   for _, row in fp_attr.iterrows()}
    rconn = nhd_network.reverse_network(connections)
    reach_list = nhd_network.dfs_decomposition(rconn, partial(nhd_network.split_at_junction, rconn))
    param_df = pd.DataFrame([{
        "seg_id": seg_int(r["link"]), "dt": DT,
        "bw": float(r["BtmWdth"]), "tw": float(r["TopWdth"]), "twcc": float(r["TopWdthCC"]),
        "dx": float(r["Length_m"]), "n": float(r["n"]),  "ncc": float(r["nCC"]),
        "cs": float(r["ChSlp"]),   "s0": float(r["So"]), "alt": 0.0,
    } for _, r in fp_attr.iterrows()]).set_index("seg_id").sort_index().astype("float32")
    seg_ids  = param_df.index.values
    area_map = {int(r["divide_id"][4:]): r["areasqkm"]*1e6 for _, r in fp.iterrows()
                if r["divide_id"] and str(r["divide_id"]).startswith("cat-")}
    valid = [c for c in cats if int(c[4:]) in link_set]
    return {
        "reaches_wTypes": [(r, 0) for r in reach_list],
        "upstreams":  dict(rconn),
        "param_df":   param_df,
        "n_segs":     len(param_df),
        "cats":       valid,
        "area_km2":   area_km2,
        "area_m2":    area_km2 * 1e6,
        "terminal_pos": int(np.where(seg_ids == terminal_int)[0][0]),
        "seg_pos":  {int(c[4:]): int(np.where(seg_ids == int(c[4:]))[0][0]) for c in valid},
        "factor":   {int(c[4:]): area_map.get(int(c[4:]), 0) / 1000.0 / 3600.0 for c in valid},
        "q0_df":    pd.DataFrame(np.zeros((len(param_df), 3), dtype="float32"),
                                 index=param_df.index, columns=["qu0","qd0","h0"]),
    }

def route_flow(net, qlat, nts, dt):
    e1i = np.zeros(0, dtype="int32");  e1f  = np.zeros(0, dtype="float32")
    e2f = np.zeros((0, nts), dtype="float32"); e00f32 = np.zeros((0,0), dtype="float32")
    e00f64 = np.zeros((0,0), dtype="float64"); e00i32 = np.zeros((0,0), dtype="int32")
    res = compute_network_structured(
        nts, dt, QTS_SUBDIVISIONS, net["reaches_wTypes"], net["upstreams"],
        net["param_df"].index.values.astype("int64"), net["param_df"].columns.values,
        net["param_df"].values, net["q0_df"].values.astype("float32"), qlat.astype("float32"),
        [], e00f64, {}, e00i32, False, "2020-01-01_00:00:00",
        e2f, e1i, e1i, e1i, e1f, e1f, 0.0, e2f, e1i, e1f, e1f, e1f, e1f, e1f,
        e2f, e1i, e1f, e1f, e1f, e1f, e1f, e2f, e1i, e1i, [], e1i, e1i, e1f, e1i, e1i,
        e1i, e1i, e1f, e1i, e1f, e1i, e1i, e00f32)
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
    return 1 - np.sum((o - s)**2) / np.sum((o - np.mean(o))**2)

def kge_score(obs, sim):
    mask = ~(np.isnan(obs) | np.isnan(sim)) & (obs > 0)
    if mask.sum() < 2: return np.nan
    o, s = obs[mask], sim[mask]
    r = np.corrcoef(o, s)[0, 1]
    return 1 - np.sqrt((r-1)**2 + (s.std()/o.std()-1)**2 + (s.mean()/o.mean()-1)**2)

# ── Load helpers ──────────────────────────────────────────────────────────────
def load_cal_and_route(net, cal_dir, times, nts):
    sim_dict = {}
    for cat in net["cats"]:
        f = cal_dir / cat / f"{cat}_cal_results.csv"
        if f.exists():
            try:
                df = pd.read_csv(f, parse_dates=["date"]).set_index("date").reindex(times)
                sim_dict[cat] = df["sim_mm_h"].fillna(0).values
            except Exception:
                pass
    if not sim_dict:
        return None
    print(f"    cal loaded: {len(sim_dict)}/{len(net['cats'])} cats from {cal_dir.name}")
    q_m3s = route_ts(net, sim_dict, nts)
    return q_m3s * 1000.0 * 3600.0 / net["area_m2"]

def load_krig_and_route(net, krig_dir, col, times, nts):
    obs_dict = {}
    for cat in net["cats"]:
        f = krig_dir / f"{cat}.csv"
        if f.exists():
            try:
                df = pd.read_csv(f, parse_dates=["time"]).set_index("time").reindex(times)
                obs_dict[cat] = df[col].fillna(0).values
            except Exception:
                pass
    if not obs_dict:
        return None
    print(f"    krig loaded: {len(obs_dict)}/{len(net['cats'])} cats from {krig_dir.name}")
    q_m3s = route_ts(net, obs_dict, nts)
    return q_m3s * 1000.0 * 3600.0 / net["area_m2"]

# ── Per-gauge ─────────────────────────────────────────────────────────────────
def process_gauge(gauge, times, nts):
    print(f"\nGauge {gauge}")

    # 1. Gauge-level lumped cal CSV
    gauge_csv = GAUGE_CAL_DIR / f"calibration_gauge_{gauge}" / f"{gauge}_cal_results.csv"
    if not gauge_csv.exists():
        print(f"  SKIP — no gauge-level cal CSV"); return
    df_g = pd.read_csv(gauge_csv, parse_dates=["date"]).set_index("date").reindex(times)
    sim_lumped = df_g["sim_mm_h"].values
    obs_usgs   = df_g["obs_mm_h"].values

    # 2. Build routing network
    gpkg_matches = sorted(GPKG_DIR.glob(f"*{gauge}*.gpkg"))
    if not gpkg_matches:
        print(f"  SKIP — no GPKG"); return
    try:
        net = build_network(gpkg_matches[0])
    except Exception as e:
        print(f"  SKIP — network build failed: {e}"); return

    # 3. Route: held-in v2, held-out v2, krig obs
    heldin_dir  = SPECIAL_HELDIN_DIR.get(gauge,  HELDIN_V2_DIR)
    heldout_dir = SPECIAL_HELDOUT_DIR.get(gauge, HELDOUT_V2_DIR)
    krig_dir    = SPECIAL_KRIG_DIR.get(gauge,    KRIG_OBS_DIR)
    krig_col    = SPECIAL_KRIG_COL.get(gauge,    "qkrig_mm_hr")

    # Override area when only a subset of catchments have results
    if gauge in SPECIAL_AREA_KM2:
        old_area = net["area_km2"]
        net["area_km2"] = SPECIAL_AREA_KM2[gauge]
        net["area_m2"]  = SPECIAL_AREA_KM2[gauge] * 1e6
        print(f"  Area override: {old_area:.2f} km2 → {net['area_km2']:.2f} km2")

    if gauge in SPECIAL_TERMINAL_INT:
        t_int = SPECIAL_TERMINAL_INT[gauge]
        seg_ids = net["param_df"].index.values
        idx = int(np.where(seg_ids == t_int)[0][0])
        net["terminal_pos"] = idx
        print(f"  Terminal override: wb-{t_int}")

    q_heldin  = load_cal_and_route(net, heldin_dir,  times, nts)
    q_heldout = load_cal_and_route(net, heldout_dir, times, nts)
    q_krig    = load_krig_and_route(net, krig_dir, krig_col, times, nts)

    # 4. Metrics
    def m(arr):
        if arr is None: return float("nan"), float("nan")
        return nse(obs_usgs, arr), kge_score(obs_usgs, arr)

    nse_lu, kge_lu   = m(sim_lumped)
    nse_hi, kge_hi   = m(q_heldin)
    nse_ho, kge_ho   = m(q_heldout)
    nse_kr, kge_kr   = m(q_krig)

    print(f"  Lumped CFE:         NSE={nse_lu:.3f}  KGE={kge_lu:.3f}")
    print(f"  Held-in  v2+route:  NSE={nse_hi:.3f}  KGE={kge_hi:.3f}")
    print(f"  Held-out v2+route:  NSE={nse_ho:.3f}  KGE={kge_ho:.3f}")
    print(f"  Krig+route:         NSE={nse_kr:.3f}  KGE={kge_kr:.3f}")

    # 5. Plot — 2-panel layout
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), sharex=False)
    zoom = (times >= pd.Timestamp("2021-01-01")) & (times < pd.Timestamp("2022-01-01"))

    for ax, mask in [(ax1, np.ones(nts, dtype=bool)), (ax2, zoom)]:
        t = times[mask]
        ax.plot(t, obs_usgs[mask], 'k.', markersize=1.5, alpha=0.6, label="USGS obs")
        if q_krig is not None:
            ax.plot(t, q_krig[mask],   color="steelblue", lw=0.9, alpha=0.85,
                    label=f"Krig + T-Route  (NSE={nse_kr:.3f} | KGE={kge_kr:.3f})")
        if q_heldin is not None:
            ax.plot(t, q_heldin[mask], color="goldenrod",  lw=0.9, alpha=0.85,
                    label=f"Held-in CFE + T-Route  (NSE={nse_hi:.3f} | KGE={kge_hi:.3f})")
        if q_heldout is not None:
            ax.plot(t, q_heldout[mask], color="tomato",   lw=0.9, alpha=0.85,
                    label=f"Held-out CFE + T-Route  (NSE={nse_ho:.3f} | KGE={kge_ho:.3f})")
        ax.plot(t, sim_lumped[mask], color="seagreen", lw=0.9, alpha=0.85,
                label=f"Lumped CFE  (NSE={nse_lu:.3f} | KGE={kge_lu:.3f})")
        ax.set_ylabel("Discharge (mm/h)")
        ax.legend(fontsize=9)

    ax1.set_title(f"Gauge {gauge} | Cal 2020–2022")
    ax2.set_title(f"Gauge {gauge} | Zoom: 2021")
    plt.tight_layout()

    out_dir = OUT_DIR / gauge
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / f"{gauge}_cal_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gauge", default=None)
    parser.add_argument("--all",   action="store_true")
    args = parser.parse_args()
    if not args.gauge and not args.all:
        parser.error("Specify --gauge XXXXXXXX or --all")

    times = pd.date_range(CAL_START, CAL_END, freq="h")
    nts   = len(times)
    gauges = GAUGES if args.all else [args.gauge]
    for gauge in gauges:
        try:
            process_gauge(gauge, times, nts)
        except Exception as e:
            import traceback
            print(f"  FAILED {gauge}: {e}"); traceback.print_exc()
    print("\nDone.")


if __name__ == "__main__":
    main()
