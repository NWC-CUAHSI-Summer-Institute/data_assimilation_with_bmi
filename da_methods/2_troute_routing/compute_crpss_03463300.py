#!/usr/bin/env python3
"""
compute_crpss_03463300.py

Gauge-level CRPSS for South Toe River (03463300), held-out and held-in.

Reads routed parquets (m^3/s at gauge) produced by route_forecast.py:
    <fc_dir>/routed_5min/routed_forecast_leadtime_da.parquet
    <fc_dir>/routed_5min/routed_forecast_leadtime_openloop.parquet
  columns: issue_time, lead_hour, [valid_time], member_0000..member_NNNN

USGS obs (Sep 20-29 2024): tries these in order:
  1. /mnt/disk2/1400_sites_helene/nwm_benchmark/usgs_03463300_helene_m3s.csv
  2. /mnt/disk2/1400_sites_helene/heldout_kv_03463300/ (KV files)
  3. NWIS download (requires internet; saves a local copy)

crps_ensemble() is the exact O(m log m) NRG estimator from
  run_cfe_enkf_forecast_hrrr_multihour.py (Suma's commit).

Usage:
  ~/miniconda3/envs/troute/bin/python3 ~/compute_crpss_03463300.py
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAUGE       = "03463300"
FC_WINDOW   = ("2024-09-20", "2024-09-29")

# ── Routed parquet dirs ───────────────────────────────────────────────────────
HELD_OUT_FC_DIR = Path("/mnt/disk2/1400_sites_helene/da_forecast_heldout_no03463300_leadindexed/routed_5min")
HELD_IN_FC_DIR  = Path("/mnt/disk2/suma_helen_poster/da_forecast_heldin_range100_leadindexed/routed_5min")

OUT_DIR = Path("/home/svyas/qkrig_DA/results/crpss_03463300")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── USGS obs loading ──────────────────────────────────────────────────────────
USGS_CSV_SEARCH = [
    Path(f"/mnt/disk2/1400_sites_helene/nwm_benchmark/usgs_{GAUGE}_helene_m3s.csv"),
    Path(f"/mnt/disk2/1400_sites_helene/nwm_benchmark/usgs_{GAUGE}_hourly_full.csv"),
    Path(f"/mnt/disk2/1400_sites_helene/nwm_benchmark/usgs_{GAUGE}_hourly.csv"),
    Path(f"/mnt/disk2/suma_helen_poster/{GAUGE}_usgs_hourly_2018_2024.csv"),
    Path(f"/home/svyas/{GAUGE}_usgs_hourly.csv"),
]

# KV directory that includes 03463300 obs (excluded from main training KV)
HELDOUT_KV_DIR = Path("/mnt/disk2/1400_sites_helene/heldout_kv_03463300")


def load_usgs_csv(path, area_m2):
    df = pd.read_csv(path)
    date_col = next(c for c in df.columns
                    if c.lower() in ("date","time","datetime","timestamp","time_utc"))
    q_col    = next(c for c in df.columns
                    if any(k in c.lower() for k in ("m3","cms","discharge","flow","q_m3s","q")))
    df[date_col] = pd.to_datetime(df[date_col])
    q = df.set_index(date_col)[q_col].astype(float)
    if "mm" in q_col.lower():
        q = q * area_m2 / 1000.0 / 3600.0
        print(f"  Converted {q_col} mm/h -> m^3/s")
    else:
        print(f"  Using '{q_col}' as m^3/s")
    freq = pd.infer_freq(q.index[:20])
    if freq is not None and freq not in ("H", "h", "60min"):
        q = q.resample("h").mean()
        print(f"  Resampled from {freq} to hourly ({len(q)} hourly values)")
    return q


def load_usgs_from_kv(kv_dir, start, end):
    """Build hourly USGS obs for 03463300 from KV txt files (heldout_kv_03463300)."""
    records = {}
    for kv in sorted(kv_dir.glob("*.kv.txt")):
        stem = kv.stem.replace(".kv", "")
        try:
            ts = pd.Timestamp(stem.replace("_", " ") + ":00:00")
        except Exception:
            continue
        if ts < pd.Timestamp(start) or ts > pd.Timestamp(end):
            continue
        with open(kv) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                try:
                    gid, rest = line.split("=", 1)
                    if gid != GAUGE:
                        continue
                    parts = rest.split(",")
                    if parts[0] != "OK":
                        continue
                    records[ts] = float(parts[3])
                except Exception:
                    continue
    if not records:
        return None
    q = pd.Series(records).sort_index()
    print(f"  Loaded {len(q)} hourly USGS obs from KV files")
    return q


def download_usgs_nwis(gauge_id, start="2024-09-01", end="2024-09-30"):
    """Download hourly streamflow from USGS NWIS, return Series (m^3/s)."""
    import urllib.request
    url = (f"https://waterservices.usgs.gov/nwis/iv/?sites={gauge_id}"
           f"&parameterCd=00060&startDT={start}&endDT={end}"
           f"&siteStatus=all&format=rdb")
    print(f"  Downloading USGS NWIS: {url}")
    with urllib.request.urlopen(url, timeout=30) as r:
        raw = r.read().decode()
    lines = [l for l in raw.splitlines() if not l.startswith("#") and l.strip()]
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("agency_cd"))
    rows = [l.split("\t") for l in lines[header_idx+2:] if l.strip()]
    cols = lines[header_idx].split("\t")
    df = pd.DataFrame(rows, columns=cols)
    q_col = next(c for c in df.columns if "00060" in c and not c.endswith("cd"))
    df["datetime"] = pd.to_datetime(df["datetime"])
    q = df.set_index("datetime")[q_col].apply(pd.to_numeric, errors="coerce")
    q = q * 0.0283168  # cfs -> m^3/s
    q = q.resample("h").mean()
    out_csv = OUT_DIR / f"{gauge_id}_usgs_helene_window.csv"
    q.reset_index().to_csv(out_csv, index=False)
    print(f"  Saved USGS download: {out_csv}")
    return q


def get_usgs_obs(area_m2):
    for p in USGS_CSV_SEARCH:
        if p.exists():
            print(f"  Found USGS CSV: {p}")
            return load_usgs_csv(p, area_m2)

    # try KV files (Sep window)
    if HELDOUT_KV_DIR.exists():
        print(f"  Trying KV dir: {HELDOUT_KV_DIR}")
        q = load_usgs_from_kv(HELDOUT_KV_DIR, "2024-09-01", "2024-09-30")
        if q is not None and len(q) > 0:
            return q

    print("  No USGS CSV or KV found; downloading from NWIS...")
    try:
        return download_usgs_nwis(GAUGE, start="2024-09-01", end="2024-09-30")
    except Exception as e:
        print(f"  NWIS download failed: {e}")
        return None


# ── Basin area ────────────────────────────────────────────────────────────────
import sqlite3

def get_area_m2():
    candidates = sorted(Path("/mnt/disk2/1400_sites_helene/gpkg").glob(f"*{GAUGE}*.gpkg"))
    if not candidates:
        print("WARNING: no GPKG found; using hardcoded 92.5 km² (South Toe)")
        return 92.5e6
    con = sqlite3.connect(str(candidates[0]))
    div = pd.read_sql("SELECT areasqkm FROM divides", con)
    con.close()
    a = float(div["areasqkm"].sum()) * 1e6
    print(f"Basin area: {a/1e6:.2f} km²  ({candidates[0].name})")
    return a

area_m2 = get_area_m2()

print("\nLoading USGS obs...")
usgs_q = get_usgs_obs(area_m2)
if usgs_q is not None:
    print(f"  {len(usgs_q)} USGS obs, {usgs_q.index.min()} – {usgs_q.index.max()}")


# ── CRPS (exact NRG estimator from Suma's run_cfe_enkf_forecast_hrrr_multihour.py) ──
def crps_ensemble(obs, ens):
    """obs: (n,)  ens: (n, m) — both must be finite (filter before calling).
    Returns (n,) CRPS values in same units as obs."""
    n, m = ens.shape
    obs = np.asarray(obs, dtype=float).reshape(n, 1)
    term1 = np.abs(ens - obs).mean(axis=1)
    ens_sorted = np.sort(ens, axis=1)
    i = np.arange(1, m + 1)
    term2 = (2 * i - m - 1) @ ens_sorted.T / (m ** 2)
    return term1 - term2

def crpss(crps_fc, crps_ref):
    if crps_ref is None or not np.isfinite(crps_ref) or crps_ref == 0.0:
        return float("nan")
    return 1.0 - crps_fc / crps_ref

def compute_crps_by_lead(df, usgs_q, label):
    """Returns DataFrame: lead_hour, crps, n."""
    member_cols = [c for c in df.columns if c.startswith("member_")]
    if not member_cols:
        print(f"  ERROR: no member_ columns in {label}"); return None
    df = df.copy()
    df["issue_time"] = pd.to_datetime(df["issue_time"])
    if "valid_time" not in df.columns:
        df["valid_time"] = df["issue_time"] + pd.to_timedelta(df["lead_hour"], unit="h")
    else:
        df["valid_time"] = pd.to_datetime(df["valid_time"])
    df["obs"] = usgs_q.reindex(df["valid_time"]).values

    print(f"  {label}: {len(member_cols)} members, "
          f"{df['lead_hour'].nunique()} leads, "
          f"{df['issue_time'].nunique()} issues, "
          f"{df['obs'].notna().sum()}/{len(df)} obs matches")

    rows = []
    for lead, g in df.groupby("lead_hour"):
        o = g["obs"].values.astype(float)
        ens = g[member_cols].values.astype(float)
        finite = np.isfinite(o) & np.all(np.isfinite(ens), axis=1)
        if not finite.any():
            rows.append({"lead_hour": int(lead), "crps": np.nan, "n": 0})
            continue
        crps_vals = crps_ensemble(o[finite], ens[finite])
        rows.append({"lead_hour": int(lead), "crps": float(np.mean(crps_vals)), "n": int(finite.sum())})
    return pd.DataFrame(rows).sort_values("lead_hour").reset_index(drop=True)


# ── Process both experiments ──────────────────────────────────────────────────
experiments = {"held_out": HELD_OUT_FC_DIR, "held_in": HELD_IN_FC_DIR}
results = {}

for exp_name, fc_dir in experiments.items():
    print(f"\n{'='*60}\nExperiment: {exp_name}  →  {fc_dir}\n{'='*60}")
    if not fc_dir.exists():
        print(f"  SKIP — dir not found"); continue

    da_pq  = fc_dir / "routed_forecast_leadtime_da.parquet"
    ol_pq  = fc_dir / "routed_forecast_leadtime_openloop.parquet"

    # Fallback: first matching parquet
    if not da_pq.exists():
        m = sorted(fc_dir.glob("*da*.parquet"))
        da_pq = m[0] if m else None
    if not ol_pq.exists():
        m = sorted(fc_dir.glob("*openloop*.parquet")) + sorted(fc_dir.glob("*ol*.parquet"))
        ol_pq = m[0] if m else None

    if not da_pq or not da_pq.exists():
        print("  SKIP — DA parquet not found"); continue
    if not ol_pq or not ol_pq.exists():
        print("  SKIP — OL parquet not found"); continue
    if usgs_q is None:
        print("  SKIP — no USGS obs"); continue

    print(f"  DA : {da_pq.name}")
    print(f"  OL : {ol_pq.name}")

    df_da = pd.read_parquet(da_pq)
    df_ol = pd.read_parquet(ol_pq)

    crps_da = compute_crps_by_lead(df_da, usgs_q, f"DA  ({exp_name})")
    crps_ol = compute_crps_by_lead(df_ol, usgs_q, f"OL  ({exp_name})")
    if crps_da is None or crps_ol is None: continue

    merged = crps_da.merge(crps_ol, on="lead_hour", suffixes=("_da","_ol"))
    merged["crpss"] = merged.apply(
        lambda r: crpss(r["crps_da"], r["crps_ol"]), axis=1)

    results[exp_name] = {"crps_da": crps_da, "crps_ol": crps_ol, "merged": merged}

    print(f"\n  Lead | CRPS_DA  | CRPS_OL  | CRPSS   | n")
    print(f"  -----|----------|----------|---------|----")
    for _, r in merged.iterrows():
        print(f"   {int(r['lead_hour']):2d}  | {r['crps_da']:8.4f} | "
              f"{r['crps_ol']:8.4f} | {r['crpss']:+7.4f} | {int(r['n_da'])}")

    out_csv = OUT_DIR / f"crpss_03463300_{exp_name}.csv"
    merged.to_csv(out_csv, index=False)
    print(f"\n  Saved: {out_csv}")

if not results:
    # ── Fallback: aggregate per-catchment CRPSS from JSON summaries ──────────
    print("\n\nRouted parquets not found — falling back to catchment-level CRPSS from JSON summaries")
    for exp_name, fc_dir in experiments.items():
        base = fc_dir.parent
        json_files = sorted(base.glob("cat-*/cat-*_forecast_summary.json"))
        if not json_files:
            print(f"  {exp_name}: no JSON summaries found in {base}"); continue
        print(f"  {exp_name}: {len(json_files)} catchment JSON summaries")
        by_lead = {}
        for jf in json_files:
            with open(jf) as f:
                s = json.load(f)
            for lead_str, v in s.get("skill_by_lead_ens_mean_da", {}).items():
                lead = int(lead_str)
                by_lead.setdefault(lead, []).append(v.get("crpss_vs_openloop", np.nan))
        if not by_lead: continue
        rows = [{"lead_hour": lead, "crpss_mean": np.nanmean(vals),
                 "crpss_median": np.nanmedian(vals), "n_cats": len(vals)}
                for lead, vals in sorted(by_lead.items())]
        df = pd.DataFrame(rows)
        results[f"{exp_name}_catchment"] = df
        print(f"\n  Lead | CRPSS_mean | n_cats")
        for _, r in df.iterrows():
            print(f"   {int(r['lead_hour']):2d}  | {r['crpss_mean']:+.4f}     | {int(r['n_cats'])}")
        df.to_csv(OUT_DIR / f"crpss_03463300_{exp_name}_catchment_level.csv", index=False)

if not results:
    print("\nNo results to plot."); raise SystemExit(0)

# ── Combined plot ─────────────────────────────────────────────────────────────
colors = {"held_out": "tomato", "held_in": "steelblue"}
labels = {"held_out": "Held-out (03463300 excluded from krig)",
          "held_in":  "Held-in  (03463300 included in krig)"}

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: CRPSS
ax = axes[0]
for exp_name, res in results.items():
    if "merged" in res:
        df = res["merged"]
        ax.plot(df["lead_hour"], df["crpss"], "o-",
                color=colors.get(exp_name, "grey"), lw=2, ms=6,
                label=labels.get(exp_name, exp_name))
    else:
        ax.plot(df["lead_hour"], df["crpss_mean"], "o--",
                color=colors.get(exp_name.replace("_catchment",""), "grey"),
                lw=1.5, ms=5, alpha=0.7, label=f"{exp_name} (catchment avg)")
ax.axhline(0, color="grey", lw=0.8, ls="--")
ax.set_xlabel("Lead hour"); ax.set_ylabel("CRPSS (vs open loop)")
ax.set_title(f"Gauge {GAUGE} — CRPSS by lead hour")
ax.set_xticks(range(1, 19)); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_ylim(bottom=-0.1)

# Panel 2: Absolute CRPS
ax = axes[1]
for exp_name, res in results.items():
    if "crps_da" not in res: continue
    c = colors.get(exp_name, "grey")
    ax.plot(res["crps_da"]["lead_hour"], res["crps_da"]["crps"],
            "o-", color=c, lw=2, ms=5, label=f"DA ({exp_name.replace('_',' ')})")
    ax.plot(res["crps_ol"]["lead_hour"], res["crps_ol"]["crps"],
            "s--", color=c, lw=1.5, ms=4, alpha=0.55,
            label=f"OL ({exp_name.replace('_',' ')})")
ax.set_xlabel("Lead hour"); ax.set_ylabel("Mean CRPS (m³/s)")
ax.set_title(f"Gauge {GAUGE} — Mean CRPS by lead hour")
ax.set_xticks(range(1, 19)); ax.legend(fontsize=7); ax.grid(alpha=0.3)

plt.suptitle(f"CRPSS — South Toe River ({GAUGE})  |  Helene Sep 20-29 2024",
             fontsize=12, fontweight="bold")
plt.tight_layout()
out_png = OUT_DIR / "crpss_03463300_held_in_vs_held_out.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()
print(f"\nSaved: {out_png}")

# Individual bar charts
for exp_name, res in results.items():
    if "merged" not in res: continue
    df = res["merged"]
    fig, ax = plt.subplots(figsize=(10, 5))
    bar_colors = [colors.get(exp_name, "steelblue") if v >= 0 else "#aaaaaa"
                  for v in df["crpss"]]
    ax.bar(df["lead_hour"], df["crpss"], color=bar_colors, alpha=0.85, edgecolor="white")
    ax.axhline(0, color="black", lw=0.8)
    ax.fill_between(df["lead_hour"], 0, df["crpss"].clip(lower=0),
                    color=colors.get(exp_name,"steelblue"), alpha=0.2)
    ax.set_xlabel("Lead hour"); ax.set_ylabel("CRPSS (vs open loop)")
    ax.set_xticks(range(1, 19)); ax.grid(axis="y", alpha=0.3)
    ax.set_title(f"Gauge {GAUGE} — CRPSS vs open loop  [{labels.get(exp_name,exp_name)}]\n"
                 f"n_issues={int(df['n_da'].max())}, routed ensemble vs USGS")
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"crpss_03463300_{exp_name}_bar.png",
                dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved: {OUT_DIR}/crpss_03463300_{exp_name}_bar.png")

print("\nAll done.")
