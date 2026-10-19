#!/usr/bin/env python3
"""
run_forecast_al_storm.py -- 18-hour ensemble forecast for Alabama storm 2026.

Starts from the 600-member DA states saved by run_cfe_enkf_al_storm.py
(Oct 7-8 member_states_al_storm.parquet), drives each member forward 18 hours
with perturbed HRRR forecast forcing, and saves the ensemble streamflow.

WORKFLOW
    1. Open dynamical.org HRRR, select FORECAST_INIT_TIME + lead_times 1h..18h
       for the Alabama bbox, extract per-catchment centroids.
    2. Load member states from <da-dir>/<cat>/<cat>_member_states_al_storm.parquet
       at FORECAST_STATE_TIME (last hour of DA window = Oct 8 17:00).
    3. Reinitialize 600 CFE members from those states (soil, GW, Nash cascades).
    4. Step each member forward 18h using perturbed HRRR forecast precip/PET.
    5. Save ensemble runoff: <out-dir>/<cat>/<cat>_forecast18h.parquet
       columns: valid_time, lead_h, member_0000 .. member_0599  (mm/h, float32)

USAGE
    # smoke test (4 members, one catchment)
    python run_forecast_al_storm.py --cat-id cat-1011916 --n-members 4 \
        --out-dir /tmp/fc_smoke

    # full run (all 1072 catchments, 600 members)
    python run_forecast_al_storm.py --cat-id all

OUTPUTS
    <cat>_forecast18h.parquet   -- ensemble runoff (mm/h) for each of 18 leads
    <cat>_forecast18h_mean.csv  -- ensemble mean + std for quick inspection
    <cat>_forecast18h.png       -- plot: ensemble fan + mean
"""

import argparse
import json
import sqlite3
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xarray as xr

bmi_cfe = None

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DEFAULT_DA_DIR      = "/mnt/disk2/svyas_clean_run/da_al_storm_2026"
DEFAULT_PARAMS_DIR  = "/mnt/disk2/svyas_clean_run/calibration_al_storm_2026"
DEFAULT_OBS_DIR     = "/mnt/disk2/al_storm_2026/catchment_obs_r100"
DEFAULT_GPKG_DIR    = "/mnt/disk2/al_storm_2026/gpkg"
DEFAULT_CFE_DIR     = "/mnt/disk2/suma_helen_poster/cfe_py"
DEFAULT_CONFIG_FILE = "/mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json"
DEFAULT_OUT_DIR     = "/mnt/disk2/svyas_clean_run/forecast_al_storm_2026"

# ---------------------------------------------------------------------------
# Forecast settings
# ---------------------------------------------------------------------------
# Which HRRR initialization to use for the 18h forecast
# Use the most recent available: Oct 8 12z gives valid Oct 8 13:00 -> Oct 9 06:00
# Adjust to the latest HRRR run available at your run time.
FORECAST_INIT_TIME  = "2026-10-08T12"   # HRRR init time (UTC)
LEAD_HOURS          = 18

# Which DA member states to start from (last hour saved in the storm window)
FORECAST_STATE_TIME = "2026-10-08 17:00:00"

# Alabama bounding box
LAT_MIN, LAT_MAX = 29.5, 35.5
LON_MIN, LON_MAX = -89.0, -84.0

# Perturbation (same as DA run -- HRRR-tier)
PRECIP_SIGMA = 0.15
PET_SIGMA    = 0.10
ALBEDO       = 0.20
ALPHA_PT     = 1.26

AL_GAUGES = [
    "02361000", "02369800", "02371500", "02372250", "02374500", "02408540",
    "02422500", "02450250", "02464000", "02465493", "02469800", "03574500",
]

DYNAMICAL_DATASET = "noaa-hrrr-forecast-18-hour-virtual"

HRRR_VAR_MAP = {
    "total_precipitation_surface":               "APCP_1hr_acc_fcst",
    "downward_short_wave_radiation_flux_surface": "DSWRF",
    "temperature_2m":                            "TMP",   # °C -> saved as K
}


# ---------------------------------------------------------------------------
# HRRR forecast download  (f01..f18 from one init time)
# ---------------------------------------------------------------------------
def priestley_taylor_pet(srad_wm2, T_kelvin):
    T = T_kelvin - 273.15
    delta = 4098 * (0.6108 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3) ** 2
    gamma = 0.0638
    Rn_mj = np.maximum((1.0 - ALBEDO) * srad_wm2, 0.0) * 0.0036
    lam = 2.501 - 0.002361 * T
    return np.maximum(ALPHA_PT * (delta / (delta + gamma)) * Rn_mj / lam, 0.0)


def get_centroids(gpkg_dir):
    try:
        import geopandas as gpd
        centroids = {}
        for gauge in AL_GAUGES:
            for gpkg in sorted(Path(gpkg_dir).glob(f"*{gauge}*.gpkg")):
                div = gpd.read_file(str(gpkg), layer="divides").to_crs("EPSG:4326")
                for _, row in div.iterrows():
                    cat_id = str(row["divide_id"]).strip()
                    if cat_id.startswith("cat-"):
                        c = row.geometry.centroid
                        centroids[cat_id] = (float(c.y), float(c.x))
        return centroids
    except ImportError:
        pass
    # sqlite fallback
    centroids = {}
    for gauge in AL_GAUGES:
        for gpkg in sorted(Path(gpkg_dir).glob(f"*{gauge}*.gpkg")):
            con = sqlite3.connect(str(gpkg))
            try:
                rows = pd.read_sql(
                    "SELECT divide_id, "
                    "(ST_MinX(geom)+ST_MaxX(geom))/2 AS lon, "
                    "(ST_MinY(geom)+ST_MaxY(geom))/2 AS lat FROM divides", con)
                for _, r in rows.iterrows():
                    cat_id = str(r["divide_id"]).strip()
                    if cat_id.startswith("cat-"):
                        centroids[cat_id] = (float(r["lat"]), float(r["lon"]))
            except Exception:
                pass
            finally:
                con.close()
    return centroids


def download_hrrr_forecast(init_time_str, gpkg_dir, cache_dir):
    """Download HRRR f01-f18 from dynamical.org for Alabama catchments.

    Returns dict: cat_id -> DataFrame with columns
        valid_time (str), lead_h (int), total_precipitation (mm),
        potential_evaporation (mm/h)
    """
    cache_file = Path(cache_dir) / f"hrrr_fc_{init_time_str.replace(':', '').replace('-', '')}.parquet"
    if cache_file.exists():
        print(f"  Reusing cached HRRR forecast: {cache_file.name}")
        df = pd.read_parquet(cache_file)
        return {cat: grp.reset_index(drop=True)
                for cat, grp in df.groupby("cat_id")}

    print(f"  Opening dynamical.org HRRR for init_time={init_time_str} ...")
    try:
        import dynamical_catalog
        ds_full = dynamical_catalog.open(DYNAMICAL_DATASET, chunks=None)
    except ImportError:
        import icechunk, pystac
        catalog = pystac.Catalog.from_file("https://stac.dynamical.org/catalog.json")
        collection = catalog.get_child(DYNAMICAL_DATASET)
        asset = collection.assets["icechunk-https"]
        authorize = icechunk.containers_credentials(
            {"s3://noaa-hrrr-bdp-pds/": icechunk.s3_anonymous_credentials()})
        repo = icechunk.Repository.open(
            icechunk.http_storage(asset.href),
            authorize_virtual_chunk_access=authorize)
        ds_full = xr.open_zarr(repo.readonly_session("main").store, chunks=None)

    ds = ds_full[list(HRRR_VAR_MAP.keys())]
    ds = ds.rename(HRRR_VAR_MAP)
    ds = ds.sel(init_time=init_time_str)   # single init time
    # select lead_times 1h..18h
    leads = [np.timedelta64(h, "h") for h in range(1, LEAD_HOURS + 1)]
    ds = ds.sel(lead_time=leads)

    lat_name = next(c for c in ds.coords if "lat" in c.lower())
    lon_name = next(c for c in ds.coords if "lon" in c.lower())
    lat = ds[lat_name].values
    lon = ds[lon_name].values
    is_2d = lat.ndim == 2

    # Alabama subset
    if is_2d:
        mask = ((lat >= LAT_MIN) & (lat <= LAT_MAX) & (lon >= LON_MIN) & (lon <= LON_MAX))
        ys, xs = np.where(mask)
        ds = ds.isel(y=slice(int(ys.min()), int(ys.max()) + 1),
                     x=slice(int(xs.min()), int(xs.max()) + 1))
        lat = ds[lat_name].values
        lon = ds[lon_name].values
    else:
        ds = ds.isel(**{lat_name: (lat >= LAT_MIN) & (lat <= LAT_MAX),
                        lon_name: (lon >= LON_MIN) & (lon <= LON_MAX)})
        lat = ds[lat_name].values
        lon = ds[lon_name].values

    print(f"  Loading forecast slice into memory ...")
    ds = ds.load()

    centroids = get_centroids(gpkg_dir)
    init_ts = pd.Timestamp(init_time_str)
    lead_hours = list(range(1, LEAD_HOURS + 1))

    rows = []
    for cat_id, (clat, clon) in centroids.items():
        if is_2d:
            dist = (lat - clat) ** 2 + (lon - clon) ** 2
            yi, xi = np.unravel_index(dist.argmin(), dist.shape)
        else:
            yi = int(np.argmin(np.abs(lat - clat)))
            xi = int(np.argmin(np.abs(lon - clon)))

        for li, lh in enumerate(lead_hours):
            valid_ts = (init_ts + pd.Timedelta(hours=lh))
            if is_2d:
                prec = float(ds["APCP_1hr_acc_fcst"].isel(lead_time=li, y=yi, x=xi).values)
                dswrf = float(ds["DSWRF"].isel(lead_time=li, y=yi, x=xi).values)
                tmp_c = float(ds["TMP"].isel(lead_time=li, y=yi, x=xi).values)
            else:
                prec = float(ds["APCP_1hr_acc_fcst"].isel(lead_time=li, **{lat_name: yi, lon_name: xi}).values)
                dswrf = float(ds["DSWRF"].isel(lead_time=li, **{lat_name: yi, lon_name: xi}).values)
                tmp_c = float(ds["TMP"].isel(lead_time=li, **{lat_name: yi, lon_name: xi}).values)
            pet = float(priestley_taylor_pet(np.array([dswrf]), np.array([tmp_c + 273.15]))[0])
            rows.append({"cat_id": cat_id, "valid_time": valid_ts.strftime("%Y-%m-%d %H:%M:%S"),
                         "lead_h": lh, "precip_mm": max(prec, 0.0), "pet_mm": pet})

    df = pd.DataFrame(rows)
    df.to_parquet(cache_file, index=False)
    print(f"  Cached forecast forcing: {cache_file.name}")
    return {cat: grp.reset_index(drop=True) for cat, grp in df.groupby("cat_id")}


# ---------------------------------------------------------------------------
# CFE helpers
# ---------------------------------------------------------------------------
def write_bmi_config(template_file, best_params, forcing_csv, out_path):
    with open(template_file) as f:
        cfg = json.load(f)
    cfg["forcing_file"]          = str(forcing_csv)
    cfg["soil_params"]["bb"]     = best_params["bb"]
    cfg["soil_params"]["smcmax"] = best_params["smcmax"]
    cfg["soil_params"]["satdk"]  = best_params["satdk"]
    cfg["slop"]                  = best_params["slop"]
    cfg["max_gw_storage"]        = best_params["max_gw_storage"]
    cfg["expon"]                 = best_params["expon"]
    cfg["Cgw"]                   = best_params["Cgw"]
    cfg["K_lf"]                  = best_params["K_lf"]
    cfg["K_nash"]                = best_params["K_nash"]
    cfg["partition_scheme"]      = ("Schaake" if best_params["scheme"] <= 0.5
                                    else "Xinanjiang")
    with open(out_path, "w") as f:
        json.dump(cfg, f)
    return out_path


def init_member_from_state(cfg_path, forcing_df, state_row, rng):
    """Initialize one CFE member from a saved state row."""
    bmi_forcing = forcing_df.rename(columns={"date": "time"})
    m = bmi_cfe.BMI_CFE(cfg_file=cfg_path)
    m.load_forcing_file = (lambda self=m: setattr(self, "forcing_data", bmi_forcing))
    m.initialize()
    sm_max = m.soil_reservoir["storage_max_m"]
    gw_max = m.gw_reservoir["storage_max_m"]
    m.soil_reservoir["storage_m"] = float(np.clip(state_row["soil_m"], 0.0, sm_max))
    m.gw_reservoir["storage_m"]   = float(np.clip(state_row["gw_m"], 0.0, gw_max))
    m.nash_storage[0] = float(max(state_row["nash0"], 0.0))
    m.nash_storage[1] = float(max(state_row["nash1"], 0.0))
    return m


def plot_forecast(valid_times, q_members, cat_id, init_time_str, out_png):
    q = np.array(q_members)   # shape: (18, n_members)
    mean = q.mean(axis=1)
    std  = q.std(axis=1, ddof=1)
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.fill_between(valid_times, np.clip(mean - 2*std, 0, None), mean + 2*std,
                    color="tomato", alpha=0.25, label="±2 std")
    ax.plot(valid_times, mean, "tomato", lw=2, label="ensemble mean")
    ax.set_ylabel("Runoff (mm/h)")
    ax.set_title(f"{cat_id} | 18h forecast | HRRR init {init_time_str}")
    ax.legend()
    twin = ax.twinx()
    # plot nothing on twin but keep the axis for future precip overlay
    plt.tight_layout()
    plt.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Main per-catchment forecast
# ---------------------------------------------------------------------------
def run_catchment_forecast(cat_id, forecast_forcing, args):
    t0 = time.time()
    out_dir = Path(args.out_dir) / cat_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load DA member states at the requested state time
    state_time = getattr(args, "forecast_state_time", FORECAST_STATE_TIME)
    states_pq = Path(args.da_dir) / cat_id / f"{cat_id}_member_states_al_storm.parquet"
    if not states_pq.exists():
        raise FileNotFoundError(f"No member states: {states_pq}")
    states_df = pd.read_parquet(states_pq)
    states_t  = states_df[states_df["time"] == state_time]
    if len(states_t) == 0:
        avail = sorted(states_df["time"].unique())
        print(f"  {cat_id}: state time '{state_time}' not found, using {avail[-1]}")
        states_t = states_df[states_df["time"] == avail[-1]]
    states_t = states_t.sort_values("member").reset_index(drop=True)
    n_members = min(args.n_members, len(states_t))
    states_t = states_t.iloc[:n_members]
    print(f"  Loaded {n_members} member states at {states_t['time'].iloc[0]}")

    # Calibrated params
    params_file = Path(args.params_dir) / cat_id / f"{cat_id}_best_params.json"
    if not params_file.exists():
        params_file = Path(args.params_dir) / f"{cat_id}_best_params.json"
    with open(params_file) as f:
        best_params = json.load(f)["best_parameters"]

    # Forecast forcing for this catchment
    if cat_id not in forecast_forcing:
        raise KeyError(f"No forecast forcing for {cat_id}")
    fc_df = forecast_forcing[cat_id]
    valid_times = pd.to_datetime(fc_df["valid_time"])

    # Need a dummy forcing CSV for BMI config (path must exist; we set values via set_value)
    # Reuse the DA hrrr_combined CSV if it exists, else create minimal dummy
    dummy_csv = Path(args.da_dir) / cat_id / f"{cat_id}_hrrr_combined.csv"
    if not dummy_csv.exists():
        dummy_csv = out_dir / f"{cat_id}_dummy_forcing.csv"
        dummy_csv.write_text("time,APCP_1hr_acc_fcst,DSWRF,TMP\n2026-10-01 01:00:00,0,0,280\n")

    cfg_path = write_bmi_config(args.config_file, best_params, dummy_csv,
                                out_dir / f"{cat_id}_bmi_config_fc.json")

    # Build a minimal forcing_df for ensemble initialization (not actually used for stepping)
    forcing_stub = pd.DataFrame({"time": [str(t) for t in valid_times],
                                 "total_precipitation": fc_df["precip_mm"].values,
                                 "potential_evaporation": fc_df["pet_mm"].values})

    rng = np.random.default_rng(args.rng_seed)
    print(f"  Initializing {n_members} members from DA states...")
    models = [init_member_from_state(str(cfg_path), forcing_stub, states_t.iloc[i], rng)
              for i in range(n_members)]

    # Run 18h forecast
    q_all = np.zeros((LEAD_HOURS, n_members), dtype=np.float32)
    for li, row in enumerate(fc_df.itertuples()):
        precip = float(row.precip_mm)
        pet    = float(row.pet_mm)
        p_pert = precip * rng.lognormal(-0.5 * PRECIP_SIGMA**2, PRECIP_SIGMA, n_members)
        e_pert = np.maximum(pet * (1.0 + PET_SIGMA * rng.standard_normal(n_members)), 0.0)
        for i, m in enumerate(models):
            m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                        float(p_pert[i]) / 1000.0)
            m.set_value("water_potential_evaporation_flux",
                        float(e_pert[i]) / 1000.0 / 3600.0)
            m.update()
            q_all[li, i] = float(m.get_value("land_surface_water__runoff_depth") * 1000.0)

    for m in models:
        m.finalize()

    # Save ensemble runoff
    cols = [f"member_{i:04d}" for i in range(n_members)]
    fc_out = pd.DataFrame(q_all, columns=cols, dtype=np.float32)
    fc_out.insert(0, "valid_time", [t.strftime("%Y-%m-%d %H:%M:%S") for t in valid_times])
    fc_out.insert(1, "lead_h", list(range(1, LEAD_HOURS + 1)))
    tag = f"_{args.out_tag}" if getattr(args, "out_tag", "") else ""
    fc_out.to_parquet(out_dir / f"{cat_id}_forecast18h{tag}.parquet", index=False)

    # Mean/std summary CSV
    mean = q_all.mean(axis=1)
    std  = q_all.std(axis=1, ddof=1)
    summary = pd.DataFrame({"valid_time": fc_out["valid_time"],
                             "lead_h": fc_out["lead_h"],
                             "q_mean_mm_h": mean, "q_std_mm_h": std})
    summary.to_csv(out_dir / f"{cat_id}_forecast18h_mean.csv", index=False)

    plot_forecast(valid_times, q_all, cat_id, FORECAST_INIT_TIME,
                  out_dir / f"{cat_id}_forecast18h.png")

    print(f"  peak q_mean={mean.max():.3f} mm/h at lead {int(np.argmax(mean))+1}h | "
          f"{time.time()-t0:.1f}s")
    return {"cat_id": cat_id, "peak_q_mean_mm_h": float(mean.max()),
            "peak_lead_h": int(np.argmax(mean)) + 1}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    global bmi_cfe, FORECAST_INIT_TIME
    parser = argparse.ArgumentParser(
        description="18h ensemble forecast for Alabama storm 2026.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cat-id",       required=True,
                        help="catchment id or 'all'")
    parser.add_argument("--n-members",    type=int, default=600)
    parser.add_argument("--rng-seed",     type=int, default=99)
    parser.add_argument("--da-dir",       default=DEFAULT_DA_DIR)
    parser.add_argument("--params-dir",   default=DEFAULT_PARAMS_DIR)
    parser.add_argument("--gpkg-dir",     default=DEFAULT_GPKG_DIR)
    parser.add_argument("--cfe-dir",      default=DEFAULT_CFE_DIR)
    parser.add_argument("--config-file",  default=DEFAULT_CONFIG_FILE)
    parser.add_argument("--out-dir",      default=DEFAULT_OUT_DIR)
    parser.add_argument("--forecast-init-time", default=FORECAST_INIT_TIME,
                        help="HRRR init time, e.g. 2026-10-09T12")
    parser.add_argument("--forecast-state-time", default=FORECAST_STATE_TIME,
                        help="DA member state timestamp to restart from, e.g. '2026-10-09 12:00:00'")
    parser.add_argument("--out-tag", default="",
                        help="suffix appended to output filenames, e.g. '2026100912' → _forecast18h_2026100912.parquet")
    args = parser.parse_args()

    FORECAST_INIT_TIME = args.forecast_init_time
    sys.path.insert(0, args.cfe_dir)
    import bmi_cfe as _bmi_cfe
    bmi_cfe = _bmi_cfe

    if args.cat_id == "all":
        cat_ids = sorted(
            p.name for p in Path(args.da_dir).iterdir()
            if p.is_dir() and (p / f"{p.name}_member_states_al_storm.parquet").exists())
        print(f"Running forecast for {len(cat_ids)} catchments")
    else:
        cat_ids = [args.cat_id]

    # Download forecast forcing once for all catchments
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    print(f"Downloading HRRR forecast (init={FORECAST_INIT_TIME}, f01-f{LEAD_HOURS})...")
    forecast_forcing = download_hrrr_forecast(
        FORECAST_INIT_TIME, args.gpkg_dir, Path(args.out_dir))

    failed = []
    summaries = []
    for cat_id in cat_ids:
        try:
            s = run_catchment_forecast(cat_id, forecast_forcing, args)
            summaries.append(s)
        except Exception:
            print(f"\nFAILED: {cat_id}")
            traceback.print_exc()
            failed.append(cat_id)

    if summaries:
        df = pd.DataFrame(summaries)
        df.to_csv(Path(args.out_dir) / "forecast18h_summary.csv", index=False)
        print(f"\nPeak runoff summary saved to {args.out_dir}/forecast18h_summary.csv")
        print(f"Median peak q_mean: {df['peak_q_mean_mm_h'].median():.4f} mm/h "
              f"at lead {df['peak_lead_h'].median():.0f}h")

    if failed:
        print(f"\n{len(failed)}/{len(cat_ids)} FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"\nAll {len(cat_ids)} catchment(s) done.")


if __name__ == "__main__":
    main()
