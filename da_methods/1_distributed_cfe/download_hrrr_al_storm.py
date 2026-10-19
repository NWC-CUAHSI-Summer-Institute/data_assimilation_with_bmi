#!/usr/bin/env python3
"""
download_hrrr_al_storm.py -- Download HRRR from dynamical.org and extract
per-catchment forcing CSVs for the Alabama storm 2026 DA experiment.

SOURCE   : dynamical.org  --  NOAA HRRR forecast, 18-hour (Icechunk/zarr)
           Catalog page: https://dynamical.org/catalog/noaa-hrrr-forecast-18-hour-virtual/
DOMAIN   : Alabama storm catchments (1072 cats across 12 gauges)
PERIOD   : 2026-10-01 to 2026-10-08 (init_time), lead_time = 1h
VARIABLES (dynamical.org name → our column):
           total_precipitation_surface              → APCP_1hr_acc_fcst  (kg m⁻² = mm, 1-hr acc)
           downward_short_wave_radiation_flux_surface → DSWRF            (W m⁻²)
           temperature_2m                           → TMP               (°C → converted to K)

OUTPUT (matches build_hrrr_forcing() in run_cfe_enkf_al_storm.py):
  <OUT_DIR>/camels_YYYYMMDD/cat-<id>.csv   (24 rows per day)
  Columns: time, APCP_1hr_acc_fcst, DSWRF, TMP   (TMP in Kelvin)

HOW TO RUN
    # install deps once:
    pip install dynamical-catalog icechunk geopandas

    # dry-run: open catalog, print variables, don't download
    python3 download_hrrr_al_storm.py --dry-run

    # full download (writes ~1072 x 8 = ~8576 CSVs)
    python3 download_hrrr_al_storm.py
"""

import argparse
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
GPKG_DIR = Path("/mnt/disk2/al_storm_2026/gpkg")
OUT_DIR  = Path("/mnt/disk2/svyas_clean_run/hrrr_al_storm_2026")

# ---------------------------------------------------------------------------
# dynamical.org dataset slug (18-hour HRRR forecast, virtual Icechunk store)
# Docs: https://dynamical.org/catalog/noaa-hrrr-forecast-18-hour-virtual/
# ---------------------------------------------------------------------------
DYNAMICAL_DATASET = "noaa-hrrr-forecast-18-hour-virtual"

# Mapping: dynamical.org variable name → output CSV column name
HRRR_VAR_MAP = {
    "total_precipitation_surface":               "APCP_1hr_acc_fcst",  # kg/m² = mm (1-hr acc at f01)
    "downward_short_wave_radiation_flux_surface": "DSWRF",             # W/m²
    "temperature_2m":                            "TMP",                # °C in source → saved as K
}

# ---------------------------------------------------------------------------
# Time window  (init_time range; valid_time = init_time + 1h)
# ---------------------------------------------------------------------------
INIT_START = "2026-10-09T00"   # extend: only download new dates
INIT_END   = "2026-10-12T23"   # storm peaks Oct 10-12

# Alabama bounding box (degrees, slight padding)
LAT_MIN, LAT_MAX = 29.5, 35.5
LON_MIN, LON_MAX = -89.0, -84.0

# Gauge list
AL_GAUGES = [
    "02361000", "02369800", "02371500", "02372250", "02374500", "02408540",
    "02422500", "02450250", "02464000", "02465493", "02469800", "03574500",
]


# ---------------------------------------------------------------------------
# Catchment centroids
# ---------------------------------------------------------------------------
def get_centroids_geopandas(gpkg_dir):
    import geopandas as gpd
    centroids = {}
    for gauge in AL_GAUGES:
        matches = sorted(gpkg_dir.glob(f"*{gauge}*.gpkg"))
        if not matches:
            print(f"  WARNING: no GPKG for gauge {gauge}")
            continue
        try:
            div = gpd.read_file(str(matches[0]), layer="divides").to_crs("EPSG:4326")
            for _, row in div.iterrows():
                cat_id = str(row["divide_id"]).strip()
                if cat_id.startswith("cat-"):
                    c = row.geometry.centroid
                    centroids[cat_id] = (float(c.y), float(c.x))
        except Exception as e:
            print(f"  WARNING: {matches[0].name}: {e}")
    return centroids


def get_centroids_sqlite(gpkg_dir):
    """Fallback when geopandas is unavailable."""
    centroids = {}
    for gauge in AL_GAUGES:
        matches = sorted(gpkg_dir.glob(f"*{gauge}*.gpkg"))
        if not matches:
            continue
        con = sqlite3.connect(str(matches[0]))
        try:
            bbox = pd.read_sql(
                "SELECT i.min_x, i.max_x, i.min_y, i.max_y, "
                "d.divide_id FROM rtree_divides_geom i "
                "JOIN divides d ON d.fid = i.id", con)
            for _, row in bbox.iterrows():
                cat_id = str(row["divide_id"]).strip()
                if cat_id.startswith("cat-"):
                    centroids[cat_id] = (
                        (row["min_y"] + row["max_y"]) / 2,
                        (row["min_x"] + row["max_x"]) / 2,
                    )
        except Exception as e:
            print(f"  WARNING: sqlite centroid fallback failed for {matches[0].name}: {e}")
        finally:
            con.close()
    return centroids


def load_centroids(gpkg_dir):
    try:
        import geopandas  # noqa: F401
        c = get_centroids_geopandas(gpkg_dir)
    except ImportError:
        print("  geopandas not available, using sqlite fallback")
        c = get_centroids_sqlite(gpkg_dir)
    print(f"  {len(c)} catchment centroids loaded")
    return c


# ---------------------------------------------------------------------------
# Open dynamical.org HRRR dataset
# ---------------------------------------------------------------------------
def open_dynamical_hrrr():
    """Open HRRR 18-h forecast from dynamical.org via dynamical-catalog."""
    try:
        import dynamical_catalog  # pip install dynamical-catalog
        ds = dynamical_catalog.open(DYNAMICAL_DATASET, chunks=None)
        print(f"  Opened via dynamical_catalog: {DYNAMICAL_DATASET}")
        return ds
    except ImportError:
        pass

    # Fallback: open via icechunk + pystac directly
    try:
        import icechunk
        import pystac
        print("  dynamical_catalog not found, trying icechunk + pystac...")
        catalog = pystac.Catalog.from_file("https://stac.dynamical.org/catalog.json")
        collection = catalog.get_child(DYNAMICAL_DATASET)
        asset = collection.assets["icechunk-https"]
        authorize = icechunk.containers_credentials(
            {"s3://noaa-hrrr-bdp-pds/": icechunk.s3_anonymous_credentials()}
        )
        repo = icechunk.Repository.open(
            icechunk.http_storage(asset.href),
            authorize_virtual_chunk_access=authorize,
        )
        session = repo.readonly_session("main")
        ds = xr.open_zarr(session.store, chunks=None)
        print(f"  Opened via icechunk pystac")
        return ds
    except ImportError:
        pass

    raise ImportError(
        "Cannot open dynamical.org dataset. Install dependencies:\n"
        "  pip install dynamical-catalog\n"
        "  OR: pip install icechunk pystac"
    )


# ---------------------------------------------------------------------------
# Download + extract
# ---------------------------------------------------------------------------
def nearest_idx_2d(lat2d, lon2d, clat, clon):
    dist = (lat2d - clat) ** 2 + (lon2d - clon) ** 2
    return np.unravel_index(dist.argmin(), dist.shape)


def nearest_idx_1d(lat1d, lon1d, clat, clon):
    yi = int(np.argmin(np.abs(lat1d - clat)))
    xi = int(np.argmin(np.abs(lon1d - clon)))
    return yi, xi


def download_and_extract_all_leads(centroids, out_dir, lead_max=18, dry_run=False):
    """Download HRRR leads 1..lead_max and write fh01..fhNN/camels_YYYYMMDD/ structure.

    Each fh{L}/camels_{valid_date}/{cat}.csv contains all rows where
    valid_time = init_time + L falls on that calendar day, i.e. the row's
    'time' column is the VALID TIME (not the init time).  This is exactly the
    structure expected by run_cfe_enkf_forecast_hrrr_multihour_al_storm.py.
    """
    ds_full = open_dynamical_hrrr()
    missing = [v for v in HRRR_VAR_MAP if v not in ds_full.data_vars]
    if missing:
        print(f"WARNING: variables not found: {missing}")
        if not dry_run:
            sys.exit(1)

    if dry_run:
        print(f"\n[DRY RUN --all-leads] Would write fh01..fh{lead_max:02d}/camels_YYYYMMDD/")
        print(f"  init_time [{INIT_START} .. {INIT_END}], lead_time 1h-{lead_max}h")
        print(f"  {len(centroids)} catchments")
        return

    ds_vars = ds_full[list(HRRR_VAR_MAP.keys())]
    # select all needed lead_times at once
    lead_tds = [np.timedelta64(L, "h") for L in range(1, lead_max + 1)]
    available_leads = ds_vars.lead_time.values
    sel_leads = [lt for lt in lead_tds if lt in available_leads]
    if not sel_leads:
        raise ValueError(f"No lead_times 1h-{lead_max}h found in dataset "
                         f"(available: {available_leads[:5]}...)")
    ds = ds_vars.sel(lead_time=sel_leads)
    ds = ds.sel(init_time=slice(INIT_START, INIT_END))

    lat_name = next((c for c in ds.coords if "lat" in c.lower()), None)
    lon_name = next((c for c in ds.coords if "lon" in c.lower()), None)
    lat = ds[lat_name].values
    lon = ds[lon_name].values
    is_2d = lat.ndim == 2

    if is_2d:
        mask = ((lat >= LAT_MIN) & (lat <= LAT_MAX) &
                (lon >= LON_MIN) & (lon <= LON_MAX))
        ys, xs = np.where(mask)
        y_sl = slice(int(ys.min()), int(ys.max()) + 1)
        x_sl = slice(int(xs.min()), int(xs.max()) + 1)
        ds = ds.isel(y=y_sl, x=x_sl)
        lat = ds[lat_name].values
        lon = ds[lon_name].values
    else:
        ds = ds.isel(**{lat_name: (lat >= LAT_MIN) & (lat <= LAT_MAX),
                        lon_name: (lon >= LON_MIN) & (lon <= LON_MAX)})
        lat = ds[lat_name].values
        lon = ds[lon_name].values

    print(f"  Loading all leads 1h-{lead_max}h into memory...")
    ds = ds.load()
    print(f"  Loaded. {dict(ds.dims)}")

    init_times  = pd.to_datetime(ds.init_time.values)
    lead_values = ds.lead_time.values                 # timedelta64 array
    out_col_names = list(HRRR_VAR_MAP.values())

    # For each (init, lead) → valid_time; organize into fh{L}/camels_{valid_date}/
    for lidx, ltd in enumerate(lead_values):
        lead_h = int(ltd / np.timedelta64(1, "h"))
        lead_dir_name = f"fh{lead_h:02d}"
        print(f"\n  Writing {lead_dir_name} ({len(init_times)} init times x {len(centroids)} catchments)...")

        valid_times = init_times + pd.Timedelta(hours=lead_h)
        # group init indices by valid_date
        valid_dates = [vt.date() for vt in valid_times]
        unique_dates = sorted(set(valid_dates))

        for vdate in unique_dates:
            vday_str = vdate.strftime("%Y%m%d")
            vday_dir = out_dir / lead_dir_name / f"camels_{vday_str}"
            vday_dir.mkdir(parents=True, exist_ok=True)
            day_idxs = [i for i, d in enumerate(valid_dates) if d == vdate]

            n_written = 0
            for cat_id, (clat, clon) in centroids.items():
                out_csv = vday_dir / f"{cat_id}.csv"
                if out_csv.exists():
                    continue
                if is_2d:
                    yi, xi = nearest_idx_2d(lat, lon, clat, clon)
                else:
                    yi, xi = nearest_idx_1d(lat, lon, clat, clon)

                rows_txt = [f"time,{','.join(out_col_names)}"]
                for tidx in day_idxs:
                    valid_t = valid_times[tidx].strftime("%Y-%m-%d %H:%M:%S")
                    vals = []
                    for zarr_var, col in HRRR_VAR_MAP.items():
                        if is_2d:
                            v = float(ds[zarr_var].isel(
                                lead_time=lidx, init_time=tidx, y=yi, x=xi).values)
                        else:
                            v = float(ds[zarr_var].isel(
                                lead_time=lidx, init_time=tidx,
                                **{lat_name: yi, lon_name: xi}).values)
                        if col == "APCP_1hr_acc_fcst":
                            v = max(v, 0.0)
                        elif col == "TMP":
                            v = v + 273.15
                        vals.append(f"{v:.7g}")
                    rows_txt.append(f"{valid_t},{','.join(vals)}")

                out_csv.write_text("\n".join(rows_txt) + "\n")
                n_written += 1
            print(f"    {lead_dir_name}/camels_{vday_str}: {n_written} CSVs")

    print(f"\nAll-leads download complete. Output: {out_dir}")
    print(f"Next step: run_cfe_enkf_forecast_hrrr_multihour_al_storm.py "
          f"--forecast-forcing-dir {out_dir}")


def download_and_extract(centroids, out_dir, dry_run=False):
    ds_full = open_dynamical_hrrr()

    print(f"\n  Dataset variables: {list(ds_full.data_vars)}")
    print(f"  Dimensions: {dict(ds_full.dims)}")

    # Verify all expected variables exist
    missing = [v for v in HRRR_VAR_MAP if v not in ds_full.data_vars]
    if missing:
        print(f"\n  WARNING: variables not found: {missing}")
        print(f"  Available: {list(ds_full.data_vars)}")
        print("  Update HRRR_VAR_MAP at the top of this script and re-run.")
        if dry_run:
            return
        sys.exit(1)

    if dry_run:
        print(f"\n[DRY RUN] Variable check passed. Would extract:")
        print(f"  {len(centroids)} catchments x init_time [{INIT_START} .. {INIT_END}]")
        print(f"  lead_time = 1h  =>  valid_time = init_time + 1h")
        print(f"  Output: {out_dir}/camels_YYYYMMDD/cat-*.csv")
        print(f"  TMP will be converted: °C + 273.15 → K")
        return

    # Select: lead_time=1h, init_time window, Alabama bbox
    ds = ds_full[list(HRRR_VAR_MAP.keys())]
    ds = ds.sel(lead_time=np.timedelta64(1, "h"))
    ds = ds.sel(init_time=slice(INIT_START, INIT_END))

    # Identify lat/lon coordinate names
    lat_name = next((c for c in ds.coords if "lat" in c.lower()), None)
    lon_name = next((c for c in ds.coords if "lon" in c.lower()), None)
    if lat_name is None or lon_name is None:
        raise KeyError(f"Cannot find lat/lon coords. Available: {list(ds.coords)}")

    lat = ds[lat_name].values
    lon = ds[lon_name].values
    is_2d = lat.ndim == 2

    # Spatial subset to Alabama bbox
    if is_2d:
        mask = ((lat >= LAT_MIN) & (lat <= LAT_MAX) &
                (lon >= LON_MIN) & (lon <= LON_MAX))
        ys, xs = np.where(mask)
        if len(ys) == 0:
            raise ValueError(f"No grid points in Alabama bbox. Check lat/lon range.")
        y_sl = slice(int(ys.min()), int(ys.max()) + 1)
        x_sl = slice(int(xs.min()), int(xs.max()) + 1)
        ds = ds.isel(y=y_sl, x=x_sl)
        lat = ds[lat_name].values
        lon = ds[lon_name].values
        print(f"  Alabama subset: {lat.shape[0]} x {lat.shape[1]} grid points")
    else:
        lat_mask = (lat >= LAT_MIN) & (lat <= LAT_MAX)
        lon_mask = (lon >= LON_MIN) & (lon <= LON_MAX)
        ds = ds.isel(**{lat_name: lat_mask, lon_name: lon_mask})
        lat = ds[lat_name].values
        lon = ds[lon_name].values
        print(f"  Alabama subset: {len(lat)} lat x {len(lon)} lon")

    print("  Loading HRRR Alabama subset into memory...")
    ds = ds.load()
    print("  Loaded.")

    # init_time dimension → group by calendar day
    init_times = pd.to_datetime(ds.init_time.values)
    out_col_names = list(HRRR_VAR_MAP.values())
    days = sorted(set(init_times.date))
    print(f"  {len(init_times)} init hours across {len(days)} days")

    for day in days:
        day_str = day.strftime("%Y%m%d")
        day_dir = out_dir / f"camels_{day_str}"
        day_dir.mkdir(parents=True, exist_ok=True)

        day_mask = init_times.date == day
        day_inits = init_times[day_mask]
        n_written = 0

        for cat_id, (clat, clon) in centroids.items():
            out_csv = day_dir / f"{cat_id}.csv"
            if out_csv.exists():
                continue

            if is_2d:
                yi, xi = nearest_idx_2d(lat, lon, clat, clon)
            else:
                yi, xi = nearest_idx_1d(lat, lon, clat, clon)

            rows = [f"time,{','.join(out_col_names)}"]
            for tidx, init_t in zip(np.where(day_mask)[0], day_inits):
                # valid_time = init_time + 1h  (that's the timestamp we store)
                valid_t = (init_t + pd.Timedelta("1h")).strftime("%Y-%m-%d %H:%M:%S")
                vals = []
                for zarr_var, col in HRRR_VAR_MAP.items():
                    if is_2d:
                        v = float(ds[zarr_var].isel(init_time=tidx, y=yi, x=xi).values)
                    else:
                        v = float(ds[zarr_var].isel(
                            init_time=tidx,
                            **{lat_name: yi, lon_name: xi}).values)
                    if col == "APCP_1hr_acc_fcst":
                        v = max(v, 0.0)          # kg/m² = mm, clip negative
                    elif col == "TMP":
                        v = v + 273.15           # °C → K (priestley_taylor_pet expects K)
                    vals.append(f"{v:.7g}")
                rows.append(f"{valid_t},{','.join(vals)}")

            out_csv.write_text("\n".join(rows) + "\n")
            n_written += 1

        print(f"  {day_str}: wrote {n_written} CSVs -> {day_dir.name}/")

    print(f"\nDone. Output: {out_dir}")
    print(f"Next step: python3 run_cfe_enkf_al_storm.py --cat-id all --hrrr-dir {out_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    global INIT_START, INIT_END
    parser = argparse.ArgumentParser(
        description="Download HRRR (dynamical.org) for Alabama storm 2026.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gpkg-dir",   default=str(GPKG_DIR))
    parser.add_argument("--out-dir",    default=str(OUT_DIR))
    parser.add_argument("--init-start", default=INIT_START,
                        help="first HRRR init time to download, e.g. 2026-10-09T00")
    parser.add_argument("--init-end",   default=INIT_END,
                        help="last HRRR init time to download, e.g. 2026-10-09T12")
    parser.add_argument("--dry-run",    action="store_true",
                        help="open catalog and check variables without downloading")
    parser.add_argument("--all-leads",  action="store_true",
                        help="download leads 1h-18h and write fh01..fh18/camels_YYYYMMDD/ "
                             "(needed by run_cfe_enkf_forecast_hrrr_multihour_al_storm.py)")
    parser.add_argument("--lead-max",   type=int, default=18,
                        help="highest lead to download with --all-leads (default 18)")
    args = parser.parse_args()

    INIT_START = args.init_start
    INIT_END   = args.init_end

    gpkg_dir = Path(args.gpkg_dir)
    out_dir  = Path(args.out_dir)

    print("Reading catchment centroids...")
    centroids = load_centroids(gpkg_dir)
    if not centroids:
        sys.exit("No catchments found -- check --gpkg-dir")

    out_dir.mkdir(parents=True, exist_ok=True)
    if args.all_leads:
        download_and_extract_all_leads(centroids, out_dir,
                                       lead_max=args.lead_max,
                                       dry_run=args.dry_run)
    else:
        download_and_extract(centroids, out_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
