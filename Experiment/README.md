## Setup — build the environment (once)
Requires **Python 3.10 with t-route**. 
`troute-network` / `troute-routing` / `troute-config` and `bmipy`.
```bash
conda create -n t-route python=3.10 -y && conda activate t-route
pip install numpy==1.26.4 pandas==2.2.0 pyarrow matplotlib netcdf4 pyyaml bmipy
pip install -r requirements.txt                 # repo extras, run from the repo root

# t-route routing engine — clone and build per its README:
#   https://github.com/NOAA-OWP/t-route
```
The CFE model code is needed for step 1 (the EnKF/forecast). On this server it is at
`/mnt/disk2/suma_helen_poster/cfe_py`.

## Each shell (before running)
```bash
conda activate t-route
PY=python                                       # your t-route env's python
export MPLBACKEND=Agg OMP_NUM_THREADS=1
cd <repo>/da_methods/forecast_experiment        # where you cloned it
```

> Paths below are gauge **03463300** (also the built-in defaults). For a different gauge,
> change the paths in each command, and edit `CATS` / `TERMINAL_INT` / `WATERSHED_AREA_KM2`
> at the top of `route_forecast.py`.

## Step 1 — EnKF + 18 h forecasts (DA + open loop) **and** member_q trajectories, all 21 catchments
Writes per-catchment forecast parquets to `--out-dir/<cat>/`, and (via `--save-members`) the
per-member runoff trajectories to `--da-traj-dir` / `--ol-traj-dir` (read by the routing warm-up).
Assimilation always uses the **NWM** operational forcing (`--test-forcing-dir1/2`); the 18 h
forecast leads use **HRRR** (`--forecast-forcing-dir`) — drop it to fall back to NWM perfect-forcing.
```bash
$PY run_cfe_enkf_forecast_wide_spread.py \
    --cat-id all --n-members 600 --rng-seed 42 --lead-hours 18 \
    --issue-start "2024-09-20 00:00:00" --issue-end "2024-09-29 23:00:00" \
    --obs-dir      /mnt/disk2/1400_sites_helene/catchment_ts_03463300_dynamic_variance_rekrig \
    --params-dir   /mnt/disk2/1400_sites_helene/da_results_dynamic_novrugt_seeded \
    --cfe-dir      /mnt/disk2/suma_helen_poster/cfe_py \
    --config-file  /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \
    --test-forcing-dir1 /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2023_2024_feb/forcings \
    --test-forcing-dir2 /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_nwmoperational/03463300/2024_feb_2025_sep/forcings \
    --forecast-forcing-dir /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/output_03463300_hrrr \
    --save-members \
    --da-traj-dir  /mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_da \
    --ol-traj-dir  /mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_ol \
    --out-dir      /mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_hrrr
```

## Step 2 — route every member to the gauge
Writes `routed_forecast_leadtime_{da,openloop}.parquet`
(columns: issue_time, lead_hour, member_0000..0599; discharge in m³/s) to `--out-dir`.
The warm-up trajectories are the NWM-assimilation `member_q` (same for NWM or HRRR forecasts).
```bash
$PY route_forecast.py \
    --scenario both --warmup-mode trajectory --warmup-h 48 --procs 20 \
    --gpkg        /mnt/disk1/usgs_streamflow_allgauges/subdaily_15min/test/gage-03463300_subset.gpkg \
    --fc-dir      /mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_hrrr \
    --da-traj-dir /mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_da \
    --ol-traj-dir /mnt/disk2/1400_sites_helene/da_results_enkf600_wide_spread_traj_ol \
    --out-dir     /mnt/disk2/1400_sites_helene/da_results_enkf600_forecast18h_hrrr/routed
```

---
- To run unattended/in the background, prefix a command with `nohup ... > run.log 2>&1 &`.

## CFE parameters — where they come from / how to re-calibrate

The `--params-dir` in Step 1 points to pre-calibrated CFE best-parameter JSON files
(`<cat-id>_best_params.json`, one per catchment).  Two sets are available:

| Set | Calibration | Use when |
|-----|-------------|----------|
| `held_in_gauge/params/` | gauge 03463300 **included** in calibration | default — best skill at the target gauge |
| `held_out_gauge/params/` | gauge 03463300 **withheld** | evaluating out-of-sample generalization |

Both sets live in the companion calibration repo.  To **re-run or extend calibration**,
follow the workflow in:

> **[calibrate-cfe PR #25](https://github.com/NWC-CUAHSI-Summer-Institute/calibrate-cfe/pull/25)**
> (`multicatchment-calibration-held-in-out` branch)

That repo uses DDS (Dynamically Dimensioned Search, 1 000 iterations via spotpy) against
the 2020–2022 retrospective period and writes one `<cat-id>_best_params.json` per catchment.
Clone it, checkout the branch above, and follow its README.

