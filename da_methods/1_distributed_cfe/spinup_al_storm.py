#!/usr/bin/env python3
"""
spinup_al_storm.py
Runs a single deterministic 1-year CFE spin-up for every Alabama catchment
using WY2026 NWM operational forcing (Oct 2025 → Oct 2026).

Saves final soil/gw/Nash states to <states-dir>/<cat_id>_spinup_states.json
so the DA + OL ensembles can start from proper, dry October 2026 conditions
rather than the wet Pigeon River template (soil=0.667).

USAGE
  python spinup_al_storm.py \\
      --nwm-wy2026-dir /mnt/disk2/al_storm_2026/nwm_operational/wy2026/forcings \\
      --params-dir     /mnt/disk2/svyas_clean_run/calibration_al_storm_2026 \\
      --config-file    /mnt/disk2/suma_helen_poster/run_gpu/cat_03463300_bmi_config_cfe.json \\
      --states-dir     /mnt/disk2/svyas_clean_run/spinup_states_al_storm_2026 \\
      --cat-id         all \\
      --procs          20

  # single catchment smoke test
  python spinup_al_storm.py --cat-id cat-1011916 --procs 1 ...
"""
import argparse
import json
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

# reuse Alabama DA engine (bmi_cfe, build_nwm_forcing, write_bmi_config, etc.)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cfe_enkf_al_storm as enkf

SPINUP_WY26_START = "2025-10-01 00:00:00"
SPINUP_WY26_END   = "2026-10-01 00:00:00"   # save states at this point

ALL_CATS_FILE = Path(__file__).parent / "al_storm_cat_ids.txt"

G = {}   # multiprocessing shared args


def run_spinup_catchment(cat_id):
    args       = G["args"]
    states_dir = Path(args.states_dir)
    out_json   = states_dir / f"{cat_id}_spinup_states.json"
    if out_json.exists() and not args.force:
        return cat_id, "skip"

    t0 = time.time()
    try:
        # ── scratch dir for BMI config + NWM cache ──────────────────────────
        scratch = states_dir / "scratch" / cat_id
        scratch.mkdir(parents=True, exist_ok=True)

        # ── load calibrated params ────────────────────────────────────────────
        pf = Path(args.params_dir) / cat_id / f"{cat_id}_best_params.json"
        if not pf.exists():
            pf = Path(args.params_dir) / f"{cat_id}_best_params.json"
        with open(pf) as f:
            best_params = json.load(f)["best_parameters"]

        # ── build WY2026 NWM forcing CSV ─────────────────────────────────────
        forcing_csv = enkf.build_nwm_forcing(cat_id, args.nwm_wy2026_dir, scratch)
        forcing     = enkf.load_forcing(forcing_csv)

        # ── slice to spin-up window ───────────────────────────────────────────
        spin = forcing[(forcing["date"] >= SPINUP_WY26_START)
                       & (forcing["date"] <= SPINUP_WY26_END)].reset_index(drop=True)
        if len(spin) < 100:
            return cat_id, f"ERROR: only {len(spin)} rows in wy2026 forcing"

        # ── write BMI config with calibrated params ───────────────────────────
        cfg_path = enkf.write_bmi_config(
            args.config_file, best_params, forcing_csv,
            scratch / f"{cat_id}_bmi_config_spinup.json")

        # ── initialise a single deterministic CFE model ───────────────────────
        bmi_forcing = spin.rename(columns={"date": "time"})
        m = enkf.bmi_cfe.BMI_CFE(cfg_file=cfg_path)
        m.load_forcing_file = (lambda self=m: setattr(self, "forcing_data", bmi_forcing))
        m.initialize()

        # ── step forward: one deterministic run, no DA, no noise ─────────────
        for _, row in spin.iterrows():
            precip_mmh = float(row["total_precipitation"])
            pet_mmh    = float(row["potential_evaporation"])
            m.set_value("atmosphere_water__time_integral_of_precipitation_mass_flux",
                        precip_mmh / 1000.0)
            m.set_value("water_potential_evaporation_flux",
                        pet_mmh / 1000.0 / 3600.0)
            m.update()

        # ── extract final states ─────────────────────────────────────────────
        states = {
            "cat_id":             cat_id,
            "spinup_end":         SPINUP_WY26_END,
            "soil_storage_m":     float(m.soil_reservoir["storage_m"]),
            "gw_storage_m":       float(m.gw_reservoir["storage_m"]),
            "nash_storage_0":     float(m.nash_storage[0]),
            "nash_storage_1":     float(m.nash_storage[1]),
            "runoff_queue":       [float(x) for x in m.runoff_queue_m_per_timestep],
        }

        with open(out_json, "w") as f:
            json.dump(states, f, indent=2)

        dt = time.time() - t0
        sm_frac = states["soil_storage_m"] / float(m.soil_reservoir["storage_max_m"])
        gw_frac = states["gw_storage_m"]   / float(m.gw_reservoir["storage_max_m"])
        return cat_id, f"ok  soil={sm_frac:.2%}  gw={gw_frac:.2%}  {dt:.0f}s"

    except Exception as e:
        return cat_id, f"ERROR: {e}\n{traceback.format_exc()}"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="1-year deterministic CFE spin-up for all Alabama catchments.")
    parser.add_argument("--nwm-wy2026-dir", required=True,
                        help="NWM operational WY2026 forcing dir "
                             "(e.g. /mnt/disk2/al_storm_2026/nwm_operational/wy2026/forcings)")
    parser.add_argument("--params-dir",   default=enkf.DEFAULT_PARAMS_DIR)
    parser.add_argument("--config-file",  default=enkf.DEFAULT_CONFIG_FILE)
    parser.add_argument("--states-dir",   required=True,
                        help="output dir for <cat>_spinup_states.json files")
    parser.add_argument("--cat-id",       default="all",
                        help="single cat-XXXXXX or 'all'")
    parser.add_argument("--procs",        type=int, default=20)
    parser.add_argument("--force",        action="store_true",
                        help="rerun even if states JSON already exists")
    args = parser.parse_args()

    # ── import bmi_cfe (must be inside main after sys.path set) ──────────────
    import importlib
    cfe_dir = Path(enkf.DEFAULT_CFE_DIR)
    if str(cfe_dir) not in sys.path:
        sys.path.insert(0, str(cfe_dir))
    enkf.bmi_cfe = importlib.import_module("bmi_cfe")

    Path(args.states_dir).mkdir(parents=True, exist_ok=True)

    # ── cat list ───────────────────────────────────────────────────────────────
    if args.cat_id == "all":
        params_dir = Path(args.params_dir)
        cats = sorted([d.name for d in params_dir.iterdir()
                       if d.is_dir() and d.name.startswith("cat-")]
                      if any(params_dir.iterdir()) else [])
        if not cats:
            cats = sorted([f.stem.replace("_best_params", "")
                           for f in params_dir.glob("cat-*_best_params.json")])
        if not cats and ALL_CATS_FILE.exists():
            cats = [l.strip() for l in ALL_CATS_FILE.read_text().splitlines() if l.strip()]
        if not cats:
            print("ERROR: could not enumerate catchments. Use --cat-id cat-XXXXXX or "
                  f"create {ALL_CATS_FILE}")
            sys.exit(1)
    else:
        cats = [args.cat_id]

    already = sum(1 for c in cats
                  if (Path(args.states_dir) / f"{c}_spinup_states.json").exists())
    print(f"Spin-up: {len(cats)} catchments  |  {already} already done  "
          f"|  {args.procs} procs")

    G["args"] = args

    t_all = time.time()
    errors = []
    if args.procs == 1:
        for cat in cats:
            cid, msg = run_spinup_catchment(cat)
            print(f"  {cid}: {msg}")
            if msg.startswith("ERROR"):
                errors.append(cid)
    else:
        with Pool(args.procs) as pool:
            for cid, msg in pool.imap_unordered(run_spinup_catchment, cats):
                print(f"  {cid}: {msg}")
                if msg.startswith("ERROR"):
                    errors.append(cid)

    elapsed = time.time() - t_all
    print(f"\nDone: {len(cats)} catchments in {elapsed/60:.1f} min")
    if errors:
        print(f"ERRORS ({len(errors)}): {', '.join(errors[:10])}")
        sys.exit(1)
    print(f"States saved to: {args.states_dir}")
    print("Next: re-run run_cfe_enkf_forecast_hrrr_multihour_al_storm.py "
          "with --spinup-states-dir pointing here.")


if __name__ == "__main__":
    main()
