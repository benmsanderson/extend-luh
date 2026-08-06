#!/usr/bin/env python3
"""Surgically regenerate ONLY the wood-harvest ramp files for given scenarios.

Replicates the wood-harvest block of src.pipeline.extend_scenario
(pipeline.py:435-471) so the delivered VL/H files use the current
linear ramp-to-zero logic without regenerating every other output.
"""
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, ".")
from src import config as cfg


def regen(key: str) -> None:
    sc = cfg.get_scenario(key)

    EXT_YEARS = np.arange(cfg.YR_END_INPUT + 1, cfg.YR_END_OUTPUT + 1)

    ds_wh_in = xr.open_dataset(sc.woodharvest_path)
    wh_2100 = ds_wh_in["woodharvest"].isel(time=-1).values
    ccodes = ds_wh_in.coords["country_code"].values
    ds_wh_in.close()

    # Linear ramp from 1 at 2100 to 0 at YR_AFOLU_ZERO (common endpoint
    # across scenarios), matching src.pipeline.extend_scenario.
    yr_zero = cfg.YR_AFOLU_ZERO
    wh_ramp = np.clip((yr_zero - EXT_YEARS) / (yr_zero - cfg.YR_END_INPUT), 0, 1)
    wh_ext = wh_ramp[:, np.newaxis] * wh_2100[np.newaxis, :]

    global_attrs = {
        "Conventions": "CF-1.6",
        "activity_id": "ScenarioMIP",
        "source_model": sc.model,
        "source_scenario": sc.scenario,
        "institution": "CICERO Center for International Climate Research",
        "source_states_file": sc.states_file,
        "source_csv_file": cfg.CSV_FILE.name,
        "creation_date": datetime.utcnow().isoformat() + "Z",
        "grid_label": "gn",
        "frequency": "yr",
        "nominal_resolution": "50 km",
    }

    ds_wh_out = xr.Dataset(
        {"woodharvest": xr.DataArray(
            wh_ext.astype(np.float32),
            dims=["time", "country_code"],
            coords={"time": EXT_YEARS, "country_code": ccodes},
            attrs={"units": "MgC",
                   "long_name": "wood harvest carbon demand (ramp × IAM 2100)"},
        )},
        attrs={**global_attrs,
               "title": f"Wood harvest demand ramp ({EXT_YEARS[0]}-{EXT_YEARS[-1]})",
               "description": ("IAM 2100 wood harvest demand scaled by a "
                               "linear ramp from 1 (at 2100) to 0 (at "
                               f"{yr_zero}). Intended for use "
                               "as: h(t) = max(h_maintenance(t), this_file).")},
    )

    WH_STEM = Path(sc.woodharvest_file).stem
    fname = f"{WH_STEM}_extension_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"
    out = sc.output_dir / fname

    # Remove any older-named variant (e.g. the _extension_ramp_ files)
    for old in sc.output_dir.glob(f"{WH_STEM}_extension*_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"):
        if old.name != fname:
            old.unlink()
            print(f"[{key}] removed {old.name}")

    ds_wh_out.to_netcdf(out)
    print(f"[{key}] wrote {out}  (yr_zero={yr_zero}, "
          f"ramp {len(EXT_YEARS)} yr x {len(ccodes)} countries)")


if __name__ == "__main__":
    keys = sys.argv[1:] or ["VL", "H"]
    for k in keys:
        regen(k)
