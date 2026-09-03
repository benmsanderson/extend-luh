"""
Batch pipeline: calibration, gridded extension, and verification.

Extracts the core computation from notebooks 02, 03, and 04 into callable
functions that can be looped over scenarios by run_all.py.
"""

import gc
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from . import carbon_model as cm
from . import config as cfg
from .config import ScenarioConfig
from .data import (
    beccs_scaling_factors,
    filter_scenario,
    get_variable,
    load_biof,
    load_csv,
    load_derived_ramp,
    load_management,
    load_states,
    rates_2100,
    state_2100,
)
from .extend import extend_biofuel, extend_states


# ── Shared helpers ──────────────────────────────────────────────────────

def _cell_area(lat, lon):
    """Grid-cell area (km²) from latitude / longitude vectors."""
    R_EARTH = 6371.0
    dlat = np.abs(np.diff(lat[:2]))[0]
    dlon = np.abs(np.diff(lon[:2]))[0]
    lat_rad = np.deg2rad(lat)
    area_1d = (R_EARTH ** 2) * np.deg2rad(dlat) * np.deg2rad(dlon) * np.cos(lat_rad)
    return np.broadcast_to(area_1d[:, np.newaxis], (len(lat), len(lon)))


# =====================================================================
# Step 1 — Calibration  (≈ notebook 02)
# =====================================================================

def calibrate_scenario(sc: ScenarioConfig, *, verbose: bool = True) -> dict:
    """Run AFOLU calibration for one scenario.

    Returns a dict summarising calibration quality and the path to the
    saved .npz file.
    """
    log = print if verbose else (lambda *a, **k: None)
    sc.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load IAM target ─────────────────────────────────────────────
    df = load_csv()
    sc_data = filter_scenario(df, scen=sc)
    afolu_target = get_variable(sc_data, "Emissions|CO2|AFOLU")

    # ── Load gridded data ───────────────────────────────────────────
    try:
        ds_input = load_states(scen=sc)
        has_gridded = True
    except FileNotFoundError:
        log(f"  ⚠ {sc.key}: no gridded states file – skipping")
        return {"key": sc.key, "status": "no_gridded_data"}

    lat = ds_input.lat.values
    lon = ds_input.lon.values
    cell_area_2d = _cell_area(lat, lon)

    REGROW_VARS = cm.regrow_vars_for(sc.state_vars)
    TC = cm.TC_TO_MTCO2

    # ── Pre-2100 transition flux & area increments ──────────────────
    input_years = ds_input.time.values.astype(int)
    input_flux_years = input_years[1:]
    n_flux = len(input_flux_years)

    F_trans_pre = np.zeros(n_flux)
    area_incr_pre = {v: np.zeros(n_flux) for v in REGROW_VARS}

    for ti in range(1, len(input_years)):
        flux_t = 0.0
        for v in sc.state_vars:
            delta = ds_input[v].isel(time=ti).values - ds_input[v].isel(time=ti - 1).values
            flux_t += np.nansum(-delta * cm.CARBON_DENSITY[v] * cell_area_2d) * TC
        F_trans_pre[ti - 1] = flux_t

        for v in REGROW_VARS:
            delta = ds_input[v].isel(time=ti).values - ds_input[v].isel(time=ti - 1).values
            area_incr_pre[v][ti - 1] = np.nansum(delta * cell_area_2d)

    years_pre = input_flux_years
    afolu_cal = afolu_target.loc[input_flux_years].values

    log(f"  {sc.key}: F_trans {input_flux_years[0]}–{input_flux_years[-1]}, "
        f"range [{F_trans_pre.min():.0f}, {F_trans_pre.max():.0f}]")

    # ── Calibrate ───────────────────────────────────────────────────
    cal_kwargs = dict(
        transition_flux=F_trans_pre,
        area_incr=area_incr_pre,
        afolu_target=afolu_cal,
        years=input_flux_years,
        tau_grid=np.arange(20, 201, 2),
        regrow_vars=REGROW_VARS,
    )
    if sc.tau_fixed is not None:
        cal_kwargs["tau_fixed"] = sc.tau_fixed

    result = cm.calibrate_afolu_model(**cal_kwargs)

    beta = result["beta"]
    gamma = result["gamma"]
    alpha_trans = result["alpha_trans"]
    alpha_stock = result["alpha_stock"]
    tau_opt = result["tau"]
    r2 = result["r2"]
    fitted = result["fitted"]
    model_type = "4-predictor"

    # Physical consistency: fall back to 3-predictor if α_stock < 0
    if alpha_stock < 0:
        log(f"  {sc.key}: α_stock < 0 ({alpha_stock:.3f}), refitting 3-predictor")
        n_cal = len(afolu_cal)
        time_trend_cal = (input_flux_years - input_flux_years.mean()) / 10.0
        X3 = np.column_stack([np.ones(n_cal), time_trend_cal, F_trans_pre])
        coefs3, _, _, _ = np.linalg.lstsq(X3, afolu_cal, rcond=None)
        beta, gamma, alpha_trans = coefs3
        alpha_stock = 0.0
        fitted = X3 @ coefs3
        ss_res = np.sum((afolu_cal - fitted) ** 2)
        ss_tot = np.sum((afolu_cal - afolu_cal.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot
        model_type = "3-predictor"

    rmse_cal = np.sqrt(np.mean((afolu_cal - fitted) ** 2))

    # ── Rate normalisation ──────────────────────────────────────────
    vals_2100 = state_2100(ds_input, scen=sc)
    rates_2100_dict = rates_2100(ds_input, scen=sc)

    rate_sum = np.zeros_like(rates_2100_dict[sc.state_vars[0]], dtype=np.float64)
    for v in sc.state_vars:
        rate_sum += np.where(np.isfinite(rates_2100_dict[v]), rates_2100_dict[v], 0.0)
    land_mask = np.isfinite(vals_2100[sc.residual_var])
    rates_2100_dict[sc.residual_var] = np.where(
        land_mask,
        rates_2100_dict[sc.residual_var] - rate_sum,
        rates_2100_dict[sc.residual_var],
    )

    # Unit-rate fluxes at 2100
    F_trans_unit = sum(
        np.nansum(-rates_2100_dict[v] * cm.CARBON_DENSITY[v] * cell_area_2d) * TC
        for v in sc.state_vars
    )
    dA_unit = {}
    for v in REGROW_VARS:
        dA_unit[v] = np.nansum(rates_2100_dict[v] * cell_area_2d)

    # ── Forward solve for ramp ──────────────────────────────────────
    years_ext = np.arange(cfg.YR_END_INPUT + 1, cfg.YR_END_OUTPUT + 1)
    afolu_target_ext = afolu_target.loc[years_ext].values

    cal_residual_2100 = fitted[-1] - afolu_cal[-1]

    ramp_result = cm.forward_solve_ramp(
        beta=beta, gamma=gamma,
        alpha_trans=alpha_trans, alpha_stock=alpha_stock, tau=tau_opt,
        area_incr_pre=area_incr_pre,
        F_trans_unit=F_trans_unit,
        dA_unit=dA_unit,
        afolu_target_ext=afolu_target_ext,
        years_pre=years_pre,
        years_ext=years_ext,
        regrow_vars=REGROW_VARS,
        calibration_residual_2100=cal_residual_2100,
    )

    ramp = ramp_result["ramp"]

    # ── Save ────────────────────────────────────────────────────────
    out_path = sc.output_dir / "afolu_calibration.npz"
    np.savez(
        out_path,
        years=years_ext,
        ramp=ramp,
        tau=tau_opt,
        beta=beta,
        gamma=gamma,
        alpha_trans=alpha_trans,
        alpha_stock=alpha_stock,
        r2=r2,
        afolu_target=afolu_target_ext,
        afolu_reconstructed=ramp_result["afolu_reconstructed"],
        residual=ramp_result["residual"],
    )

    ds_input.close()

    rmse_ramp = np.sqrt(np.mean(ramp_result["residual"] ** 2))
    ramp_end_yr = int(years_ext[np.where(ramp < 1e-6)[0][0]]) if np.any(ramp < 1e-6) else None

    log(f"  {sc.key}: R²={r2:.4f}, τ={tau_opt:.0f}, "
        f"RMSE_cal={rmse_cal:.0f}, RMSE_ramp={rmse_ramp:.0f}, "
        f"model={model_type}")

    return {
        "key": sc.key,
        "status": "ok",
        "model_type": model_type,
        "tau": float(tau_opt),
        "r2": float(r2),
        "rmse_cal": float(rmse_cal),
        "rmse_ramp": float(rmse_ramp),
        "r_2101": float(ramp[0]),
        "ramp_zero_yr": ramp_end_yr,
        "alpha_trans": float(alpha_trans),
        "alpha_stock": float(alpha_stock),
        "beta": float(beta),
        "gamma": float(gamma),
        "cal_residual_2100": float(cal_residual_2100),
        "calibration_file": str(out_path),
    }


# =====================================================================
# Step 2 — Gridded extension  (≈ notebook 03)
# =====================================================================

def extend_scenario(sc: ScenarioConfig, *, verbose: bool = True) -> dict:
    """Produce gridded extension NetCDF files for one scenario.

    Requires calibrate_scenario() to have been run first (reads the
    saved afolu_calibration.npz).
    """
    log = print if verbose else (lambda *a, **k: None)
    sc.output_dir.mkdir(parents=True, exist_ok=True)

    # Remove old extension NetCDFs to avoid stale files with different dates
    for old in sc.output_dir.glob("*_extension_*.nc"):
        old.unlink()
        log(f"  {sc.key}: removed stale {old.name}")

    EXT_YEARS = np.arange(cfg.YR_END_INPUT + 1, cfg.YR_END_OUTPUT + 1)

    # ── Load AFOLU-consistent ramp ──────────────────────────────────
    derived_ramp_full = load_derived_ramp(EXT_YEARS, scen=sc)
    ramp_active = derived_ramp_full > 1e-8
    last_active = np.where(ramp_active)[0][-1] if np.any(ramp_active) else 0
    RAMP_YEARS = EXT_YEARS[: last_active + 1]
    derived_ramp = derived_ramp_full[: len(RAMP_YEARS)]

    log(f"  {sc.key}: ramp {RAMP_YEARS[0]}–{RAMP_YEARS[-1]} "
        f"({len(RAMP_YEARS)} yr), r=[{derived_ramp.min():.4f}, {derived_ramp.max():.4f}]")

    # ── Load gridded inputs ─────────────────────────────────────────
    ds_states = load_states(scen=sc)
    ds_mgmt = load_management(scen=sc)
    ds_biof = load_biof(scen=sc)

    lat = ds_states.lat.values
    lon = ds_states.lon.values
    nlat, nlon = len(lat), len(lon)

    vals = state_2100(ds_states, scen=sc)
    rates = rates_2100(ds_states, scen=sc)
    crpbiof_2100 = ds_biof["crpbiof"].isel(time=-1).values

    # ── Extend states ───────────────────────────────────────────────
    ext_states = extend_states(
        vals, rates, RAMP_YEARS, ramp=derived_ramp,
        state_vars=sc.state_vars, residual_var=sc.residual_var,
    )

    # ── Extend biofuel ──────────────────────────────────────────────
    df = load_csv()
    sc_data = filter_scenario(df, scen=sc)
    beccs = beccs_scaling_factors(sc_data)
    beccs_dict = beccs.to_dict()
    crpbiof_ext = extend_biofuel(crpbiof_2100, beccs_dict, EXT_YEARS)

    # ── Extend management ───────────────────────────────────────────
    # Management/biofuel variables may live in the management file, the
    # biofuel file, or (for scenarios with embedded management) the states
    # file, depending on the IAM.  Fetch each variable's 2100 slice from
    # whichever input contains it, preferring the management dataset so the
    # existing per-scenario sources are unchanged.
    def _get_2100(name):
        for ds in (ds_mgmt, ds_biof, ds_states):
            if name in ds:
                return ds[name].isel(time=-1).values
        return None

    ext_mgmt = {}
    for v in sc.mgmt_hold_constant:
        arr = _get_2100(v)
        if arr is not None:
            ext_mgmt[v] = arr
    for v in sc.mgmt_biofuel_vars:
        if v == "crpbiof":
            continue
        arr = _get_2100(v)
        if arr is not None:
            ext_mgmt[v] = extend_biofuel(arr, beccs_dict, EXT_YEARS)

    # ── Build & write NetCDF files ──────────────────────────────────
    _ALL_VAR_META = {
        "primf": {"long_name": "forested primary land", "units": "share of carea"},
        "secdf": {"long_name": "potentially forested secondary land", "units": "share of carea"},
        "primn": {"long_name": "non-forested primary land", "units": "share of carea"},
        "secdn": {"long_name": "potentially non-forested secondary land", "units": "share of carea"},
        "forest": {"long_name": "forest land", "units": "share of carea"},
        "c3ann": {"long_name": "C3 annual crops", "units": "share of carea"},
        "c3nfx": {"long_name": "C3 nitrogen-fixing crops", "units": "share of carea"},
        "c3per": {"long_name": "C3 perennial crops", "units": "share of carea"},
        "c4ann": {"long_name": "C4 annual crops", "units": "share of carea"},
        "c4per": {"long_name": "C4 perennial crops", "units": "share of carea"},
        "pastr": {"long_name": "managed pasture", "units": "share of carea"},
        "range": {"long_name": "rangeland", "units": "share of carea"},
        "urban": {"long_name": "urban land", "units": "share of carea"},
        "pltns": {"long_name": "plantation forests", "units": "share of carea"},
        "crpbiof": {"long_name": "biofuel fraction of cropland", "units": "fraction of cropland area"},
        "crpbf_c3per": {"long_name": "C3 perennial crops grown as second-generation biofuels", "units": "share of c3per"},
        "crpbf_c4per": {"long_name": "C4 perennial crops grown as second-generation biofuels", "units": "share of c4per"},
        "rndwd": {"long_name": "industrial roundwood fraction of wood harvest biomass carbon", "units": "share of bioh"},
        "fulwd": {"long_name": "fuelwood fraction of wood harvest biomass carbon", "units": "share of bioh"},
        "pltns_wdprd": {"long_name": "fraction of harvested plantations biomass used for wood products", "units": "share of pltns"},
        "pltns_bfuel": {"long_name": "fraction of harvested plantations biomass used for bioenergy", "units": "share of pltns"},
    }
    for crop in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]:
        _ALL_VAR_META[f"irrig_{crop}"] = {"long_name": f"irrigated fraction of {crop}", "units": f"share of {crop}"}
        _ALL_VAR_META[f"fertl_{crop}"] = {"long_name": f"fertilization rate for {crop}", "units": "kg ha-1 yr-1"}
        _ALL_VAR_META[f"cpbf1_{crop}"] = {"long_name": f"first-generation biofuel fraction of {crop}", "units": f"share of {crop}"}
    _ALL_VAR_META["flood"] = {"long_name": "flooded fraction of rice cropland", "units": "fraction of rice area"}

    def _make_da(varname, data):
        meta = _ALL_VAR_META.get(varname, {})
        return xr.DataArray(
            data.astype(np.float32),
            dims=["time", "lat", "lon"],
            attrs={"_FillValue": np.float32(1e20), "cell_methods": "time:mean",
                   "standard_name": "area_fraction", **meta},
        )

    global_attrs = {
        "Conventions": "CF-1.6",
        "activity_id": "ScenarioMIP",
        "source_model": sc.model,
        "source_scenario": sc.scenario,
        "institution": "CICERO Center for International Climate Research",
        "extension_method": (
            f"Per-cell rates from 2100 × AFOLU-consistent ramp r(t). "
            f"Negatives clamped, deficit absorbed by {sc.residual_var}. "
            f"Biofuel scaled by BECCS ratio. Management held constant."
        ),
        "source_states_file": sc.states_file,
        "source_csv_file": cfg.CSV_FILE.name,
        "creation_date": datetime.utcnow().isoformat() + "Z",
        "grid_label": "gn",
        "frequency": "yr",
        "nominal_resolution": "50 km",
    }
    if sc.mgmt_file is not None:
        global_attrs["source_management_file"] = sc.mgmt_file
    if sc.fertl_file is not None:
        global_attrs["source_fertl_file"] = sc.fertl_file
    if sc.flood_file is not None:
        global_attrs["source_flood_file"] = sc.flood_file
    if sc.var_renames:
        global_attrs["variable_renames"] = str(sc.var_renames)

    STATES_STEM = Path(sc.states_file).stem
    BIOF_STEM = Path(sc.biof_file).stem
    if sc.mgmt_file:
        MGMT_STEM = Path(sc.mgmt_file).stem
    elif sc.fertl_file:
        MGMT_STEM = Path(sc.fertl_file).stem
    else:
        MGMT_STEM = f"mgmt_{STATES_STEM}"

    # File 1 — states (transient, active ramp)
    state_dvars = {}
    for v in sc.state_vars:
        out_name = sc.var_renames.get(v, v)
        state_dvars[out_name] = _make_da(out_name, ext_states[v])
    ds_states_out = xr.Dataset(
        state_dvars,
        coords={"time": RAMP_YEARS.astype(np.float64), "lat": lat, "lon": lon},
    )
    ds_states_out.attrs = {
        **global_attrs,
        "title": f"LUH state extension ({RAMP_YEARS[0]}-{RAMP_YEARS[-1]})",
    }
    fname_states = f"{STATES_STEM}_extension_transient_{RAMP_YEARS[0]}-{RAMP_YEARS[-1]}.nc"
    del ext_states
    gc.collect()

    # File 2 — biofuel (transient, full 2101-2500)
    ds_biof_out = xr.Dataset(
        {"crpbiof": _make_da("crpbiof", crpbiof_ext)},
        coords={"time": EXT_YEARS.astype(np.float64), "lat": lat, "lon": lon},
    )
    ds_biof_out.attrs = {
        **global_attrs,
        "title": f"LUH biofuel extension ({EXT_YEARS[0]}-{EXT_YEARS[-1]})",
    }
    fname_biof = f"{BIOF_STEM}_extension_transient_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"

    # File 3 — management transient (BECCS-scaled vars)
    mgmt_trans_dvars = {}
    for v in ("crpbf_c3per", "crpbf_c4per"):
        if v in ext_mgmt:
            mgmt_trans_dvars[v] = _make_da(v, np.asarray(ext_mgmt[v]))
    ds_mgmt_trans_out = xr.Dataset(
        mgmt_trans_dvars,
        coords={"time": EXT_YEARS.astype(np.float64), "lat": lat, "lon": lon},
    )
    ds_mgmt_trans_out.attrs = {
        **global_attrs,
        "title": f"LUH management extension – transient ({EXT_YEARS[0]}-{EXT_YEARS[-1]})",
    }
    fname_mgmt_trans = f"{MGMT_STEM}_extension_transient_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"

    # File 4 — management static (hold-constant)
    mgmt_static_dvars = {}
    for v in sorted(sc.mgmt_hold_constant):
        if v not in ext_mgmt:
            continue
        mgmt_static_dvars[v] = _make_da(v, ext_mgmt[v][np.newaxis])
    ds_mgmt_static_out = xr.Dataset(
        mgmt_static_dvars,
        coords={"time": EXT_YEARS[0:1].astype(np.float64), "lat": lat, "lon": lon},
    )
    ds_mgmt_static_out.attrs = {
        **global_attrs,
        "title": "LUH management extension – static (constant at 2100)",
    }
    fname_mgmt_static = f"{MGMT_STEM}_extension_static.nc"

    # File 5 — wood harvest demand ramp (country-level)
    ds_wh_out = None
    fname_wh = None
    ds_fw_out = None
    fname_fw = None
    try:
        ds_wh_in = xr.open_dataset(sc.woodharvest_path)
        wh_2100 = ds_wh_in["woodharvest"].isel(time=-1).values  # (country_code,)
        ccodes = ds_wh_in.coords["country_code"].values
        ds_wh_in.close()

        # Linear ramp from 1 at 2100 to 0 at YR_AFOLU_ZERO, the year the
        # AFOLU emissions target reaches zero.  We deliberately pin the
        # endpoint to YR_AFOLU_ZERO (common across scenarios) rather than
        # to the forward-solved ramp's numerical zero-crossing
        # (RAMP_YEARS[-1]), which lands at a scenario-dependent year.
        # Unlike the AFOLU ramp (which may start < 1 due to boundary-blend
        # correction), wood harvest demand must be continuous at 2100.
        yr_zero = cfg.YR_AFOLU_ZERO                # common endpoint (r==0)
        wh_ramp = np.clip((yr_zero - EXT_YEARS) / (yr_zero - cfg.YR_END_INPUT), 0, 1)
        wh_ext = wh_ramp[:, np.newaxis] * wh_2100[np.newaxis, :]

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
        fname_wh = f"{WH_STEM}_extension_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"
        log(f"  {sc.key}: built wood harvest ramp ({len(EXT_YEARS)} yr × {len(ccodes)} countries)")
    except Exception as exc:
        log(f"  {sc.key}: skipping wood harvest extension ({exc})")

    try:
        ds_fw_in = xr.open_dataset(sc.fuelwood_path)
        fw_2100 = ds_fw_in["fuelwood"].isel(time=-1).values  # (country_code,)
        fw_ccodes = ds_fw_in.coords["country_code"].values
        ds_fw_in.close()

        # Hold fuelwood fraction constant at 2100 value
        fw_ext = np.broadcast_to(fw_2100[np.newaxis, :], (len(EXT_YEARS), len(fw_ccodes)))

        ds_fw_out = xr.Dataset(
            {"fuelwood": xr.DataArray(
                fw_ext.astype(np.float32),
                dims=["time", "country_code"],
                coords={"time": EXT_YEARS, "country_code": fw_ccodes},
                attrs={"units": "fraction",
                       "long_name": "fuelwood fraction (constant at 2100 value)"},
            )},
            attrs={**global_attrs,
                   "title": f"Fuelwood fraction extension ({EXT_YEARS[0]}-{EXT_YEARS[-1]})",
                   "description": "Fuelwood fraction held constant at 2100 values."},
        )
        FW_STEM = Path(sc.fuelwood_file).stem
        fname_fw = f"{FW_STEM}_extension_static_{EXT_YEARS[0]}-{EXT_YEARS[-1]}.nc"
        log(f"  {sc.key}: built fuelwood extension (constant at 2100)")
    except Exception as exc:
        log(f"  {sc.key}: skipping fuelwood extension ({exc})")

    # Write
    output_plan = [(ds_states_out, fname_states), (ds_biof_out, fname_biof)]
    if mgmt_trans_dvars:
        output_plan.append((ds_mgmt_trans_out, fname_mgmt_trans))
    if mgmt_static_dvars:
        output_plan.append((ds_mgmt_static_out, fname_mgmt_static))
    if ds_wh_out is not None:
        output_plan.append((ds_wh_out, fname_wh))
    if ds_fw_out is not None:
        output_plan.append((ds_fw_out, fname_fw))

    comp = {"zlib": True, "complevel": 4, "dtype": "float32"}
    total_bytes = 0
    written_files = []
    for ds_out, fname in output_plan:
        path = sc.output_dir / fname
        enc = {}
        for v in ds_out.data_vars:
            if "lat" in ds_out[v].dims:
                enc[v] = {**comp, "chunksizes": (1, nlat, nlon)}
            else:
                enc[v] = {**comp}
        ds_out.to_netcdf(path, encoding=enc)
        sz = os.path.getsize(path)
        total_bytes += sz
        written_files.append(fname)
        log(f"  {sc.key}: wrote {fname} ({sz / 1e6:.1f} MB)")
        ds_out.close()

    ds_states.close()
    ds_mgmt.close()
    ds_biof.close()
    gc.collect()

    return {
        "key": sc.key,
        "status": "ok",
        "files": written_files,
        "total_mb": total_bytes / 1e6,
        "ramp_years": f"{RAMP_YEARS[0]}-{RAMP_YEARS[-1]}",
        "n_ramp_years": len(RAMP_YEARS),
    }


# =====================================================================
# Step 3 — Verification  (≈ notebook 04)
# =====================================================================

def verify_scenario(sc: ScenarioConfig, *, verbose: bool = True) -> dict:
    """Recompute AFOLU fluxes from gridded extension and compare to IAM target.

    Returns verification metrics (RMSE, jump at 2100 boundary, etc.).
    """
    import glob

    log = print if verbose else (lambda *a, **k: None)

    REGROW_VARS = cm.regrow_vars_for(sc.state_vars)
    TC = cm.TC_TO_MTCO2

    # ── Load calibration params ─────────────────────────────────────
    cal = np.load(sc.output_dir / "afolu_calibration.npz")
    tau_opt = float(cal["tau"])
    beta = float(cal["beta"])
    gamma = float(cal["gamma"])
    alpha_trans = float(cal["alpha_trans"])
    alpha_stock = float(cal["alpha_stock"])
    r_derived = cal["ramp"]
    ext_years_saved = cal["years"]

    # ── IAM target ──────────────────────────────────────────────────
    df = load_csv()
    sc_data = filter_scenario(df, scen=sc)
    afolu_target = get_variable(sc_data, "Emissions|CO2|AFOLU")

    # ── Find extension file ─────────────────────────────────────────
    ext_files = sorted(glob.glob(str(sc.output_dir / "*_extension_transient_2101-*.nc")))
    states_ext_file = [
        f for f in ext_files
        if "states" in f.lower() or "step4" in f.lower() or "timeseries" in f.lower()
    ]
    if not states_ext_file:
        return {"key": sc.key, "status": "no_extension_file"}
    states_ext_file = states_ext_file[0]

    ds_ext = xr.open_dataset(states_ext_file)
    ds_input = load_states(scen=sc)

    lat = ds_ext.lat.values
    lon = ds_ext.lon.values
    cell_area_2d = _cell_area(lat, lon)
    ramp_years = ds_ext.time.values.astype(int)

    # ── Pre-2100 fluxes ─────────────────────────────────────────────
    input_years = ds_input.time.values.astype(int)
    input_flux_years = input_years[1:]
    n_pre = len(input_flux_years)

    time_trend_pre = (input_flux_years - input_flux_years.mean()) / 10.0
    time_2100 = (2100 - input_flux_years.mean()) / 10.0

    transition_flux_pre = np.zeros(n_pre)
    area_incr_pre = {v: np.zeros(n_pre) for v in REGROW_VARS}

    for ti in range(1, len(input_years)):
        af = 0.0
        for v in sc.state_vars:
            delta = ds_input[v].isel(time=ti).values - ds_input[v].isel(time=ti - 1).values
            af += np.nansum(-delta * cm.CARBON_DENSITY[v] * cell_area_2d) * TC
        transition_flux_pre[ti - 1] = af
        for v in REGROW_VARS:
            delta = ds_input[v].isel(time=ti).values - ds_input[v].isel(time=ti - 1).values
            area_incr_pre[v][ti - 1] = np.nansum(delta * cell_area_2d)

    # ── Extension period fluxes ─────────────────────────────────────
    vals_2100 = {v: ds_input[v].isel(time=-1).values for v in sc.state_vars}
    n_ramp = len(ramp_years)
    full_ext_years = np.arange(cfg.YR_END_INPUT + 1, cfg.YR_END_OUTPUT + 1)
    n_ext = len(full_ext_years)

    ext_transition = np.zeros(n_ext)
    area_incr_ext = {v: np.zeros(n_ext) for v in REGROW_VARS}

    for ti in range(n_ramp):
        af = 0.0
        for v in sc.state_vars:
            ext_v = sc.var_renames.get(v, v)
            prev = vals_2100[v] if ti == 0 else ds_ext[ext_v].isel(time=ti - 1).values
            curr = ds_ext[ext_v].isel(time=ti).values
            delta = curr - prev
            af += np.nansum(-delta * cm.CARBON_DENSITY[v] * cell_area_2d) * TC
        ext_transition[ti] = af
        for v in REGROW_VARS:
            ext_v = sc.var_renames.get(v, v)
            prev = vals_2100[v] if ti == 0 else ds_ext[ext_v].isel(time=ti - 1).values
            curr = ds_ext[ext_v].isel(time=ti).values
            area_incr_ext[v][ti] = np.nansum((curr - prev) * cell_area_2d)

    # ── Stock-change flux ───────────────────────────────────────────
    decay = np.exp(-1.0 / tau_opt)
    G_state = {v: 0.0 for v in REGROW_VARS}

    stock_flux_pre = np.zeros(n_pre)
    for t in range(n_pre):
        for v in REGROW_VARS:
            ceq = cm.C_EQ[v]
            if ceq == 0:
                continue
            G_state[v] = G_state[v] * decay + area_incr_pre[v][t]
            stock_flux_pre[t] += -(ceq / tau_opt) * G_state[v] * TC

    ext_stock = np.zeros(n_ext)
    for t in range(n_ext):
        for v in REGROW_VARS:
            ceq = cm.C_EQ[v]
            if ceq == 0:
                continue
            G_state[v] = G_state[v] * decay + area_incr_ext[v][t]
            ext_stock[t] += -(ceq / tau_opt) * G_state[v] * TC

    # ── AFOLU from gridded data ─────────────────────────────────────
    ramp_for_ext = np.zeros(n_ext)
    ramp_for_ext[: len(r_derived)] = r_derived[:n_ext]

    baseline = beta + gamma * time_2100
    calibrated_pre = (
        beta + gamma * time_trend_pre
        + alpha_trans * transition_flux_pre
        + alpha_stock * stock_flux_pre
    )
    ext_calib = (
        baseline * ramp_for_ext
        + alpha_trans * ext_transition
        + alpha_stock * ext_stock
    )

    # ── Metrics ─────────────────────────────────────────────────────
    iam_ext_vals = np.array([
        afolu_target.loc[yr] if yr in afolu_target.index else 0.0
        for yr in full_ext_years
    ])
    residual = ext_calib - iam_ext_vals

    # Jump at 2100 boundary
    cal_at_2100 = calibrated_pre[-1]
    ext_at_2101 = ext_calib[0]
    jump = ext_at_2101 - cal_at_2100

    # RMSE over ramp period
    ramp_mask = ramp_for_ext > 1e-8
    if np.any(ramp_mask):
        rmse_ext = np.sqrt(np.mean(residual[ramp_mask] ** 2))
        max_resid = np.max(np.abs(residual[ramp_mask]))
    else:
        rmse_ext = 0.0
        max_resid = 0.0

    ds_ext.close()
    ds_input.close()

    log(f"  {sc.key}: RMSE={rmse_ext:.1f}, jump={jump:+.1f}, "
        f"max|resid|={max_resid:.1f}")

    return {
        "key": sc.key,
        "status": "ok",
        "rmse": float(rmse_ext),
        "jump_2100": float(jump),
        "max_residual": float(max_resid),
        "cal_at_2100": float(cal_at_2100),
        "ext_at_2101": float(ext_at_2101),
        "iam_at_2100": float(afolu_target.loc[2100]),
        "iam_at_2101": float(afolu_target.loc[2101]),
    }
