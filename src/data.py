"""Load and filter CSV / NetCDF data for the target scenario."""

import pandas as pd
import xarray as xr
import numpy as np
from . import config as cfg
from .config import ScenarioConfig


def _resolve_scenario(scenario):
    """Return a ScenarioConfig from a key string, ScenarioConfig, or None."""
    if scenario is None:
        return cfg.SCENARIOS[cfg.DEFAULT_SCENARIO_KEY]
    if isinstance(scenario, str):
        return cfg.get_scenario(scenario)
    return scenario


def _normalize_coords(ds):
    """Rename latitude/longitude → lat/lon if the short forms are absent."""
    renames = {}
    if 'latitude' in ds.dims and 'lat' not in ds.dims:
        renames['latitude'] = 'lat'
    if 'longitude' in ds.dims and 'lon' not in ds.dims:
        renames['longitude'] = 'lon'
    return ds.rename(renames) if renames else ds


# ---------- CSV helpers ----------

def load_csv():
    """Load the full emissions timeseries CSV."""
    return pd.read_csv(cfg.CSV_FILE)


def filter_scenario(df, model=None, scenario=None, scen=None):
    """Return rows for a single model+scenario.

    Parameters
    ----------
    df : DataFrame
    model, scenario : str, optional
        Explicit model/scenario names. If omitted, taken from *scen*.
    scen : str or ScenarioConfig, optional
        Scenario key or config. Default: VL.
    """
    sc = _resolve_scenario(scen)
    if model is None:
        model = sc.model
    if scenario is None:
        scenario = sc.scenario
    return df[(df["model"] == model) & (df["scenario"] == scenario)]


def get_variable(df, variable, region="World"):
    """Extract a single timeseries as a pandas Series indexed by year."""
    row = df[(df["variable"] == variable) & (df["region"] == region)]
    if len(row) == 0:
        raise KeyError(f"No data for {variable}, region={region}")
    yr_cols = [c for c in df.columns if c not in
               ("model", "scenario", "region", "workflow", "variable", "unit")]
    ts = row.iloc[0][yr_cols].astype(float)
    ts.index = ts.index.astype(float).astype(int)
    ts.name = variable
    return ts


# ---------- BECCS scaling curve ----------

def beccs_scaling_factors(df):
    """Compute global BECCS scaling factor timeseries (relative to 2100).

    Returns a Series indexed by integer year (2100–2500), values = ratio.
    """
    beccs = df[df["variable"] == "Emissions|CO2|BECCS"]
    yr_cols = [c for c in df.columns if c not in
               ("model", "scenario", "region", "workflow", "variable", "unit")]
    total = beccs[yr_cols].sum(axis=0).astype(float)
    total.index = total.index.astype(float).astype(int)
    val_2100 = total[2100]
    factors = total / val_2100
    return factors.loc[2100:]


# ---------- AFOLU ramp multiplier ----------

def afolu_ramp(years=None):
    """Return the AFOLU rate multiplier: 1 at 2100, linear to 0 at YR_AFOLU_ZERO.

    Parameters
    ----------
    years : array-like, optional
        Years to evaluate. Default: 2100–2500.

    Returns
    -------
    np.ndarray of multipliers, same length as *years*.
    """
    if years is None:
        years = np.arange(cfg.YR_END_INPUT, cfg.YR_END_OUTPUT + 1)
    years = np.asarray(years, dtype=float)
    n = cfg.YR_AFOLU_ZERO - cfg.YR_END_INPUT  # 49
    dt = years - cfg.YR_END_INPUT
    mult = np.clip(1.0 - dt / n, 0.0, 1.0)
    return mult


def load_derived_ramp(years=None, scen=None):
    """Load the AFOLU-consistent ramp derived by notebook 02.

    Parameters
    ----------
    years : array-like, optional
        Years to return. Must be a subset of the saved years (2101-2500).
        Default: all saved years.
    scen : str or ScenarioConfig, optional
        Scenario key or config. Default: VL.

    Returns
    -------
    np.ndarray of multipliers, same length as *years*.
    """
    sc = _resolve_scenario(scen)
    ramp_path = sc.output_dir / 'afolu_calibration.npz'
    data = np.load(ramp_path)
    saved_years = data['years']
    saved_ramp = data['ramp']
    if years is None:
        return saved_ramp
    years = np.asarray(years, dtype=int)
    indices = np.searchsorted(saved_years, years)
    return saved_ramp[indices]


# ---------- NetCDF helpers ----------

def load_states(scen=None):
    """Open the states dataset (lazy), normalising lat/lon coord names."""
    sc = _resolve_scenario(scen)
    return _normalize_coords(xr.open_dataset(sc.states_path))


def load_management(scen=None):
    """Open the management dataset (lazy).

    Three cases:
    - mgmt_in_states=True (H/GCAM): management vars live in the states file.
    - mgmt_file set (VL/MAgPIE): single combined management NetCDF.
    - component files only (M/IMAGE): fertl_file + flood_file merged on the fly.
      Note: close() on the merged result does not release the underlying file
      handles — they are freed when the Dataset goes out of scope.
    """
    sc = _resolve_scenario(scen)
    if sc.mgmt_in_states:
        return xr.open_dataset(sc.states_path)
    if sc.mgmt_path is not None:
        return xr.open_dataset(sc.mgmt_path)
    # Multi-file scenario (e.g. M/IMAGE): merge component management files
    parts = []
    if sc.fertl_path is not None:
        parts.append(xr.open_dataset(sc.fertl_path))
    if sc.flood_path is not None:
        parts.append(xr.open_dataset(sc.flood_path))
    if not parts:
        raise FileNotFoundError(
            f"Scenario '{sc.key}' has no management file and no component "
            f"files (fertl_file, flood_file).")
    merged = xr.merge(parts) if len(parts) > 1 else parts[0]
    return _normalize_coords(merged)


def load_biof(scen=None):
    """Open the biofuel dataset."""
    sc = _resolve_scenario(scen)
    return xr.open_dataset(sc.biof_path)


def state_2100(ds, scen=None):
    """Return all state variables at 2100 as a dict of 2-D arrays."""
    sc = _resolve_scenario(scen)
    return {v: ds[v].isel(time=-1).values for v in sc.state_vars}


def rates_2100(ds, scen=None):
    """Per-cell rate of change at 2100 (2099→2100 difference)."""
    sc = _resolve_scenario(scen)
    return {v: ds[v].isel(time=-1).values - ds[v].isel(time=-2).values
            for v in sc.state_vars}
