"""Load and filter CSV / NetCDF data for the target scenario."""

import pandas as pd
import xarray as xr
import numpy as np
from . import config as cfg


# ---------- CSV helpers ----------

def load_csv():
    """Load the full emissions timeseries CSV."""
    return pd.read_csv(cfg.CSV_FILE)


def filter_scenario(df, model=cfg.MODEL, scenario=cfg.SCENARIO):
    """Return rows for a single model+scenario."""
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


# ---------- NetCDF helpers ----------

def load_states():
    """Open the states dataset (lazy)."""
    return xr.open_dataset(cfg.STATES_FILE)


def load_management():
    """Open the management dataset (lazy)."""
    return xr.open_dataset(cfg.MGMT_FILE)


def load_biof():
    """Open the biofuel dataset."""
    return xr.open_dataset(cfg.BIOF_FILE)


def state_2100(ds):
    """Return all state variables at 2100 as a dict of 2-D arrays."""
    return {v: ds[v].isel(time=-1).values for v in cfg.STATE_VARS}


def rates_2100(ds):
    """Per-cell rate of change at 2100 (2099→2100 difference)."""
    return {v: ds[v].isel(time=-1).values - ds[v].isel(time=-2).values
            for v in cfg.STATE_VARS}
