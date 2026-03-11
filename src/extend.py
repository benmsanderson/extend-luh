"""Core extension logic: ramp state fields and scale biofuel fractions."""

import numpy as np
from . import config as cfg


def extend_states(values_2100, rates_2100, years, ramp=None,
                  state_vars=None, residual_var=None):
    """Extend land-use state fields from 2100 to *years*.

    Strategy:
    - Multiply each cell's rate of change by the ramp multiplier at each year.
    - Clamp declining variables at zero; redirect excess to residual_var.
    - Hold constant once ramp reaches zero.

    Parameters
    ----------
    values_2100 : dict[str, ndarray]  — 2-D (lat, lon) arrays at 2100.
    rates_2100  : dict[str, ndarray]  — per-cell rates at 2100.
    years       : array-like of output years (first should be 2101).
    ramp        : array-like, optional — per-year multipliers, same length as
                  *years*. If None, uses a linear ramp (1 at 2101 → 0 at
                  YR_AFOLU_ZERO).
    state_vars  : list[str], optional — variables to extend.
                  Default: cfg.STATE_VARS.
    residual_var : str, optional — variable used as residual for per-cell
                   conservation. Default: "secdf".

    Returns
    -------
    dict[str, ndarray] — 3-D arrays (time, lat, lon) for each state var.
    """
    if state_vars is None:
        state_vars = cfg.STATE_VARS
    if residual_var is None:
        residual_var = "secdf"

    years = np.asarray(years)
    nt = len(years)
    first_var = state_vars[0]
    shape2d = values_2100[first_var].shape
    n_ramp = cfg.YR_AFOLU_ZERO - cfg.YR_END_INPUT  # 49

    # Build or validate ramp
    if ramp is not None:
        ramp = np.asarray(ramp, dtype=np.float64)
        assert len(ramp) == nt, f"ramp length {len(ramp)} != years length {nt}"
    else:
        # Default: linear ramp
        ramp = np.clip(1.0 - (years - cfg.YR_END_INPUT) / n_ramp, 0.0, 1.0)

    # Work and store entirely in float64 to avoid float32 rounding artefacts
    # in global sums.  Convert to float32 only when writing final NetCDF.
    out = {v: np.empty((nt,) + shape2d, dtype=np.float64) for v in state_vars}

    # Previous timestep values — start from 2100
    prev = {v: values_2100[v].copy().astype(np.float64) for v in state_vars}

    # Per-cell total to preserve: sum of all state vars at 2100.
    # On ocean cells (all NaN) this is 0, but residual guard below keeps them NaN.
    cell_total = np.zeros(shape2d, dtype=np.float64)
    for v in state_vars:
        cell_total += np.where(np.isfinite(values_2100[v]), values_2100[v], 0.0)

    # Land mask: cells where residual_var (and hence all states) is finite
    land = np.isfinite(values_2100[residual_var])

    # Normalise rates so they sum to zero per cell.  Input datasets may have
    # small per-cell imbalances (especially with aggregated variables like
    # "forest"); without this correction, the residual variable would absorb
    # the full imbalance, causing a rate discontinuity at the 2100 boundary.
    rates = {v: rates_2100[v].astype(np.float64) for v in state_vars}
    rate_sum = np.zeros(shape2d, dtype=np.float64)
    for v in state_vars:
        rate_sum += np.where(np.isfinite(rates[v]), rates[v], 0.0)
    # Spread the imbalance onto the residual variable (largest pool)
    rates[residual_var] = np.where(land,
                                   rates[residual_var] - rate_sum,
                                   rates[residual_var])

    for ti, yr in enumerate(years):
        mult = ramp[ti]

        # Apply ramped rates to all variables except the residual
        new = {}
        for v in state_vars:
            if v == residual_var:
                continue
            tentative = prev[v] + rates[v] * mult
            new[v] = np.where(tentative < 0, 0.0, tentative)

        # Residual absorbs only the small clamping deficit, not rate imbalance
        non_residual_sum = np.zeros(shape2d, dtype=np.float64)
        for v in state_vars:
            if v != residual_var:
                non_residual_sum += np.where(np.isfinite(new[v]), new[v], 0.0)
        residual_new = cell_total - non_residual_sum
        new[residual_var] = np.where(land, residual_new, np.nan)

        # Store and advance
        for v in state_vars:
            out[v][ti] = new[v]
            prev[v] = new[v]

    return out


def extend_biofuel(crpbiof_2100, beccs_factors, years):
    """Scale crpbiof by BECCS factor for each year, capped at 1.0.

    Parameters
    ----------
    crpbiof_2100 : ndarray (lat, lon)
    beccs_factors : Series or dict mapping year→factor
    years : array-like

    Returns
    -------
    ndarray (time, lat, lon)
    """
    years = np.asarray(years)
    nt = len(years)
    out = np.empty((nt,) + crpbiof_2100.shape, dtype=np.float32)
    for ti, yr in enumerate(years):
        f = beccs_factors.get(int(yr), beccs_factors.get(yr, 1.0))
        scaled = crpbiof_2100 * f
        out[ti] = np.minimum(scaled, 1.0).astype(np.float32)
    return out


def extend_constant(field_2100, n_years):
    """Repeat a 2-D field for n_years (for hold-constant variables)."""
    return np.broadcast_to(field_2100[np.newaxis], (n_years,) + field_2100.shape)
