"""Core extension logic: ramp state fields and scale biofuel fractions."""

import numpy as np
from . import config as cfg


def extend_states(values_2100, rates_2100, years):
    """Extend land-use state fields from 2100 to *years*.

    Strategy:
    - Linearly ramp each cell's rate of change to zero by YR_AFOLU_ZERO.
    - Clamp declining variables at zero; redirect excess to secdf.
    - Hold constant from YR_AFOLU_ZERO onward.

    Parameters
    ----------
    values_2100 : dict[str, ndarray]  — 2-D (lat, lon) arrays at 2100.
    rates_2100  : dict[str, ndarray]  — per-cell rates at 2100.
    years       : array-like of output years (first should be 2101).

    Returns
    -------
    dict[str, ndarray] — 3-D arrays (time, lat, lon) for each state var.
    """
    years = np.asarray(years)
    nt = len(years)
    shape2d = values_2100["primf"].shape
    n_ramp = cfg.YR_AFOLU_ZERO - cfg.YR_END_INPUT  # 49

    # Initialise output arrays
    out = {v: np.empty((nt,) + shape2d, dtype=np.float32) for v in cfg.STATE_VARS}

    # Previous timestep values — start from 2100
    prev = {v: values_2100[v].copy().astype(np.float64) for v in cfg.STATE_VARS}

    for ti, yr in enumerate(years):
        dt = yr - cfg.YR_END_INPUT
        # Rate multiplier: linear ramp from 1→0 over n_ramp years
        mult = max(0.0, 1.0 - dt / n_ramp)

        # Apply rates (except secdf which absorbs residuals)
        delta = {}
        for v in cfg.STATE_VARS:
            if v == "secdf":
                continue
            delta[v] = rates_2100[v] * mult

        # Compute tentative new values and clamp at zero
        excess = np.zeros(shape2d, dtype=np.float64)
        new = {}
        for v in cfg.STATE_VARS:
            if v == "secdf":
                continue
            tentative = prev[v] + delta[v]
            # Where value was valid (not NaN) and goes negative, clamp
            neg = np.isfinite(tentative) & (tentative < 0)
            excess += np.where(neg, -tentative, 0.0)
            new[v] = np.where(neg, 0.0, tentative)

        # secdf absorbs its own rate plus any excess from clamped variables
        secdf_tentative = prev["secdf"] + rates_2100["secdf"] * mult + excess
        new["secdf"] = np.where(np.isfinite(secdf_tentative), secdf_tentative,
                                prev["secdf"])

        # Store and advance
        for v in cfg.STATE_VARS:
            out[v][ti] = new[v].astype(np.float32)
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
