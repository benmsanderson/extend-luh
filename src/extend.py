"""Core extension logic: ramp state fields and scale biofuel fractions."""

import numpy as np
from . import config as cfg


def extend_states(values_2100, rates_2100, years, ramp=None,
                  state_vars=None, residual_var=None, renormalize_overfull=True,
                  overshoot_mode="scale"):
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

    # Land mask: a cell is land if *any* state variable is finite there.
    # We must not key this off a single variable (e.g. the residual): some
    # IAMs (IMAGE) store NaN — meaning "zero" — for sparse classes such as
    # secondary forest over genuine land cells, so keying off secdf alone
    # would drop most of the land and leave those cells un-conserved.
    land = np.zeros(shape2d, dtype=bool)
    for v in state_vars:
        land |= np.isfinite(values_2100[v])

    # On land, treat a NaN state (or rate) as 0 — "class absent here" — so the
    # residual can absorb and per-cell conservation holds.  Off land, keep NaN.
    def _fill(arr):
        a = np.where(np.isfinite(arr), arr, 0.0).astype(np.float64)
        return np.where(land, a, np.nan)

    # Previous timestep values — start from 2100 (NaN→0 on land)
    prev = {v: _fill(values_2100[v]) for v in state_vars}

    # Per-cell total to preserve: sum of all state vars at 2100.
    # On ocean cells (all NaN) this is 0, but residual guard below keeps them NaN.
    cell_total = np.zeros(shape2d, dtype=np.float64)
    for v in state_vars:
        cell_total += np.where(np.isfinite(values_2100[v]), values_2100[v], 0.0)

    # Safeguard: a few input cells have a 2100 state sum > 1 (e.g. IMAGE cells
    # where primary land and cropland overlap to >100% of the cell).  Fractions
    # cannot exceed the cell, so renormalise those cells to sum to 1 at 2100;
    # the extension then preserves a physical total.  Cells with sum <= 1 are
    # left untouched (their unaccounted fraction, e.g. inland water, is kept).
    over_unit = land & (cell_total > 1.0)
    if renormalize_overfull and over_unit.any():
        for v in state_vars:
            prev[v] = np.where(over_unit, prev[v] / cell_total, prev[v])
        cell_total = np.where(over_unit, 1.0, cell_total)

    # Normalise rates so they sum to zero per cell.  Input datasets may have
    # small per-cell imbalances (especially with aggregated variables like
    # "forest"); without this correction, the residual variable would absorb
    # the full imbalance, causing a rate discontinuity at the 2100 boundary.
    rates = {v: _fill(rates_2100[v]) for v in state_vars}
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

        # If non-residual vars overshot the available cell area (which can happen
        # when crops have positive rates right up to 2100 and the residual pool is
        # nearly exhausted), scale all non-residual values down proportionally so
        # the residual never goes negative.  Without this, the residual absorbs
        # the excess as a negative value, producing individual crop fractions > 1.
        overshoot = land & (non_residual_sum > cell_total) & (non_residual_sum > 0)
        if overshoot.any() and overshoot_mode == "limit_expansion":
            # Only the expanding classes are limited, to the room left after
            # this step's declines; non-expanding classes keep their prescribed
            # path instead of being squeezed to make room.
            base = {v: np.minimum(new[v], prev[v]) for v in new}
            incr = {v: np.maximum(new[v] - prev[v], 0.0) for v in new}
            base_sum = sum(np.nan_to_num(base[v]) for v in new)
            incr_sum = sum(np.nan_to_num(incr[v]) for v in new)
            room = np.maximum(cell_total - base_sum, 0.0)
            frac = np.where(incr_sum > 0, np.minimum(1.0, room / np.where(incr_sum > 0, incr_sum, 1.0)), 1.0)
            for v in new:
                new[v] = np.where(overshoot, base[v] + incr[v] * frac, new[v])
            non_residual_sum = np.where(overshoot, base_sum + incr_sum * frac, non_residual_sum)
        elif overshoot.any():
            safe_sum = np.where(non_residual_sum > 0, non_residual_sum, 1.0)
            scale = np.where(overshoot, cell_total / safe_sum, 1.0)
            for v in state_vars:
                if v != residual_var:
                    new[v] = new[v] * scale
            non_residual_sum = np.where(overshoot, cell_total, non_residual_sum)

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


def extend_harvest(harvest_2100, ext_years):
    """Hold country-level wood harvest constant at its 2100 level for all extension years.

    Mirrors the CMIP6 LUH2 extension convention (Hurtt et al. 2020), which froze
    the end-of-period national wood harvest demand for all subsequent years.

    Parameters
    ----------
    harvest_2100 : ndarray (n_countries,)
        Per-country harvest at 2100 (e.g. MgC/yr).
    ext_years : ndarray (n_ext,)
        Full extension year range (e.g. 2101–2500).

    Returns
    -------
    ndarray (n_ext, n_countries)
    """
    n_ext = len(ext_years)
    return np.broadcast_to(harvest_2100[np.newaxis], (n_ext,) + harvest_2100.shape).copy()
