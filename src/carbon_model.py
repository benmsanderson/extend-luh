"""
Carbon cycle model for AFOLU flux calibration and projection.

Model D: 4-predictor regression with cohort-based stock-change dynamics:
    AFOLU(t) = β + γ·t + α_trans·F_trans(t) + α_stock·F_stock(t, τ)

where:
    - β: baseline AFOLU emissions (Mt CO₂/yr)
    - γ: time trend coefficient (Mt CO₂/yr per decade)
    - α_trans: transition flux multiplier (dimensionless)
    - α_stock: stock-change flux multiplier (dimensionless)
    - τ: cohort relaxation timescale (years)
"""

import numpy as np
from . import config as cfg

# ============================================================================
# Carbon cycle parameters
# ============================================================================

# Carbon densities: tC per km² (for transition flux)
# Includes entries for all variable names across scenarios.
CARBON_DENSITY = {
    'primf': 15_000, 'secdf':  8_000,
    'primn':  1_500, 'secdn':  1_000,
    'c3ann':    500, 'c3nfx':    500, 'c3per':    500,
    'c4ann':    500, 'c4per':    500,
    'pastr':    800, 'range':    600,
    'urban':    200, 'pltns':  6_000,
    # Aggregated types (used by scenarios without primary/secondary split)
    'forest': 10_000,  # area-weighted blend of primf+secdf+pltns
    # IAM-specific aliases
    'timber': 6_000,   # IMAGE name for plantation forests (= pltns)
}

# Equilibrium carbon stock for regrowing land types (tC per km²)
# Represents the total carbon a regrowing stand will accumulate over its lifetime
C_EQ = {
    'primf':      0,  # mature forest, already at equilibrium
    'secdf': 12_000,  # relaxes toward primary forest stock
    'primn':      0,
    'secdn':    800,  # recovering non-forest
    'c3ann':      0, 'c3nfx':  0, 'c3per':  0,
    'c4ann':      0, 'c4per':  0,
    'pastr':      0, 'range':  0,
    'urban':      0,
    'pltns':  8_000,  # managed plantation regrowth
    # Aggregated
    'forest': 10_000,  # net regrowth in aggregate forest pool
    # IAM-specific aliases
    'timber': 8_000,   # IMAGE name for plantation forests (= pltns)
}


def regrow_vars_for(state_vars=None):
    """Return the list of regrowing variables (C_EQ > 0) for a state var list."""
    if state_vars is None:
        state_vars = cfg.STATE_VARS
    return [v for v in state_vars if C_EQ.get(v, 0) > 0]


# Default REGROW_VARS for backward compatibility
REGROW_VARS = regrow_vars_for(cfg.STATE_VARS)

# Conversion factor: tonne carbon to Megatonne CO₂
TC_TO_MTCO2 = 3.667 / 1e6


# ============================================================================
# Core model functions
# ============================================================================

def stock_flux_for_tau(tau, area_incr, c_eq=None, regrow_vars=None):
    """
    Compute stock-change flux for a given relaxation timescale.

    Uses cohort convolution with exponential approach to equilibrium:
        G_v(t) = G_v(t-1) × exp(-1/τ) + ΔA_v(t)
        flux(t) = -Σ_v (C_eq_v / τ) × G_v(t) × TC_TO_MTCO2

    Parameters
    ----------
    tau : float
        Relaxation timescale (years). Time for 63% approach to equilibrium.
    area_incr : dict of arrays
        Net area increment per year for each regrowing land type (km²).
        Keys are land-use types, values are 1-D arrays of length n_years.
    c_eq : dict, optional
        Equilibrium carbon stock (tC/km²) for each land type.
        Defaults to module-level C_EQ.
    regrow_vars : list[str], optional
        Which variables to include. Default: keys of area_incr.

    Returns
    -------
    flux : ndarray
        Stock-change flux timeseries (Mt CO₂/yr). Negative = removal.
    """
    if c_eq is None:
        c_eq = C_EQ
    if regrow_vars is None:
        regrow_vars = list(area_incr.keys())

    n = len(next(iter(area_incr.values())))
    decay = np.exp(-1.0 / tau)
    flux = np.zeros(n)

    for v in regrow_vars:
        ceq = c_eq.get(v, 0)
        if ceq == 0:
            continue
        G = 0.0
        incr = area_incr[v]
        rate = ceq / tau * TC_TO_MTCO2
        for t in range(n):
            G = G * decay + incr[t]
            flux[t] += -rate * G

    return flux


def calibrate_afolu_model(transition_flux, area_incr, afolu_target, years,
                          tau_grid=None, tau_fixed=None, regrow_vars=None):
    """
    Calibrate 4-predictor AFOLU model via profile likelihood over τ.

    Model: AFOLU = β + γ·t + α_trans·F_trans + α_stock·F_stock(τ)

    Parameters
    ----------
    transition_flux : ndarray
        Transition flux timeseries (Mt CO₂/yr), length n.
    area_incr : dict of arrays
        Net area increment per year for each regrowing type (km²).
    afolu_target : ndarray
        Target AFOLU flux timeseries (Mt CO₂/yr), length n.
    years : ndarray
        Year labels for the timeseries (for time trend).
    tau_grid : ndarray, optional
        Grid of τ values to scan (years). Default: np.arange(20, 201, 2).
    tau_fixed : float, optional
        If provided, skip profile scan and use this τ value.

    Returns
    -------
    result : dict
        Calibration results containing:
        - 'beta', 'gamma', 'alpha_trans', 'alpha_stock': fitted coefficients
        - 'tau': τ value used (from profile optimum or tau_fixed)
        - 'tau_opt': profile optimum τ (if profile scan was run)
        - 'r2': coefficient of determination
        - 'rmse': root mean squared error
        - 'fitted': fitted AFOLU flux timeseries
        - 'residual': residual timeseries
        - 'r2_profile', 'tau_profile': profile likelihood arrays (if scanned)
    """
    n = len(afolu_target)
    time_trend = (years - years.mean()) / 10.0  # decades, centered

    if tau_fixed is not None:
        # Use fixed τ without profile scan
        tau_use = tau_fixed
        tau_opt = None
        r2_profile = None
        tau_profile = None
    else:
        # Profile likelihood scan
        if tau_grid is None:
            tau_grid = np.arange(20, 201, 2)

        r2_profile = np.zeros(len(tau_grid))
        coefs_profile = np.zeros((len(tau_grid), 4))

        for i, tau_try in enumerate(tau_grid):
            sf = stock_flux_for_tau(tau_try, area_incr, regrow_vars=regrow_vars)
            X = np.column_stack([np.ones(n), time_trend, transition_flux, sf])
            c, _, _, _ = np.linalg.lstsq(X, afolu_target, rcond=None)
            pred = X @ c
            ss_res = np.sum((afolu_target - pred) ** 2)
            ss_tot = np.sum((afolu_target - afolu_target.mean()) ** 2)
            r2_profile[i] = 1 - ss_res / ss_tot
            coefs_profile[i] = c

        best_idx = np.argmax(r2_profile)
        tau_opt = tau_grid[best_idx]
        tau_use = tau_opt
        tau_profile = tau_grid

    # Final fit at chosen τ
    stock_flux = stock_flux_for_tau(tau_use, area_incr, regrow_vars=regrow_vars)
    X = np.column_stack([np.ones(n), time_trend, transition_flux, stock_flux])
    coefs, _, _, _ = np.linalg.lstsq(X, afolu_target, rcond=None)
    beta, gamma, alpha_trans, alpha_stock = coefs
    fitted = X @ coefs
    residual = afolu_target - fitted

    ss_res = np.sum(residual ** 2)
    ss_tot = np.sum((afolu_target - afolu_target.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    rmse = np.sqrt(ss_res / n)

    return {
        'beta': beta,
        'gamma': gamma,
        'alpha_trans': alpha_trans,
        'alpha_stock': alpha_stock,
        'tau': tau_use,
        'tau_opt': tau_opt,
        'r2': r2,
        'rmse': rmse,
        'fitted': fitted,
        'residual': residual,
        'time_trend': time_trend,
        'stock_flux': stock_flux,
        'r2_profile': r2_profile,
        'tau_profile': tau_profile,
    }


def forward_solve_ramp(beta, gamma, alpha_trans, alpha_stock, tau,
                       area_incr_pre, F_trans_unit, dA_unit,
                       afolu_target_ext, years_pre, years_ext,
                       regrow_vars=None, calibration_residual_2100=0.0,
                       blend_timescale=10.0):
    """
    Forward solve for AFOLU-consistent ramp r(t) from 2101 onward.

    Given calibrated coefficients and unit-rate fluxes at 2100, solve for
    r(t) such that the implied AFOLU flux matches the IAM target.

    Parameters
    ----------
    beta, gamma, alpha_trans, alpha_stock : float
        Calibrated model coefficients.
    tau : float
        Cohort relaxation timescale (years).
    area_incr_pre : dict of arrays
        Pre-2100 area increments (km²) for computing initial G state.
    F_trans_unit : float
        Transition flux at unit rate in 2100 (Mt CO₂/yr).
    dA_unit : dict
        Area increment at unit rate in 2100 (km²) for each regrowing type.
    afolu_target_ext : ndarray
        Target AFOLU flux for extension period (Mt CO₂/yr), length n_ext.
    years_pre : ndarray
        Pre-2100 years (for computing time trend at 2100).
    years_ext : ndarray
        Extension years (2101–2500).
    calibration_residual_2100 : float, optional
        Residual (calibrated - IAM) at year 2100. Used to blend smoothly
        from the calibrated model value to the IAM target. Default 0.
    blend_timescale : float, optional
        E-folding timescale (years) for blending the calibration residual
        into the post-2100 target. Default 10.

    Returns
    -------
    result : dict
        - 'ramp': r(t) timeseries (0–1)
        - 'afolu_reconstructed': reconstructed AFOLU flux from r(t)
        - 'residual': afolu_reconstructed - afolu_target_ext
        - 'baseline': frozen baseline (β + γ·t₂₁₀₀)
        - 'G_end_pre': cohort state G_v at end of pre-2100 period
    """
    if regrow_vars is None:
        regrow_vars = REGROW_VARS

    n_pre = len(years_pre)
    n_ext = len(years_ext)
    decay = np.exp(-1.0 / tau)

    # Compute G_v at end of pre-2100 period
    G_end_pre = {}
    for v in regrow_vars:
        G = 0.0
        for t in range(n_pre):
            G = G * decay + area_incr_pre[v][t]
        G_end_pre[v] = G

    # Stock-change flux from one unit of new area
    S_new_unit = sum(
        -(C_EQ[v] / tau) * dA_unit[v] * TC_TO_MTCO2
        for v in regrow_vars if C_EQ.get(v, 0) > 0
    )

    # Baseline frozen at 2100 time trend value
    time_2100 = (2100 - years_pre.mean()) / 10.0
    baseline = beta + gamma * time_2100

    # Denominator for forward solve
    denom = baseline + alpha_trans * F_trans_unit + alpha_stock * S_new_unit

    # Blended target: smooth transition from calibrated model to IAM target
    # target_blended(t) = IAM(t) + residual_2100 × exp(-(t-2100)/T_blend)
    yr0 = int(years_ext[0])
    blend_correction = np.array([
        calibration_residual_2100 * np.exp(-(yr0 + t - cfg.YR_END_INPUT) / blend_timescale)
        for t in range(n_ext)
    ])

    # Forward solve
    r_derived = np.zeros(n_ext)
    afolu_reconstructed = np.zeros(n_ext)
    G_work = {v: G_end_pre[v] for v in regrow_vars}

    for t in range(n_ext):
        target_t = afolu_target_ext[t] + blend_correction[t]

        # Stock-change flux from previous cohorts (before this year's increment)
        S_prev = sum(
            -(C_EQ[v] / tau) * G_work[v] * decay * TC_TO_MTCO2
            for v in regrow_vars if C_EQ.get(v, 0) > 0
        )

        # Solve for r(t)
        r_t = (target_t - alpha_stock * S_prev) / denom
        r_t = np.clip(r_t, 0.0, 1.0)
        r_derived[t] = r_t

        # Update G_v
        for v in regrow_vars:
            G_work[v] = G_work[v] * decay + r_t * dA_unit[v]

        # Reconstruct AFOLU flux
        S_full = sum(
            -(C_EQ[v] / tau) * G_work[v] * TC_TO_MTCO2
            for v in regrow_vars if C_EQ.get(v, 0) > 0
        )
        afolu_reconstructed[t] = (baseline * r_t +
                                  alpha_trans * F_trans_unit * r_t +
                                  alpha_stock * S_full)

    residual = afolu_reconstructed - afolu_target_ext

    return {
        'ramp': r_derived,
        'afolu_reconstructed': afolu_reconstructed,
        'residual': residual,
        'baseline': baseline,
        'G_end_pre': G_end_pre,
        'blend_correction': blend_correction,
    }
