# Extending Land-Use Harmonization Beyond 2100

This document describes the methodology for extending gridded land-use projections beyond the end of IAM input (2100) out to 2500, in a way that preserves consistency with the IAM's AFOLU CO₂ trajectory.

---

## Pipeline Overview

The workflow has four stages, each implemented as a notebook:

| Stage | Notebook | Purpose |
|-------|----------|---------|
| **1. Explore** | `01_explore` | Visualise IAM trajectories, inspect LUH states and rates at 2100, prototype the extension logic |
| **2. Calibrate** | `02_afolu_calibration` | Fit a carbon cycle model to pre-2100 data and invert it to derive an AFOLU-consistent rate ramp $r(t)$ |
| **3. Extend** | `03_gridded_extension` | Apply $r(t)$ to gridded states; scale biofuel by BECCS; hold management constant; write NetCDF output |
| **4. Verify** | `04_verify_afolu` | Re-derive AFOLU fluxes from the gridded output and confirm they match the IAM target |

Library code lives in `src/`: `config.py` (paths, years, variable lists), `data.py` (I/O helpers), `extend.py` (per-cell extension logic), and `carbon_model.py` (calibration and forward solve).

---

## The Extension Problem

The IAM provides gridded land-use states and transition rates up to 2100, plus a global AFOLU CO₂ trajectory out to 2500. We need to:

1. Continue the 13 gridded state fractions beyond 2100, ramping rates smoothly to zero.
2. Ensure the implied AFOLU emissions match the IAM's target trajectory — a naive linear ramp does not achieve this because it ignores carbon cycle inertia (regrowing forests continue to absorb CO₂ long after planting).
3. Handle management layers (wood harvest, irrigation, fertiliser, biofuel) appropriately.

The solution is to **calibrate a carbon cycle model** on pre-2100 data, **invert it** to derive the rate multiplier $r(t)$ that reproduces the IAM AFOLU curve, and then apply that ramp to the gridded extension.

---

## Model Structure: Three Flux Channels

The total AFOLU CO₂ flux is decomposed into three physically-motivated components:

### **A. Transition Flux** — Instantaneous Carbon Release/Uptake

When land changes category (e.g., forest → cropland), the carbon density difference is released or sequestered immediately:

$$F^{\text{trans}}(t) = -\sum_v \Delta A_v(t) \cdot C_v$$

where:
- $\Delta A_v(t)$: area change of land type $v$ (km2/yr)
- $C_v$: carbon density of land type $v$ (tC/km²)

**Proportional to**: Rate of land-use change ($r(t)$ in the extension)

**Example**: Deforestation releases ~15,000 tC/km² - 500 tC/km² = 14,500 tC/km² per unit area converted.

---

### **B. Stock-Change Flux** — Cohort-Based Exponential Relaxation

Regrowing forests accumulate carbon following an exponential approach to equilibrium with timescale $\tau$:

$$\frac{dS}{dt} = \frac{C^{\text{eq}} - S(t)}{\tau}$$

For a stand reforested in year $t_r$, the uptake rate decays exponentially with stand age $a = t - t_r$:

$$\text{uptake}(a) = \frac{C^{\text{eq}}}{\tau} \cdot e^{-a/\tau}$$

**Key insight**: We track each year's area increment as a **cohort**. The total stock-change flux is computed efficiently via a running sum:

$$G_v(t) = G_v(t-1) \cdot e^{-1/\tau} + \Delta A_v(t)$$

$$F^{\text{stock}}(t) = -\sum_v \frac{C^{\text{eq}}_v}{\tau} \cdot G_v(t)$$

This flux has two components:
- **Committed removals** from pre-2100 cohorts (decays exponentially, not controlled by $r(t)$)
- **New removals** from post-2100 cohorts (proportional to $r(t)$)

**Timescale $\tau$**: Fixed at 35 years (90% recovery in ~81 years), consistent with literature values for tropical secondary forest regrowth (Poorter et al. 2016).

---

### **C. Time Trend** — Secular Change in AFOLU Emissions

A linear time trend $\gamma \cdot t$ captures the secular increase in carbon removal capacity driven by **expanding reforestation area** over time. This is distinct from individual cohort aging.

$$F^{\text{trend}}(t) = \gamma \cdot t$$

where $t$ is measured in decades from the mean year of the calibration period.

**Why is this needed?** Without the time trend, the calibration inflates $\tau$ to ~200 years (unphysically slow) to explain the expanding-area signal via the stock-change channel. With the time trend absorbing this signal, $\tau$ remains at a physically defensible value (~35 years).

**In the forward solve** (post-2100): The time trend is frozen at its 2100 value and scales with $r(t)$, so it naturally ramps to zero when land-use change stops. This prevents the model from attributing the expanding-area signal to unphysically slow cohort aging (τ → 200 years).

---

## Mathematical Formulation

The complete model regresses AFOLU flux onto four predictors:

$$\boxed{\text{AFOLU}(t) = \beta + \gamma \cdot t + \alpha_{\text{trans}} \cdot F^{\text{trans}}(t) + \alpha_{\text{stock}} \cdot F^{\text{stock}}(t;\, \tau)}$$

**Parameters**:
- $\beta$: Baseline AFOLU emissions (Mt CO₂/yr)
- $\gamma$: Time trend coefficient (Mt CO₂/yr per decade)
- $\alpha_{\text{trans}}$: Transition flux multiplier (dimensionless)
- $\alpha_{\text{stock}}$: Stock-change flux multiplier (dimensionless)
- $\tau$: Cohort relaxation timescale (typically fixed at 35 years, literature range 20–50 years)

**Calibration**: Ordinary least squares (OLS) regression on pre-2100 IAM data, with optional profile likelihood scan over $\tau$ to verify the optimum is within the literature range.

### Key Assumption: Global (Zero-Dimensional) Model

The carbon cycle model is fit to **globally-aggregated** quantities: transition fluxes and area increments are summed over all grid cells before regression. The resulting ramp $r(t)$ is a single scalar applied uniformly to every grid cell at each timestep.

This means:
- The **spatial pattern** of land-use change is fixed (frozen at the 2100 pattern); only the global amplitude varies over time.
- Carbon density, equilibrium stocks, and the relaxation timescale $\tau$ are treated as globally-representative constants rather than spatially-varying fields.
- Regional differences in recovery dynamics (e.g. tropical forests recover faster than boreal forests) are not resolved.

This is a deliberate simplification. A future refinement could fit separate parameters by biome or latitude band — for example, tropical and boreal forests likely have substantially different $\tau$ and $C^{\text{eq}}$ — and derive spatially-varying ramps. The current global approach is adequate for producing a self-consistent extension but may under- or over-estimate fluxes in specific regions.

---

## Forward Solve: Deriving the AFOLU-Consistent Ramp

Post-2100, we **invert the calibrated model** to solve for $r(t)$ that makes the implied AFOLU flux match the IAM target.

### Key Insight: Proportionality

All controllable terms are proportional to $r(t)$:
- Baseline: $(β + γ · t_{2100})$ frozen at 2100, scales with $r(t)$
- Transition flux: $r(t) \times F^{\text{trans}}_{\text{unit}}$ (per-unit-rate flux at 2100)
- Stock-change from new cohorts: $r(t) \times \Delta A_{\text{unit}}$ (per-unit-rate area increment)

Only the committed stock-change from pre-2100 cohorts decays independently of $r(t)$.

### Solution at Each Year $t > 2100$

$$r(t) = \frac{\text{AFOLU}_{\text{target}}(t) - \alpha_{\text{stock}} \cdot S_{\text{prev}}(t)}{\text{baseline} + \alpha_{\text{trans}} \cdot F^{\text{trans}}_{\text{unit}} + \alpha_{\text{stock}} \cdot S_{\text{new,unit}}}$$

where:
- $\text{baseline} = \beta + \gamma \cdot t_{2100}$ (frozen at 2100)
- $S_{\text{prev}}(t)$: stock-change flux from pre-2100 cohorts (committed removals, computed by decaying $G_v$ states forward without new area)
- $S_{\text{new,unit}}$: stock-change flux from one unit of new area at 2100

**Sequential solve**: Each year's $r(t)$ adds a new cohort, updating $G_v$ states that affect subsequent years' $S_{\text{prev}}$.

**Constraint**: $r(t) \in [0, 1]$ (ramp cannot go negative or exceed full rate)

---

## Carbon Density Parameters

Each land type has two carbon parameters:
- **Transition density $C_v$** — carbon released or sequestered when land changes category.
- **Equilibrium stock $C^{\text{eq}}_v$** — target for regrowing types; zero for types already at equilibrium.

| Land Type | $C_v$ (tC/km²) | $C^{\text{eq}}_v$ (tC/km²) | Notes |
|-----------|:--:|:---:|-------|
| Primary forest | 15,000 | — | Mature, no regrowth |
| Secondary forest | 8,000 | 12,000 | Relaxes toward primary; key driver |
| Primary non-forest | 1,500 | — | At equilibrium |
| Secondary non-forest | 1,000 | 800 | Recovers toward mature non-forest |
| Cropland (all types) | 500 | — | Managed, no accumulation |
| Pasture | 800 | — | |
| Rangeland | 600 | — | |
| Urban | 200 | — | |
| Plantations | 6,000 | 8,000 | Managed timber regrowth |

**Regrowing types** (those with $C^{\text{eq}} > 0$): `secdf`, `secdn`, `pltns`.

---

## Gridded Extension: Ramped vs Static Variables

Once $r(t)$ is derived, the gridded extension is produced as follows.

### Ramped Variables (State Fractions)

The 13 land-use state fractions (`primf`, `secdf`, `primn`, `secdn`, `c3ann`, `c3nfx`, `c3per`, `c4ann`, `c4per`, `pastr`, `range`, `urban`, `pltns`) evolve during the **active ramp period** (while $r(t) > 0$).

At each grid cell and year:

$$f_v(t+1) = f_v(t) + r(t) \cdot \dot{f}_v^{2100}$$

where $\dot{f}_v^{2100}$ is the rate of change at the last input year. Negative fractions are clamped to zero, and `secdf` (secondary forest) absorbs the residual to enforce per-cell conservation: $\sum_v f_v = 1$.

Once $r(t)$ reaches zero, all state fractions are frozen at their values from the last active ramp year for the remainder of the extension (out to 2500).

### BECCS-Scaled Variables (Biofuel)

Biofuel crop fractions (`crpbiof` and sub-components) are scaled by a time-varying BECCS factor derived from the IAM's bioenergy trajectory, independently of the state ramp. Values are capped at 1.0.

### Static Variables (Management)

Wood harvest, irrigation, and fertiliser variables are held constant at their 2100 values for all extension years. These are written as 2-D fields (no time dimension) in the output.

---

## Verification

The final notebook closes the loop by reading the gridded extension output and independently recomputing AFOLU fluxes using the same carbon model:

1. **Flux match**: Compute transition and stock-change fluxes from the extended grids, apply calibrated coefficients, and compare to the IAM target. RMSE over the ramp period should be small (< 100 Mt CO₂/yr).
2. **Per-cell conservation**: Verify $\sum_v f_v = 1$ at every grid cell and timestep (deviation < $10^{-10}$).
3. **Spatial continuity**: Inspect flux maps at several years to confirm no boundary artefacts at the 2100 join.
4. **Cumulative emissions**: Check that the integral of the flux over the extension period is consistent with the IAM target.

---

## Implementation

The carbon cycle model is implemented in [`src/carbon_model.py`](../src/carbon_model.py):

| Function | Purpose |
|----------|---------|
| `stock_flux_for_tau(tau, area_incr)` | Cohort convolution → stock-change flux timeseries |
| `calibrate_afolu_model(...)` | OLS fit of 4-predictor model with optional τ profile scan |
| `forward_solve_ramp(...)` | Invert calibrated model → AFOLU-consistent $r(t)$ |

Extension logic lives in [`src/extend.py`](../src/extend.py):

| Function | Purpose |
|----------|---------|
| `extend_states(vals, rates, years, ramp)` | Apply ramp to state fractions with per-cell conservation |
| `extend_biofuel(crpbiof, beccs, years)` | Scale biofuel by BECCS factor, cap at 1.0 |
| `extend_constant(field, n_years)` | Replicate 2-D field for hold-constant management |

---

## References

- **Poorter et al. (2016)**: "Biomass resilience of Neotropical secondary forests", *Nature* 530(7589), 211-214.
- **IPCC AR6 WG3**: Carbon cycle timescales for forest regrowth.
- **Hurtt et al. (2020)**: Land-Use Harmonization (LUH2) methodology.
