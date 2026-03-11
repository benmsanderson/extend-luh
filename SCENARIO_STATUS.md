# LUH Extension Processing Tracker

Track the status of each scenario through the processing pipeline.

## Processing Stages

1. **AFOLU Calibration** (notebook 02) — Requires: emissions CSV
2. **Gridded Extension** (notebook 03) — Requires: states file, biof file, AFOLU calibration
3. **Verification** (notebook 04) — Requires: extended output from stage 2

## Scenario Status

| Key | Model | Scenario | CSV | States | Biof | Stage 1 | Stage 2 | Stage 3 | Notes |
|-----|-------|----------|:---:|:------:|:----:|:-------:|:-------:|:-------:|-------|
| **VL** | REMIND-MAgPIE 3.5-4.11 | SSP1 - Very Low Emissions | ✓ | ✓ step4/5 | ✓ | ✓ | ✓ | ✓ | Reference (complete) |
| **H** | GCAM 8s | SSP3 - High Emissions | ✓ | ✓ timeseries | ✓ | ⏳ | ⏳ | — | Ready for processing |
| **HL** | WITCH 6.0 | SSP5 - Medium-Low Emissions_a | ✓ | ✗ | ✗ | ⏳ | ⊗ | ⊗ | CSV only |
| **L** | MESSAGEix-GLOBIOM-GAINS 2.1-M-R12 | SSP2 - Low Emissions | ✓ | ✗ | ✗ | ⏳ | ⊗ | ⊗ | CSV only |
| **LN** | AIM 3.0 | SSP2 - Low Overshoot_a | ✓ | ✗ | ✗ | ⏳ | ⊗ | ⊗ | CSV only |
| **M** | IMAGE 3.4 | SSP2 - Medium Emissions | ✓ | ✗ | ✗ | ⏳ | ⊗ | ⊗ | CSV only |
| **ML** | COFFEE 1.6 | SSP2 - Medium-Low Emissions | ✓ | ✗ | ✗ | ⏳ | ⊗ | ⊗ | CSV only |

**Legend:**
- ✓ = Complete
- ⏳ = Ready to run
- ⊗ = Cannot run (missing dependencies)
- — = Not yet ready

## Data Availability by Scenario

### VL (REMIND-MAgPIE, SSP1 Very Low)
- ✓ States: `output_annual_2024_states_step4.nc` (13 vars: primf, secdf, primn, secdn, crops, pastr, range, urban, pltns)
- ✓ Management: `output_annual_summed_management_step5.nc` (separate file)
- ✓ Biofuel: `biof_MAGPIE_LUH3.nc`
- ✓ Wood/fuelwood: country-level

### H (GCAM, SSP3 High)
- ✓ States: `annual_timeseries_v3_nan0.nc` (9 vars: forest, crops, pastr, range, urban)
- ✓ Management: embedded in states file (irrig, crpbf)
- ✓ Biofuel: `biof_GCAM_LUH3.nc`
- ✓ Wood/fuelwood: country-level
- ⚠ **Different variable structure**: `forest` instead of primf/secdf/pltns split

### HL, L, LN, M, ML
- ✓ Emissions CSV (global)
- ⊗ No gridded states
- ⊗ No biofuel grids
- Expected: Gridded data may be added later

## Processing Workflow

### For scenarios with full gridded data (VL, H):
1. Run `02_afolu_calibration_multi.ipynb` with appropriate SCENARIO_KEY
2. Run `03_gridded_extension.ipynb` (needs parameterization like notebook 02)
3. Run `04_verify_afolu.ipynb` (needs parameterization)

### For CSV-only scenarios (HL, L, LN, M, ML):
1. Run `02_afolu_calibration_multi.ipynb` — will use simplified calibration
2. Output: AFOLU-consistent ramp for comparison/analysis
3. When gridded data becomes available, can proceed to stage 2

## Next Steps

1. ✅ Parameterize notebook 02 for multi-scenario support
2. ⏳ Parameterize notebook 03 for multi-scenario support
3. ⏳ Parameterize notebook 04 for multi-scenario support
4. ⏳ Test full pipeline on H scenario
5. ⏳ Run CSV-only calibrations for HL, L, LN, M, ML
6. ⏳ Document scenario-specific carbon density assumptions

## Carbon Model Considerations

**VL-specific parameters (well-calibrated)**:
- 13-variable state with primary/secondary forest split
- secdf as residual variable
- Regrowing: secdf (12k tC/km²), secdn (800), pltns (8k)

**H-specific parameters (initial estimates)**:
- 9-variable state with aggregate forest
- forest as residual variable  
- Regrowing: forest (10k tC/km² — may need tuning)
- ⚠ Forest carbon density (10k) is a rough estimate; adjust based on calibration quality

**Action**: After H calibration, review R² and residuals. If poor fit, consider adjusting:
- `CARBON_DENSITY['forest']` (currently 10k tC/km²)
- `C_EQ['forest']` (currently 10k tC/km²)
