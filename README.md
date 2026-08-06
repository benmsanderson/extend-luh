# extend-luh

Extend LUH2 land-use harmonization grids from 2100 to 2500.

## Method

1. **AFOLU calibration** (notebook 02): A cohort-based carbon cycle model is calibrated against the IAM's **globally-aggregated** AFOLU CO₂ trajectory pre-2100, then inverted to derive the rate multiplier $r(t)$ that reproduces the target post-2100 AFOLU pathway. The model is zero-dimensional: all fluxes and area increments are summed globally before fitting, yielding a single scalar $r(t)$ applied uniformly to every grid cell.
2. **Gridded extension** (notebook 03): Per-cell rates of change are multiplied by the AFOLU-consistent ramp $r(t)$, with negative fractions clamped and deficit absorbed by secondary forest.
3. Biofuel crop fractions are scaled by an IAM-derived BECCS ratio.
4. All other management fields are held constant at their 2100 values.
5. **Verification** (notebook 04): Extension output is read back and AFOLU fluxes are recomputed to confirm consistency.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running the pipeline

### Automated (recommended)

`run_notebooks.py` uses [papermill](https://papermill.readthedocs.io) to execute notebooks 02, 03, and 04 non-interactively for one or more scenarios. Executed notebooks (with all plots and diagnostic outputs) are saved to `output/{KEY}/notebooks/` as an audit trail.

```bash
# Run all available scenarios end-to-end
python run_notebooks.py

# Run specific scenarios
python run_notebooks.py VL H

# Skip calibration if already done (e.g. re-running extension only)
python run_notebooks.py VL H --skip-calibration

# Discard executed notebooks (faster, no audit trail)
python run_notebooks.py VL --no-save
```

Flags: `--skip-calibration`, `--skip-extension`, `--skip-verify`, `--no-save`

### Interactive

Open notebooks 02–04 in Jupyter, set `SCENARIO_KEY` in the first cell, and run all cells. Useful for development and detailed diagnostics.

## Repository layout

```
src/
  config.py      # paths, scenario constants, variable lists
  data.py        # CSV / NetCDF loaders and preprocessing
  extend.py      # core extension logic (states, biofuel, harvest)
  carbon_model.py # cohort-based AFOLU carbon model
  pipeline.py    # function-level equivalents of notebooks 02–04
notebooks/
  01_explore.ipynb                    # exploratory analysis of inputs
  02_afolu_calibration_multi.ipynb    # calibrate AFOLU model, derive ramp r(t)
  03_gridded_extension.ipynb          # produce extended NetCDF files using r(t)
  04_verify_afolu.ipynb               # verify AFOLU fluxes from gridded output
run_notebooks.py                      # CLI for automated multi-scenario execution
output/                               # generated files (git-ignored)
```

## Output files

Output filenames mirror the input file stems with an `_extension_` suffix.

| Category | Contents | Time axis |
|:---------|:---------|:---------:|
| States (transient) | 13 state fractions | ramp period only |
| Biofuel (transient) | crpbiof | full extension |
| Management (transient) | BECCS-scaled biofuel sub-components | full extension |
| Management (static) | hold-constant mgmt vars (wood, irrigation, fertiliser) | single timestep |

All files: zlib level 4, float32, CF-1.6 conventions, 0.25° grid (720 × 1440).
