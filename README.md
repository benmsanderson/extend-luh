# extend-luh

Extend LUH2 land-use harmonization grids from 2100 to 2500.

## Method

1. **AFOLU calibration** (notebook 02): A cohort-based carbon cycle model is calibrated against the IAM's AFOLU CO₂ trajectory pre-2100, then inverted to derive the rate multiplier $r(t)$ that reproduces the target post-2100 AFOLU pathway.
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

## Repository layout

```
src/
  config.py    # paths, scenario constants, variable lists
  data.py      # CSV / NetCDF loaders and preprocessing
  extend.py    # core extension logic (states, biofuel, constant)
notebooks/
  01_explore.ipynb               # exploratory analysis of inputs
  02_afolu_calibration.ipynb     # calibrate AFOLU model, derive ramp r(t)
  03_gridded_extension.ipynb     # produce extended NetCDF files using r(t)
  04_verify_afolu.ipynb          # verify AFOLU fluxes from gridded output
output/                          # generated files (git-ignored)
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
