# extend-luh

Extend LUH2 land-use harmonization grids from 2100 to 2500.

## Method

- Per-cell rates of change are ramped linearly to zero over a configurable period (default: 2100–2149).
- Negative fractions are clamped; the deficit is absorbed by secondary forest.
- Biofuel crop fractions are scaled by an IAM-derived BECCS ratio.
- All other management fields are held constant at their 2100 values.

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
  02_gridded_extension.ipynb     # produce extended NetCDF files
  03_afolu_flux_check.ipynb      # AFOLU flux estimation & validation
output/                          # generated NetCDF files (git-ignored)
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
