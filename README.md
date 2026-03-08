# extend-luh

Extend LUH2 land-use harmonization grids from 2100 to 2500 for
**REMIND-MAgPIE 3.5-4.11 / SSP1 – Very Low Emissions**.

## Method

- Per-cell rates of change are ramped linearly to zero between 2100 and 2149.
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

| File | Variables | Timesteps |
|:-----|:----------|:---------:|
| `output_annual_2024_states_step4_extension_transient_2101-2150.nc` | 13 state fractions | 50 |
| `biof_MAGPIE_LUH3_extension_transient_2101-2500.nc` | crpbiof | 400 |
| `output_annual_summed_management_step5_extension_transient_2101-2500.nc` | crpbf_c3per, crpbf_c4per | 400 |
| `output_annual_summed_management_step5_extension_static.nc` | 14 hold-constant mgmt vars | 1 |

All files: zlib level 4, float32, CF-1.6 conventions, 0.25° grid (720 × 1440).
