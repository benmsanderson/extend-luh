"""Configuration: paths, scenario names, and constants."""

from pathlib import Path

# --- Paths ---
DATA_DIR = Path("/div/no-backup-nac/users/masan/for_Ben")
PROJECT_DIR = Path("/div/no-backup-nac/users/bensan/extend-luh")
OUTPUT_DIR = PROJECT_DIR / "output"

# --- Input files ---
CSV_FILE = DATA_DIR / "extensions_full_emissions_timeseries_2023_2500.csv"
STATES_FILE = DATA_DIR / "output_annual_2024_states_step4.nc"
MGMT_FILE = DATA_DIR / "output_annual_summed_management_step5.nc"
BIOF_FILE = DATA_DIR / "biof_MAGPIE_LUH3.nc"
WOODHARVEST_FILE = DATA_DIR / "woodharvest_MAGPIE_LUH3_with_pltns.nc"
FUELWOOD_FILE = DATA_DIR / "fuelwood_MAGPIE_LUH3_with_pltns_v2.nc"

# --- Scenario ---
MODEL = "REMIND-MAgPIE 3.5-4.11"
SCENARIO = "SSP1 - Very Low Emissions"

# --- Time ---
YR_END_INPUT = 2100       # Last year of gridded input
YR_AFOLU_ZERO = 2149      # Year AFOLU emissions reach zero
YR_END_OUTPUT = 2500      # Final year of extension

# --- Land-use state variables (sum to 1 per cell) ---
STATE_VARS = [
    "primf", "secdf", "primn", "secdn",
    "c3ann", "c3nfx", "c3per", "c4ann", "c4per",
    "pastr", "range", "urban", "pltns",
]

# --- Management variables ---
MGMT_BIOFUEL_VARS = ["crpbiof", "crpbf_c3per", "crpbf_c4per"]
MGMT_WOOD_VARS = ["rndwd", "fulwd", "pltns_wdprd", "pltns_bfuel"]
MGMT_IRRIG_VARS = [f"irrig_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]]
MGMT_FERTL_VARS = [f"fertl_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]]
MGMT_HOLD_CONSTANT = MGMT_WOOD_VARS + MGMT_IRRIG_VARS + MGMT_FERTL_VARS
