"""Configuration: paths, scenario names, and constants."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# --- Paths ---
DATA_DIR = Path("/div/no-backup-nac/users/masan/for_Ben")
PROJECT_DIR = Path("/div/no-backup-nac/users/bensan/extend-luh")
OUTPUT_DIR = PROJECT_DIR / "output"

# --- Shared input files ---
CSV_FILE = DATA_DIR / "extensions_full_emissions_timeseries_2023_2500.csv"

# --- Time ---
YR_END_INPUT = 2100       # Last year of gridded input
YR_AFOLU_ZERO = 2149      # Year AFOLU emissions reach zero
YR_END_OUTPUT = 2500      # Final year of extension


# =====================================================================
# Per-scenario configuration
# =====================================================================

@dataclass
class ScenarioConfig:
    """All scenario-specific metadata and file paths."""
    key: str                     # short label, e.g. "VL", "H"
    model: str                   # model name in emissions CSV
    scenario: str                # scenario name in emissions CSV
    folder: str                  # subfolder under DATA_DIR

    # --- Gridded input files (relative to DATA_DIR / folder) ---
    states_file: str             # NetCDF with land-use state timeseries
    biof_file: str               # biofuel fractions
    woodharvest_file: str        # wood harvest (country-level)
    fuelwood_file: str           # fuelwood (country-level)
    mgmt_file: Optional[str] = None       # combined management file (VL-style)
    fertl_file: Optional[str] = None      # separate gridded fertilizer file (M-style)
    flood_file: Optional[str] = None      # separate gridded flood/rice file (M-style)
    protected_file: Optional[str] = None  # separate protected-areas file (M-style)

    # Variable renames applied when writing output (e.g. IMAGE "timber" → LUH2 "pltns")
    var_renames: dict = field(default_factory=dict)

    # --- State variables that sum to 1 per land cell ---
    state_vars: list = field(default_factory=list)

    # --- The residual variable used for per-cell conservation ---
    residual_var: str = "secdf"

    # --- Management variables present in the scenario ---
    mgmt_biofuel_vars: list = field(default_factory=lambda: [
        "crpbiof", "crpbf_c3per", "crpbf_c4per"])
    mgmt_irrig_vars: list = field(default_factory=lambda: [
        f"irrig_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]])
    mgmt_fertl_vars: list = field(default_factory=lambda: [
        f"fertl_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]])
    mgmt_wood_vars: list = field(default_factory=lambda: [
        "rndwd", "fulwd", "pltns_wdprd", "pltns_bfuel"])
    # First-generation biofuel fractions (held constant at 2100)
    mgmt_cpbf1_vars: list = field(default_factory=lambda: [
        f"cpbf1_{c}" for c in ["c3ann", "c4ann", "c3per", "c4per", "c3nfx"]])
    # Extra hold-constant vars that don't fit the wood/irrig/fertl categories
    mgmt_extra_vars: list = field(default_factory=list)

    # If True, management vars are embedded in the states file (no separate mgmt file)
    mgmt_in_states: bool = False

    # --- GLM-input style scenarios (e.g. HL/WITCH) ---
    # States split across several files; merged on load (states_file is then
    # only a nominal stem).
    states_files: list = field(default_factory=list)
    # If set, residual_var is synthesised as max(0, 1 - sum(residual_from)),
    # mirroring GLM's own "other land" residual.
    residual_from: list = field(default_factory=list)
    # Renormalise input cells whose 2100 state sum exceeds 1.  Off for GLM-input
    # scenarios: GLM consumes changes, so rescaling would inject a spurious
    # 2101 land-use change in those cells.
    renormalize_overfull: bool = True
    # How a full cell handles expansion: "scale" squeezes all managed classes
    # proportionally; "limit_expansion" caps only the expanding classes to the
    # room left (non-expanding classes keep their prescribed path).
    overshoot_mode: str = "scale"
    # Gridded vars extended with the ramped 2100 rate but not part of the
    # per-cell conservation (e.g. crop types, which GLM uses only as deltas
    # to split cropland).
    extra_extend_vars: list = field(default_factory=list)
    # Output-only vars built as sums of extended vars, e.g. grazing = pastr+range.
    derived_output_vars: dict = field(default_factory=dict)
    # If set, write one states file per group {input file stem: [vars]} instead
    # of a single states file, mirroring the input file layout.
    state_output_groups: dict = field(default_factory=dict)
    # Country-level management files {filename: "hold" | "zero"}; "zero" is for
    # annual-change variables where a constant level means zero change.
    country_mgmt_files: dict = field(default_factory=dict)
    # Refit without the transition predictor if its coefficient fits negative
    # (unphysical: more conversion would lower emissions).
    constrain_alpha_trans: bool = False

    # Fixed relaxation timescale (years) for calibration.
    # None → profile-likelihood scan; a number → fix τ at that value.
    tau_fixed: Optional[float] = None

    @property
    def data_dir(self) -> Path:
        return DATA_DIR / self.folder

    @property
    def states_path(self) -> Path:
        return self.data_dir / self.states_file

    @property
    def biof_path(self) -> Path:
        return self.data_dir / self.biof_file

    @property
    def woodharvest_path(self) -> Path:
        return self.data_dir / self.woodharvest_file

    @property
    def fuelwood_path(self) -> Path:
        return self.data_dir / self.fuelwood_file

    @property
    def mgmt_path(self) -> Optional[Path]:
        if self.mgmt_file is None:
            return None
        return self.data_dir / self.mgmt_file

    @property
    def fertl_path(self) -> Optional[Path]:
        if self.fertl_file is None:
            return None
        return self.data_dir / self.fertl_file

    @property
    def flood_path(self) -> Optional[Path]:
        if self.flood_file is None:
            return None
        return self.data_dir / self.flood_file

    @property
    def protected_path(self) -> Optional[Path]:
        if self.protected_file is None:
            return None
        return self.data_dir / self.protected_file

    @property
    def output_dir(self) -> Path:
        return OUTPUT_DIR / self.key

    @property
    def mgmt_hold_constant(self) -> list:
        return (self.mgmt_wood_vars + self.mgmt_irrig_vars
                + self.mgmt_fertl_vars + self.mgmt_cpbf1_vars
                + self.mgmt_extra_vars)


# =====================================================================
# Scenario definitions
# =====================================================================

SCENARIOS = {
    "VL": ScenarioConfig(
        key="VL",
        model="REMIND-MAgPIE 3.5-4.11",
        scenario="SSP1 - Very Low Emissions",
        folder="scen7-VL",
        states_file="output_annual_2024_states_step4.nc",
        mgmt_file="output_annual_summed_management_step5.nc",
        biof_file="biof_MAGPIE_LUH3.nc",
        woodharvest_file="woodharvest_MAGPIE_LUH3_with_pltns.nc",
        fuelwood_file="fuelwood_MAGPIE_LUH3_with_pltns_v2.nc",
        state_vars=[
            "primf", "secdf", "primn", "secdn",
            "c3ann", "c3nfx", "c3per", "c4ann", "c4per",
            "pastr", "range", "urban", "pltns",
        ],
        residual_var="secdf",
        tau_fixed=35,
        mgmt_in_states=False,
    ),
    "H": ScenarioConfig(
        key="H",
        model="GCAM 8s",
        scenario="SSP3 - High Emissions",
        folder="scen7-H",
        states_file="annual_timeseries_v3_nan0.nc",
        mgmt_file=None,  # management vars embedded in states file
        biof_file="biof_GCAM_LUH3.nc",
        woodharvest_file="woodharvest_GCAM_LUH3_test2.nc",
        fuelwood_file="fuelwood_GCAM_LUH3.nc",
        state_vars=[
            "forest",
            "c3ann", "c3nfx", "c3per", "c4ann", "c4per",
            "pastr", "range", "urban",
        ],
        residual_var="forest",
        mgmt_biofuel_vars=["crpbiof", "crpbf_c3per", "crpbf_c4per"],
        mgmt_irrig_vars=[
            f"irrig_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]],
        mgmt_fertl_vars=[],  # country-level only, not gridded
        mgmt_wood_vars=[],   # country-level only
        mgmt_in_states=True,
    ),
    "M": ScenarioConfig(
        key="M",
        model="IMAGE 3.4",
        scenario="SSP2 - Medium Emissions",
        folder="scen7-M",
        states_file="output_annual_2024_states_step5.nc",
        biof_file="biof_IMAGE_LUH3.nc",
        woodharvest_file="woodharvest_IMAGE_M_LUH3.nc",
        fuelwood_file="fuelwood_IMAGE_M_LUH3.nc",
        # No combined mgmt file — variables are split across separate files
        mgmt_file=None,
        fertl_file="output_annual_2024_fertilizer_step4.nc",
        flood_file="output_annual_2024_flood_step3.nc",
        protected_file="output_annual_2024_protected_step4.nc",
        # IMAGE states use "timber" for plantation; output renamed to LUH2 "pltns"
        var_renames={"timber": "pltns"},
        state_vars=[
            "primf", "secdf", "primn", "secdn",
            "c3ann", "c3nfx", "c3per", "c4ann", "c4per",
            "pastr", "range", "urban", "timber",
        ],
        residual_var="secdf",
        # No wood or irrigation management from IMAGE — hold-constant via restart
        mgmt_wood_vars=[],
        mgmt_irrig_vars=[],
        mgmt_fertl_vars=[
            f"fertl_{c}" for c in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]],
        mgmt_extra_vars=["flood"],  # rice flooding fraction from flood_file
        mgmt_in_states=False,
        tau_fixed=35,
    ),
    "HL": ScenarioConfig(
        key="HL",
        model="WITCH 6.0",
        scenario="SSP5 - Medium-Low Emissions_a",
        folder="scen7-HL",
        # Delivered as GLM inputs: managed land only, split across files.
        # GLM derives primary/secondary natural land itself, so natural land
        # is represented here as a single residual pool "othr".
        states_file="HL_WITCH_Sept28_states.nc",  # nominal stem only
        states_files=[
            "cropland_HL_WITCH_pltns_Sept28.nc",
            "grazing_HL_WITCH_Sept28.nc",
            "urban_HL_WITCH_Sept28.nc",
            "croptypes_HL_WITCH_Sept28.nc",
        ],
        residual_from=["cropland", "grazing", "urban"],
        biof_file="biof_2nd_gen_HL_WITCH_Sept28.nc",
        woodharvest_file="woodharvest_WITCH_HL_LUH3.nc",
        fuelwood_file="fuelwood_WITCH_HL_LUH3.nc",
        state_vars=["cropland", "pastr", "range", "urban", "othr"],
        residual_var="othr",
        renormalize_overfull=False,
        overshoot_mode="limit_expansion",
        extra_extend_vars=["c3ann", "c4ann", "c3per", "c4per", "c3nfx"],
        derived_output_vars={"grazing": ["pastr", "range"]},
        state_output_groups={
            "cropland_HL_WITCH_pltns_Sept28": ["cropland"],
            "croptypes_HL_WITCH_Sept28": ["c3ann", "c4ann", "c3per", "c4per", "c3nfx"],
            "grazing_HL_WITCH_Sept28": ["grazing", "pastr", "range"],
            "urban_HL_WITCH_Sept28": ["urban"],
        },
        mgmt_in_states=True,  # no gridded mgmt file; vars looked up across inputs
        mgmt_biofuel_vars=["crpbf_c3per"],
        mgmt_irrig_vars=[],
        mgmt_fertl_vars=[],
        mgmt_wood_vars=[],
        mgmt_cpbf1_vars=[],
        country_mgmt_files={
            "fertl_WITCH_LUH3_v4.nc": "zero",   # annual % change → 0 holds level
            "irrig_WITCH_LUH3_v3.nc": "hold",   # irrigated area per country
        },
        constrain_alpha_trans=True,
        tau_fixed=35,
    ),
    "LN": ScenarioConfig(
        key="LN",
        model="AIM 3.0",
        scenario="SSP2 - Low Overshoot_a",
        # Delivered as GLM inputs (like HL); files kept in the repo's
        # incoming/ folder rather than under DATA_DIR.
        folder=str(PROJECT_DIR / "incoming"),
        states_file="LN_AIM_states.nc",  # nominal stem only
        # croptypes file holds all land-use states (incl. forest); flood and
        # protected grids are merged in so their 2100 slices are found as
        # hold-constant management.
        states_files=[
            "croptypes_LN_AIM_v3.nc",
            "flood_grids_LN_AIM.nc",
            "protected_LN_AIM.nc",
        ],
        # Unlike HL, AIM provides forest explicitly, so natural land is split
        # into forest and a non-forest residual pool "othr_nf".
        residual_from=["cropland", "grazing", "urban", "forest"],
        biof_file="biof_grids_LN_AIM.nc",
        woodharvest_file="woodharvest_AIM_LN_LUH3.nc",
        fuelwood_file="fuelwood_AIM_LN_LUH3.nc",
        state_vars=["cropland", "pastr", "range", "urban", "forest", "othr_nf"],
        residual_var="othr_nf",
        renormalize_overfull=False,
        overshoot_mode="limit_expansion",
        extra_extend_vars=["c3ann", "c4ann", "c3per", "c4per", "c3nfx"],
        derived_output_vars={"grazing": ["pastr", "range"]},
        state_output_groups={
            "croptypes_LN_AIM_v3": ["c3ann", "c4ann", "c3per", "c4per", "c3nfx",
                                    "cropland", "grazing", "forest", "urban",
                                    "pastr", "range"],
        },
        mgmt_in_states=True,
        mgmt_biofuel_vars=["crpbf_c4per"],
        mgmt_irrig_vars=[],
        mgmt_fertl_vars=[],
        mgmt_wood_vars=[],
        mgmt_cpbf1_vars=[],
        mgmt_extra_vars=["flood", "primf_prot", "secdf_prot", "primn_secdn_prot"],
        country_mgmt_files={
            "AIM_LN_fert_country_v3.nc": "zero",       # annual % change
            "irrig_AIM_LUH3_croptypes_v3.nc": "hold",  # irrigated area per country
            "flood_AIM_LUH3_v3.nc": "zero",            # labelled % change (TBC with Louise)
        },
        constrain_alpha_trans=True,
        tau_fixed=35,
    ),
}


def get_scenario(key: str) -> ScenarioConfig:
    """Look up a scenario by its short key (e.g. 'VL', 'H')."""
    if key not in SCENARIOS:
        raise KeyError(
            f"Unknown scenario '{key}'. Available: {list(SCENARIOS.keys())}")
    return SCENARIOS[key]


# =====================================================================
# Default scenario (backward compatibility)
# =====================================================================
DEFAULT_SCENARIO_KEY = "VL"
_default = SCENARIOS[DEFAULT_SCENARIO_KEY]

MODEL = _default.model
SCENARIO = _default.scenario
STATES_FILE = _default.states_path
MGMT_FILE = _default.mgmt_path
BIOF_FILE = _default.biof_path
WOODHARVEST_FILE = _default.woodharvest_path
FUELWOOD_FILE = _default.fuelwood_path
STATE_VARS = _default.state_vars

# --- Management variables (backward compat) ---
MGMT_BIOFUEL_VARS = _default.mgmt_biofuel_vars
MGMT_WOOD_VARS = _default.mgmt_wood_vars
MGMT_IRRIG_VARS = _default.mgmt_irrig_vars
MGMT_FERTL_VARS = _default.mgmt_fertl_vars
MGMT_HOLD_CONSTANT = _default.mgmt_hold_constant
