"""Diagnostic plots saved as PNGs to sc.output_dir after a completed extension run."""

import glob as _glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xarray as xr

from . import config as cfg
from . import carbon_model as cm
from .data import load_states, load_csv, filter_scenario, get_variable, load_derived_ramp
from .pipeline import _cell_area


# ── Colour palette shared across plots ──────────────────────────────────────
_PALETTE = {
    'Primary forest':       '#1b5e20',
    'Secondary forest':     '#4caf50',
    'Forest':               '#2e7d32',
    'Primary non-forest':   '#5d4037',
    'Secondary non-forest': '#a1887f',
    'Cropland':             '#ff9800',
    'Pasture':              '#e91e63',
    'Rangeland':            '#f48fb1',
    'Urban':                '#757575',
    'Plantations':          '#8bc34a',
}

_GROUP_DEFS = [
    ('Primary forest',       ['primf']),
    ('Secondary forest',     ['secdf']),
    ('Forest',               ['forest']),
    ('Primary non-forest',   ['primn']),
    ('Secondary non-forest', ['secdn']),
    ('Cropland',             ['c3ann', 'c3nfx', 'c3per', 'c4ann', 'c4per']),
    ('Pasture',              ['pastr']),
    ('Rangeland',            ['range']),
    ('Urban',                ['urban']),
    ('Plantations',          ['pltns', 'timber']),
]


def _find_ext_states(sc):
    """Return path to the extension states NetCDF, or None."""
    candidates = sorted(_glob.glob(
        str(sc.output_dir / "*_extension_transient_2101-*.nc")))
    # Prefer a file that has grid dimensions (states, not biofuel-only)
    for p in candidates:
        try:
            with xr.open_dataset(p) as ds:
                if 'lat' in ds.dims and len(ds.data_vars) > 1:
                    return p
        except Exception:
            continue
    return candidates[0] if candidates else None


def _out_to_src_name(v, sc):
    """Reverse-map an output variable name back to its source name."""
    inv = {out: src for src, out in sc.var_renames.items()}
    return inv.get(v, v)


def plot_diagnostics(sc, *, verbose=True):
    """Generate diagnostic PNGs for a completed extension run.

    Saves three files to sc.output_dir:
        diag_01_land_timeseries.png  — global land-use stacked-area
        diag_02_biomass.png          — scndMeanAge biomass trajectory
        diag_03_verification.png     — IAM AFOLU target vs reconstruction
    """
    log = print if verbose else (lambda *a, **k: None)

    ext_path = _find_ext_states(sc)
    cal_path = sc.output_dir / "afolu_calibration.npz"

    if ext_path is None:
        log(f"  {sc.key}: no extension states file — skipping plots")
        return
    if not cal_path.exists():
        log(f"  {sc.key}: no calibration file — skipping plots")
        return

    calib = np.load(cal_path)
    tau = float(calib["tau"])

    ds_ext = xr.open_dataset(ext_path)
    ds_src = load_states(scen=sc)

    ext_years = ds_ext.time.values.astype(int)
    src_years = ds_src.time.values.astype(int)

    lat = ds_ext.lat.values
    lon = ds_ext.lon.values
    cell_area = _cell_area(lat, lon)  # (nlat, nlon)

    # Output var names (used in ext file); source var names (used in src file)
    # For M: out_vars has 'pltns', src_vars has 'timber'; for VL/H they match.
    out_vars = [sc.var_renames.get(v, v) for v in sc.state_vars]

    # ── Plot 1: Global land-use timeseries ───────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(13, 5))

    # Build groups present in this scenario
    all_out = set(out_vars)
    groups = {
        label: [v for v in vlist if v in all_out]
        for label, vlist in _GROUP_DEFS
    }
    groups = {k: v for k, v in groups.items() if v}

    # Pre-2100 global sums from source file (last 50 years to keep it readable)
    src_slice = src_years[-50:] if len(src_years) > 50 else src_years
    src_data = {}
    for label, out_vlist in groups.items():
        ts = np.zeros(len(src_slice))
        for out_v in out_vlist:
            src_v = _out_to_src_name(out_v, sc)
            if src_v in ds_src:
                ti_start = len(src_years) - len(src_slice)
                ts += np.array([
                    np.nansum(ds_src[src_v].isel(time=ti_start + ti).values * cell_area)
                    for ti in range(len(src_slice))
                ])
        src_data[label] = ts

    # Extension global sums
    ext_data = {}
    for label, out_vlist in groups.items():
        ts = np.zeros(len(ext_years))
        for out_v in out_vlist:
            if out_v in ds_ext:
                ts += np.array([
                    np.nansum(ds_ext[out_v].isel(time=ti).values * cell_area)
                    for ti in range(len(ext_years))
                ])
        ext_data[label] = ts

    # Stacked area
    labels = list(groups.keys())
    colors = [_PALETTE[l] for l in labels]
    all_x = np.concatenate([src_slice, ext_years])
    all_y = np.vstack([
        np.concatenate([src_data[l], ext_data[l]]) for l in labels
    ])

    ax1.stackplot(all_x, all_y, labels=labels, colors=colors, alpha=0.85)
    ax1.axvline(2100, color='k', ls='--', lw=0.8, label='Extension start')
    ax1.axvline(int(ext_years[-1]), color='gray', ls=':', lw=0.8,
                label=f'Ramp ends ({int(ext_years[-1])})')
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Area (km²)')
    ax1.set_title(f'{sc.key} — {sc.model}: Global land-use allocation')
    ax1.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=8)
    fig1.tight_layout()
    p1 = sc.output_dir / "diag_01_land_timeseries.png"
    fig1.savefig(p1, dpi=150, bbox_inches='tight')
    plt.close(fig1)
    log(f"  {sc.key}: saved {p1.name}")

    # ── Plot 2: Biomass trajectory (scndMeanAge model) ───────────────────────
    regrow_out = [v for v in out_vars if cm.C_EQ.get(v, 0) > 0]
    regrow_src = [_out_to_src_name(v, sc) for v in regrow_out]

    _ceq = {v: cm.C_EQ.get(v, 0) for v in sc.state_vars}

    def _gbm(frac, density):
        return np.nansum(frac * cell_area * density) / 1e9  # GtC

    # Pre-2100 biomass (static density — same as calibration period)
    biomass_pre = []
    for ti in range(len(src_years)):
        b = 0.0
        for sv, ov in zip(sc.state_vars, out_vars):
            frac = ds_src[sv].isel(time=ti).values if sv in ds_src else None
            if frac is None:
                continue
            b += _gbm(frac, cm.CARBON_DENSITY.get(sv, 0))
        biomass_pre.append(b)
    biomass_pre = np.array(biomass_pre)

    # Initialise scndMeanAge at 2100
    vals_2100 = {}
    for sv, ov in zip(sc.state_vars, out_vars):
        if sv in ds_src:
            vals_2100[ov] = ds_src[sv].isel(time=-1).values
        elif ov in ds_src:
            vals_2100[ov] = ds_src[ov].isel(time=-1).values

    sma = {}
    area_prev = {}
    for ov in regrow_out:
        cd = cm.CARBON_DENSITY.get(ov, 0)
        ceq = _ceq.get(ov, 0) or _ceq.get(_out_to_src_name(ov, sc), 0)
        if ceq > 0:
            sma[ov] = -tau * np.log(1.0 - min(cd / ceq, 0.9999))
        else:
            sma[ov] = 0.0
        frac0 = vals_2100.get(ov, np.zeros((len(lat), len(lon))))
        area_prev[ov] = np.nansum(frac0 * cell_area)

    biomass_ext = []
    for ti in range(len(ext_years)):
        b = 0.0
        for ov in out_vars:
            sv = _out_to_src_name(ov, sc)
            frac = ds_ext[ov].isel(time=ti).values if ov in ds_ext else None
            if frac is None:
                continue
            if ov in regrow_out:
                a_curr = np.nansum(frac * cell_area)
                a_prev_v = area_prev[ov]
                a_gained = max(0.0, a_curr - a_prev_v)
                a_retained = a_prev_v - max(0.0, a_prev_v - a_curr)
                if a_curr > 0:
                    sma[ov] = (a_retained * (sma[ov] + 1.0) + a_gained * 0.0) / a_curr
                else:
                    sma[ov] += 1.0
                area_prev[ov] = a_curr
                ceq = _ceq.get(ov, _ceq.get(sv, 0))
                density = ceq * (1.0 - np.exp(-sma[ov] / tau))
            else:
                density = cm.CARBON_DENSITY.get(ov, cm.CARBON_DENSITY.get(sv, 0))
            b += _gbm(frac, density)
        biomass_ext.append(b)
    biomass_ext = np.array(biomass_ext)

    all_years_bm = np.concatenate([src_years, ext_years])
    all_bm = np.concatenate([biomass_pre, biomass_ext])

    fig2, ax2 = plt.subplots(figsize=(12, 4))
    ax2.plot(all_years_bm, all_bm, lw=2, color='steelblue')
    ax2.axvline(2100, color='k', ls='--', lw=0.8, label='Extension start')
    ax2.axvline(int(ext_years[-1]), color='gray', ls=':', lw=0.8,
                label=f'Ramp ends ({int(ext_years[-1])})')
    ax2.set_xlabel('Year')
    ax2.set_ylabel('Vegetation biomass (GtC)')
    ax2.set_title(f'{sc.key} — {sc.model}: scndMeanAge biomass trajectory  [τ = {tau:.0f} yr]')
    ax2.legend()
    ax2.ticklabel_format(useOffset=False)
    fig2.tight_layout()
    p2 = sc.output_dir / "diag_02_biomass.png"
    fig2.savefig(p2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    log(f"  {sc.key}: saved {p2.name}")

    # ── Plot 3: AFOLU verification ────────────────────────────────────────────
    beta = float(calib["beta"])
    gamma = float(calib["gamma"])
    alpha_trans = float(calib["alpha_trans"])
    alpha_stock = float(calib["alpha_stock"])

    TC = cm.TC_TO_MTCO2
    REGROW_VARS_SRC = [v for v in sc.state_vars if cm.C_EQ.get(v, 0) > 0]

    # IAM target (full timeseries from CSV)
    df = load_csv()
    sc_data = filter_scenario(df, scen=sc)
    afolu_target = get_variable(sc_data, "Emissions|CO2|AFOLU")
    tgt_years = np.array([int(float(c)) for c in afolu_target.index])
    tgt_vals = afolu_target.values

    # Reconstruct transition flux from source states
    src_input_years = ds_src.time.values.astype(int)
    n_src = len(src_input_years)
    trans_src = np.zeros(n_src - 1)
    area_incr_src = {v: np.zeros(n_src - 1) for v in REGROW_VARS_SRC}
    for ti in range(1, n_src):
        af = 0.0
        for sv in sc.state_vars:
            if sv not in ds_src:
                continue
            delta = (ds_src[sv].isel(time=ti).values
                     - ds_src[sv].isel(time=ti - 1).values)
            af += np.nansum(-delta * cm.CARBON_DENSITY.get(sv, 0) * cell_area) * TC
        trans_src[ti - 1] = af
        for v in REGROW_VARS_SRC:
            if v not in ds_src:
                continue
            delta = (ds_src[v].isel(time=ti).values
                     - ds_src[v].isel(time=ti - 1).values)
            area_incr_src[v][ti - 1] = np.nansum(delta * cell_area)

    flux_years_src = src_input_years[1:]
    time_src = (flux_years_src - flux_years_src.mean()) / 10.0

    # Stock-change flux from source period
    stock_src = np.zeros(n_src - 1)
    cohort_stock = {v: 0.0 for v in REGROW_VARS_SRC}
    for ti in range(n_src - 1):
        for v in REGROW_VARS_SRC:
            ceq_v = _ceq.get(v, cm.C_EQ.get(v, 0))
            decay = np.exp(-1.0 / tau) if tau > 0 else 0
            cohort_stock[v] = cohort_stock[v] * decay + area_incr_src[v][ti] * ceq_v * (1 - decay)
            stock_src[ti] -= cohort_stock[v] * TC

    recon_src = beta + gamma * time_src + alpha_trans * trans_src + alpha_stock * stock_src

    # Extension period reconstruction
    EXT_ALL = np.arange(cfg.YR_END_INPUT + 1, cfg.YR_END_OUTPUT + 1)
    derived_ramp_full = load_derived_ramp(EXT_ALL, scen=sc)
    trans_ext = np.zeros(len(ext_years))
    area_incr_ext = {v: np.zeros(len(ext_years)) for v in REGROW_VARS_SRC}
    vals_2100_src = {sv: ds_src[sv].isel(time=-1).values
                     for sv in sc.state_vars if sv in ds_src}

    for ti in range(len(ext_years)):
        af = 0.0
        for sv, ov in zip(sc.state_vars, out_vars):
            if ov not in ds_ext:
                continue
            prev = (vals_2100_src.get(sv, np.zeros((len(lat), len(lon))))
                    if ti == 0 else ds_ext[ov].isel(time=ti - 1).values)
            curr = ds_ext[ov].isel(time=ti).values
            delta = curr - prev
            af += np.nansum(-delta * cm.CARBON_DENSITY.get(sv, cm.CARBON_DENSITY.get(ov, 0)) * cell_area) * TC
        trans_ext[ti] = af
        for sv in REGROW_VARS_SRC:
            ov = sc.var_renames.get(sv, sv)
            if ov not in ds_ext:
                continue
            prev = (vals_2100_src.get(sv, np.zeros((len(lat), len(lon))))
                    if ti == 0 else ds_ext[ov].isel(time=ti - 1).values)
            curr = ds_ext[ov].isel(time=ti).values
            area_incr_ext[sv][ti] = np.nansum((curr - prev) * cell_area)

    time_ext_ref = flux_years_src.mean()
    time_ext = (ext_years - time_ext_ref) / 10.0

    stock_ext = np.zeros(len(ext_years))
    for v in REGROW_VARS_SRC:
        cohort_stock[v] = 0.0  # reset for extension
    for ti in range(len(ext_years)):
        for v in REGROW_VARS_SRC:
            ceq_v = _ceq.get(v, cm.C_EQ.get(v, 0))
            decay = np.exp(-1.0 / tau) if tau > 0 else 0
            cohort_stock[v] = cohort_stock[v] * decay + area_incr_ext[v][ti] * ceq_v * (1 - decay)
            stock_ext[ti] -= cohort_stock[v] * TC

    r_ext = derived_ramp_full[:len(ext_years)]
    recon_ext = beta + gamma * time_ext + alpha_trans * trans_ext * r_ext + alpha_stock * stock_ext

    fig3, ax3 = plt.subplots(figsize=(12, 4))
    # IAM target
    mask = (tgt_years >= flux_years_src[0]) & (tgt_years <= int(ext_years[-1]))
    ax3.plot(tgt_years[mask], tgt_vals[mask], 'k-', lw=1.5, label='IAM target')
    # Reconstructed (source period)
    ax3.plot(flux_years_src, recon_src, color='steelblue', lw=1.5, ls='--',
             label='Reconstructed (source)')
    # Reconstructed (extension period)
    ax3.plot(ext_years, recon_ext, color='steelblue', lw=2, label='Reconstructed (extension)')
    ax3.axvline(2100, color='gray', ls='--', lw=0.8)
    ax3.axhline(0, color='gray', lw=0.4)
    ax3.set_xlabel('Year')
    ax3.set_ylabel('AFOLU CO₂ (Mt CO₂/yr)')
    ax3.set_title(f'{sc.key} — {sc.model}: AFOLU flux — IAM target vs reconstruction')
    ax3.legend()
    ax3.ticklabel_format(useOffset=False)
    fig3.tight_layout()
    p3 = sc.output_dir / "diag_03_verification.png"
    fig3.savefig(p3, dpi=150, bbox_inches='tight')
    plt.close(fig3)
    log(f"  {sc.key}: saved {p3.name}")

    ds_ext.close()
    ds_src.close()
