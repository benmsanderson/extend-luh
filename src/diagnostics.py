"""Diagnostic plots saved as PNGs to sc.output_dir after a completed extension run."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xarray as xr

from . import config as cfg
from . import carbon_model as cm
from .data import load_states, load_csv, filter_scenario, get_variable, load_derived_ramp
from .pipeline import _cell_area, afolu_reconstruction, load_extension_states


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
    'Other natural':        '#c5b358',
}

_GROUP_DEFS = [
    ('Primary forest',       ['primf']),
    ('Secondary forest',     ['secdf']),
    ('Forest',               ['forest']),
    ('Primary non-forest',   ['primn']),
    ('Secondary non-forest', ['secdn']),
    ('Cropland',             ['cropland', 'c3ann', 'c3nfx', 'c3per', 'c4ann', 'c4per']),
    ('Pasture',              ['pastr']),
    ('Rangeland',            ['range']),
    ('Urban',                ['urban']),
    ('Plantations',          ['pltns', 'timber']),
    ('Other natural',        ['othr', 'othr_nf']),
]



def _out_to_src_name(v, sc):
    """Reverse-map an output variable name back to its source name."""
    inv = {out: src for src, out in sc.var_renames.items()}
    return inv.get(v, v)


def _icwtr(shape, log):
    """LUH2 ice/water fraction on the 0.25° grid (zeros if unavailable)."""
    try:
        with xr.open_dataset(cfg.LUH_STATIC_FILE, decode_times=False) as ds:
            ic = np.nan_to_num(ds["icwtr"].values.astype(np.float64))
        if ic.shape == shape:
            return ic
        log(f"  icwtr grid {ic.shape} != {shape}; residual pool not masked")
    except Exception as exc:
        log(f"  could not read {cfg.LUH_STATIC_FILE} ({exc}); residual pool not masked")
    return np.zeros(shape)


def plot_diagnostics(sc, *, verbose=True):
    """Generate diagnostic PNGs for a completed extension run.

    Saves three files to sc.output_dir:
        diag_01_land_timeseries.png  — global land-use stacked-area
        diag_02_biomass.png          — scndMeanAge biomass trajectory
        diag_03_verification.png     — IAM AFOLU target vs reconstruction
    """
    log = print if verbose else (lambda *a, **k: None)

    cal_path = sc.output_dir / "afolu_calibration.npz"
    if not cal_path.exists():
        log(f"  {sc.key}: no calibration file — skipping plots")
        return

    ds_src = load_states(scen=sc)
    ds_ext = load_extension_states(sc, ds_src)
    if ds_ext is None:
        log(f"  {sc.key}: no extension states file — skipping plots")
        ds_src.close()
        return

    tau = float(np.load(cal_path)["tau"])

    ext_years = ds_ext.time.values.astype(int)
    src_years = ds_src.time.values.astype(int)

    lat = ds_ext.lat.values
    lon = ds_ext.lon.values
    cell_area = _cell_area(lat, lon)  # (nlat, nlon)

    # A synthesised residual pool (GLM-input scenarios: 1 − managed land)
    # also covers ocean and ice; drop the ice/water fraction for area and
    # biomass totals.  Fluxes are unaffected (they use changes only).
    icwtr = _icwtr(cell_area.shape, log) if sc.residual_from else None

    def _frac(ds, v, ti):
        a = ds[v].isel(time=ti).values
        if icwtr is not None and v == sc.residual_var:
            a = np.maximum(np.nan_to_num(a) - icwtr, 0.0)
        return a

    # Output var names (used in ext file); source var names (used in src file)
    # For M: out_vars has 'pltns', src_vars has 'timber'; for VL/H they match.
    out_vars = [sc.var_renames.get(v, v) for v in sc.state_vars]

    # ── Plot 1: Global land-use timeseries ───────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(13, 5))

    all_out = set(out_vars)
    groups = {
        label: [v for v in vlist if v in all_out]
        for label, vlist in _GROUP_DEFS
    }
    groups = {k: v for k, v in groups.items() if v}

    # Pre-2100 global sums from source file (last 50 years to keep it readable)
    src_slice = src_years[-50:] if len(src_years) > 50 else src_years
    ti_start = len(src_years) - len(src_slice)
    src_data = {}
    for label, out_vlist in groups.items():
        ts = np.zeros(len(src_slice))
        for out_v in out_vlist:
            src_v = _out_to_src_name(out_v, sc)
            if src_v in ds_src:
                ts += np.array([
                    np.nansum(_frac(ds_src, src_v, ti_start + ti) * cell_area)
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
                    np.nansum(_frac(ds_ext, out_v, ti) * cell_area)
                    for ti in range(len(ext_years))
                ])
        ext_data[label] = ts

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
    regrow_out = [v for v in out_vars if cm.C_EQ.get(_out_to_src_name(v, sc), 0) > 0]

    def _gbm(frac, density):
        return np.nansum(frac * cell_area * density) / 1e9  # GtC

    # Pre-2100 biomass (static density — same as calibration period)
    biomass_pre = []
    for ti in range(len(src_years)):
        b = 0.0
        for sv in sc.state_vars:
            if sv in ds_src:
                b += _gbm(_frac(ds_src, sv, ti), cm.CARBON_DENSITY.get(sv, 0))
        biomass_pre.append(b)
    biomass_pre = np.array(biomass_pre)

    # Initialise scndMeanAge at 2100 from the static density.  Regrowing
    # density is cd · f(age)/f(age₀), f = 1 − exp(−age/τ), so it equals the
    # static density at 2100 even where cd ≥ C_eq (e.g. secdn, othr_nf).
    sma, f0, area_prev = {}, {}, {}
    for ov in regrow_out:
        sv = _out_to_src_name(ov, sc)
        cd = cm.CARBON_DENSITY.get(sv, 0)
        ceq = cm.C_EQ[sv]
        f0[ov] = min(cd / ceq, 0.9999)
        sma[ov] = -tau * np.log(1.0 - f0[ov])
        area_prev[ov] = np.nansum(_frac(ds_src, sv, -1) * cell_area)

    biomass_ext = []
    for ti in range(len(ext_years)):
        b = 0.0
        for ov in out_vars:
            if ov not in ds_ext:
                continue
            sv = _out_to_src_name(ov, sc)
            frac = _frac(ds_ext, ov, ti)
            cd = cm.CARBON_DENSITY.get(sv, cm.CARBON_DENSITY.get(ov, 0))
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
                density = cd * (1.0 - np.exp(-sma[ov] / tau)) / f0[ov]
            else:
                density = cd
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

    ds_ext.close()
    ds_src.close()

    # ── Plot 3: AFOLU verification ────────────────────────────────────────────
    # Same reconstruction as verify_scenario, so the plot and the metrics agree.
    rec = afolu_reconstruction(sc)
    afolu_target = rec["afolu_target"]
    tgt_years = np.array([int(float(c)) for c in afolu_target.index])
    tgt_vals = afolu_target.values
    flux_years_src = rec["input_flux_years"]

    # Show the ramp plus the committed regrowth tail until it is ~negligible.
    yrs = rec["ext_years"]
    active = (rec["ramp"] > 1e-8) | (np.abs(rec["ext_calib"]) > 10.0)
    last = int(yrs[np.where(active)[0][-1]]) if np.any(active) else int(ext_years[-1])
    m = yrs <= last

    fig3, (ax3, ax3b) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                     gridspec_kw={'height_ratios': [2, 1]})
    mask = (tgt_years >= flux_years_src[0]) & (tgt_years <= last)
    ax3.plot(tgt_years[mask], tgt_vals[mask], 'k-', lw=1.5, label='IAM target')
    ax3.plot(flux_years_src, rec["calibrated_pre"], color='steelblue', lw=1.5, ls='--',
             label='Reconstructed (source)')
    ax3.plot(yrs[m], rec["ext_calib"][m], color='steelblue', lw=2,
             label='Reconstructed (extension)')
    ax3.axvline(2100, color='gray', ls='--', lw=0.8)
    ax3.axvline(int(ext_years[-1]), color='gray', ls=':', lw=0.8)
    ax3.axhline(0, color='gray', lw=0.4)
    ax3.set_ylabel('AFOLU CO₂ (Mt CO₂/yr)')
    ax3.set_title(f'{sc.key} — {sc.model}: AFOLU flux — IAM target vs reconstruction')
    ax3.legend()
    ax3.ticklabel_format(useOffset=False)

    # Extension components
    ax3b.plot(yrs[m], rec["baseline_term"][m], lw=1.5, label='Baseline (β+γ·t₂₁₀₀)·r(t)')
    ax3b.plot(yrs[m], rec["transition_term"][m], lw=1.5, label='Transition')
    ax3b.plot(yrs[m], rec["stock_term"][m], lw=1.5, label='Regrowth (stock)')
    ax3b.axhline(0, color='gray', lw=0.4)
    ax3b.axvline(2100, color='gray', ls='--', lw=0.8)
    ax3b.axvline(int(ext_years[-1]), color='gray', ls=':', lw=0.8)
    ax3b.set_xlabel('Year')
    ax3b.set_ylabel('Mt CO₂/yr')
    ax3b.set_title('Extension components', fontsize=10)
    ax3b.legend(fontsize=8)
    fig3.tight_layout()
    p3 = sc.output_dir / "diag_03_verification.png"
    fig3.savefig(p3, dpi=150, bbox_inches='tight')
    plt.close(fig3)
    log(f"  {sc.key}: saved {p3.name}")
