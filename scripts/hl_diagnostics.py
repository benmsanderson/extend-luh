#!/usr/bin/env python3
"""Diagnostic figures for the HL (WITCH) extension delivery.

Usage:  python scripts/hl_diagnostics.py <out_dir>

Writes five PNGs into <out_dir>:
  01_afolu_fit.png          AFOLU target vs emulator fit / extension, and ramp r(t)
  02_managed_land.png       global managed-land areas, input + extension
  03_plantation_sink.png    WITCH net AFOLU vs plantation afforestation removal
  04_overfull_cells.png     cells whose 2100 managed-land sum exceeds 100 %
  05_biofuel_beccs.png      BECCS scaling factor and 2nd-gen biofuel at the cap
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from src import carbon_model as cm
from src import config as cfg
from src.data import beccs_scaling_factors, filter_scenario, get_variable, load_csv, load_states
from src.pipeline import _cell_area

# ── Chart style (light surface, validated 3-slot categorical palette) ──
SURFACE, INK, INK2, MUTED, GRID, BASE = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"   # blue, orange, aqua (fixed order)
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "lines.linewidth": 1.6, "legend.frameon": False, "legend.fontsize": 8,
})


def label_end(ax, x, y, text, color_ink=INK2, dy=0):
    """Direct label at a line end, in text ink (the adjacent mark carries identity)."""
    ax.annotate(text, (x, y), xytext=(4, dy), textcoords="offset points",
                va="center", fontsize=8, color=color_ink)


def vline(ax, x, text):
    ax.axvline(x, color=BASE, lw=0.8, zorder=0)
    ax.text(x, 1.0, f" {text}", transform=ax.get_xaxis_transform(),
            va="top", ha="left", fontsize=7.5, color=MUTED)


def main(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    sc = cfg.get_scenario("HL")
    ds = load_states(scen=sc)
    lat, lon = ds.lat.values, ds.lon.values
    area = _cell_area(lat, lon)
    cal = np.load(sc.output_dir / "afolu_calibration.npz")
    df = load_csv()
    sub = filter_scenario(df, scen=sc)
    afolu = get_variable(sub, "Emissions|CO2|AFOLU")

    # ── 01 AFOLU fit + ramp ────────────────────────────────────────
    years = ds.time.values.astype(int)
    fy = years[1:]
    incr = {"othr": np.array([
        np.nansum((ds["othr"].sel(time=y).values - ds["othr"].sel(time=y - 1).values) * area)
        for y in fy])}
    tau = float(cal["tau"])
    sf = cm.stock_flux_for_tau(tau, incr, regrow_vars=["othr"])
    tt = (fy - fy.mean()) / 10.0
    fitted = float(cal["beta"]) + float(cal["gamma"]) * tt + float(cal["alpha_stock"]) * sf
    ey, rec, ramp = cal["years"], cal["afolu_reconstructed"], cal["ramp"]
    xmax = 2250

    fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True,
                               gridspec_kw={"height_ratios": [2.2, 1]})
    iy = np.arange(2025, xmax + 1)
    a.plot(iy, afolu.reindex(iy).values / 1000, color=S1)
    a.plot(fy, fitted / 1000, color=S2)
    m = ey <= xmax
    a.plot(ey[m], rec[m] / 1000, color=S3)
    a.axhline(0, color=BASE, lw=0.8, zorder=0)
    label_end(a, xmax, afolu.loc[xmax] / 1000, "IAM target", dy=7)
    label_end(a, fy[-1], fitted[-1] / 1000, "emulator fit", dy=-8)
    label_end(a, xmax, rec[m][-1] / 1000, "extension (emulated)", dy=-7)
    vline(a, 2100, "2100")
    vline(a, cfg.YR_AFOLU_ZERO, f"{cfg.YR_AFOLU_ZERO}: land-use change stops")
    a.set_ylabel("AFOLU CO$_2$ (Gt CO$_2$/yr)")
    a.set_title("HL AFOLU: IAM target, calibrated emulator and extension")
    a.legend(handles=a.lines[:3], labels=["IAM target", "emulator fit", "extension (emulated)"],
             loc="upper right")
    b.plot(ey[m], ramp[m], color=S1)
    vline(b, cfg.YR_AFOLU_ZERO, "")
    b.set_ylim(0, 1.05)
    b.set_ylabel("ramp r(t)")
    b.set_title("Land-use change rate as a multiple of the 2100 rate")
    b.set_xlim(2025, xmax + 40)
    fig.tight_layout()
    fig.savefig(out_dir / "01_afolu_fit.png", dpi=150)
    plt.close(fig)

    # ── 02 Managed-land areas (input + extension), small multiples ──
    ext = {}
    for stem, gv in sc.state_output_groups.items():
        f = next(sc.output_dir.glob(f"{stem}_extension_transient_*.nc"))
        e = xr.open_dataset(f)
        for v in gv:
            ext[v] = e[v]
    panels = [("cropland", "Cropland"), ("pastr", "Pasture"), ("range", "Rangeland"), ("urban", "Urban")]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True)
    for ax, (v, name) in zip(axes.flat, panels):
        gi = [(ds[v].sel(time=y).fillna(0).values * area).sum() / 1e6 for y in years]
        et = ext[v].time.values.astype(int)
        ge = [(ext[v].sel(time=y).fillna(0).values * area).sum() / 1e6 for y in et]
        ax.plot(np.r_[years, et], np.r_[gi, ge], color=S1)
        vline(ax, 2100, "extension →")
        ax.set_title(name)
        ax.set_ylabel("Mkm$^2$")
        lo, hi = min(gi + ge), max(gi + ge)
        if hi <= 0:
            ax.set_ylim(-0.1, 1.0)
            ax.text(0.5, 0.5, "zero in every year of the WITCH input", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color=INK2)
            continue
        pad = max(0.05 * (hi - lo), 0.02 * hi)
        ax.set_ylim(lo - pad, hi + pad)
    fig.suptitle("HL global managed-land area: WITCH input (2024–2100) and extension",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_dir / "02_managed_land.png", dpi=150)
    plt.close(fig)

    # ── 03 Plantation sink vs net AFOLU (WITCH World) ───────────────
    w = pd.read_csv(sc.data_dir / "WITCH_6.0_SSP5_-_Medium-Low_Emissions_a_32_afolu_variables.csv")
    w = w[w.region == "World"].set_index("variable")
    ycols = [c for c in w.columns if c.isdigit()]
    yy = np.array([int(c) for c in ycols])
    keep = yy >= 2020
    lu = -w.loc["Carbon Removal|Land Use", ycols].astype(float).values / 1000
    pl = -w.loc["Carbon Removal|Land Use|Re/Afforestation|Plantation", ycols].astype(float).values / 1000
    net = afolu.reindex(yy).values / 1000
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.plot(yy[keep], net[keep], color=S1, marker="o", ms=3.5)
    ax.plot(yy[keep], pl[keep], color=S2, marker="o", ms=3.5)
    ax.plot(yy[keep], lu[keep], color=S3, marker="o", ms=3.5)
    ax.axhline(0, color=BASE, lw=0.8, zorder=0)
    label_end(ax, 2100, net[-1], "net AFOLU CO$_2$", dy=10)
    label_end(ax, 2100, pl[-1], "plantation afforestation", dy=-4)
    label_end(ax, 2100, lu[-1], "all land-use CDR", dy=-6)
    ax.set_xlim(2020, 2130)
    ax.set_ylabel("Gt CO$_2$/yr (removals negative)")
    ax.set_title("WITCH: the late-century AFOLU sink is plantation afforestation")
    ax.legend(handles=ax.lines[:3], labels=["net AFOLU CO$_2$", "plantation afforestation",
                                            "all land-use CDR"], loc="lower left")
    fig.tight_layout()
    fig.savefig(out_dir / "03_plantation_sink.png", dpi=150)
    plt.close(fig)

    # ── 04 Overfull cells at 2100 ───────────────────────────────────
    man = sum(ds[v].sel(time=2100).fillna(0).values for v in sc.residual_from)
    excess = np.where(man > 1.001, man - 1.0, np.nan)
    land = np.where(man > 0, 1.0, np.nan)
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.grid(False)
    ext_ = [lon.min(), lon.max() + 0.25, lat.min(), lat.max() + 0.25]
    ax.imshow(land, extent=ext_, origin="upper", cmap=matplotlib.colors.ListedColormap(["#ebeae6"]),
              interpolation="nearest")
    cm_seq = matplotlib.colors.ListedColormap(SEQ[2:])
    im = ax.imshow(excess, extent=ext_, origin="upper", cmap=cm_seq, vmin=0, vmax=0.25,
                   interpolation="nearest")
    ax.set_ylim(-60, 85)
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    n = int(np.isfinite(excess).sum())
    ax.set_title(f"Managed land over 100 % of cell at 2100 ({n:,} cells)")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("excess over 100 % of cell", color=INK2)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "04_overfull_cells.png", dpi=150)
    plt.close(fig)

    # ── 05 Biofuel: BECCS factor and cells at the cap ───────────────
    bec = beccs_scaling_factors(sub)
    by = np.arange(2100, 2301)
    bf = xr.open_dataset(next(sc.output_dir.glob("biof_2nd_gen_*_extension_transient_*.nc")))["crpbf_c3per"]
    cy = np.arange(2101, 2301)
    share = []
    for y in cy:
        s = bf.sel(time=y).values
        pos = s > 0
        share.append(100 * np.sum(s[pos] >= 0.999) / max(pos.sum(), 1))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    a.plot(by, bec.reindex(by).values, color=S1)
    a.axhline(1, color=BASE, lw=0.8, zorder=0)
    a.set_title("BECCS scaling factor (× 2100)")
    b.plot(cy, share, color=S1)
    b.set_ylim(0, 100)
    b.set_title("crpbf_c3per cells at the 1.0 cap (%)")
    fig.tight_layout()
    fig.savefig(out_dir / "05_biofuel_beccs.png", dpi=150)
    plt.close(fig)
    print(f"wrote 5 figures to {out_dir}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else cfg.get_scenario("HL").output_dir / "figures")
