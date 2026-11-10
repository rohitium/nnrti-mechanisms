#!/usr/bin/env python3
"""Per-residue flexibility across the whole of p66.

Every other structural figure in this paper looks at the pocket. This one
looks at the rest of the enzyme, and it is where G190E separates itself from
the panel: it is the only genotype whose p66 fingers are resolvably more
mobile than wild type, and it is more mobile in five regions at once.

The curves are referenced to wild type residue by residue, in the bound
state, and are restricted to p66 -- the subdomain boundaries drawn along the
top are p66's, and averaging p66 with p51 at the same canonical number would
mix two different chains into one trace.

Per-residue secondary structure was drawn here as a second panel and dropped:
the WT-referenced difference is zero almost everywhere, and where it is not
it is one terminus or one frayed loop flickering between assignments, which
says nothing about resistance. The per-region summary stays in the
supplementary table.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_protein_level_figure
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

#: p66 subdomains, canonical HIV-1 RT numbering. Drawn as bands so a reader
#: can place a peak without counting along the axis. The fourth field is the
#: label row: the fingers insertion at 118-155 and the palm segment before it
#: are each only ~35 residues wide, so their labels collide on one line and
#: the narrower of the two is raised.
SUBDOMAINS = [("fingers", 1, 85, 0), ("palm", 86, 117, 0),
              ("fingers", 118, 155, 1), ("palm", 156, 236, 0),
              ("thumb", 237, 318, 0), ("connection", 319, 426, 0),
              ("RNase H", 427, 560, 0)]

HIGHLIGHT = "G190E"
C_HIGHLIGHT = "#C62828"
C_PANEL = "#1A1A1A"
C_OTHER = "#6E8FB2"

#: Genotypes left out of the panel. F227C is excluded so this figure covers
#: the same 18 mutants as the rest of the protein-level analysis.
EXCLUDE = ("F227C",)

#: The two regions the Results text discusses, marked so a reader can find
#: them. The primer grip is contiguous and gets a band; the NNIBP lining is
#: 15 scattered residues (canonical p66 numbering, as in
#: compute_protein_level_changes.NNIBP_P66) and gets a tick each.
PRIMER_GRIP = (225, 235)
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229,
             234, 318)
C_GRIP = "#E8A33D"
C_NNIBP = "#1B7F79"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 20, "axes.labelsize": 22, "axes.titlesize": 24,
        "xtick.labelsize": 19, "ytick.labelsize": 19,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


def wt_referenced(d: pd.DataFrame, col: str,
                  exclude: tuple[str, ...] = EXCLUDE) -> pd.DataFrame:
    """Per-residue genotype minus wild type, p66 only, bound state."""
    d = d[(d.structure_type == "holo") & (d.chain == 0)]
    m = d.groupby(["safe_label", "canonical"])[col].mean().unstack(0)
    if "wt" not in m.columns:
        raise SystemExit(f"no wild-type column for {col}")
    drop = ["wt"] + [g for g in exclude if g in m.columns]
    missing = [g for g in exclude if g not in m.columns]
    if missing:
        print(f"note: nothing to exclude for {missing}")
    return m.drop(columns=drop).sub(m["wt"], axis=0)


def draw(ax, delta: pd.DataFrame, ylabel: str):
    ax.axhline(0, color="0.55", lw=0.8)
    others = [c for c in delta.columns if c != HIGHLIGHT]
    for c in others:
        ax.plot(delta.index, delta[c], lw=0.6, alpha=0.35, color=C_OTHER,
                zorder=2)
    ax.plot(delta.index, delta.mean(axis=1), lw=1.8, color=C_PANEL, zorder=4)
    if HIGHLIGHT in delta.columns:
        ax.plot(delta.index, delta[HIGHLIGHT], lw=2.2, color=C_HIGHLIGHT,
                zorder=5)
    ax.set_ylabel(ylabel)
    ax.set_xlim(float(delta.index.min()), float(delta.index.max()))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def shade_subdomains(ax):
    """Alternating bands, so a peak can be placed without counting ticks."""
    for k, (_name, a, b, _row) in enumerate(SUBDOMAINS):
        if k % 2:
            ax.axvspan(a, b, color="0.92", zorder=0)


def highlight_regions(ax):
    """Mark the two regions the Results text names."""
    from matplotlib.transforms import blended_transform_factory

    a, b = PRIMER_GRIP
    ax.axvspan(a - 0.5, b + 0.5, color=C_GRIP, alpha=0.30, lw=0, zorder=1)
    # Ticks hang off the bottom spine in axes fraction, so they do not scale
    # with the y data and cannot be mistaken for a trace.
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for r in NNIBP_P66:
        ax.plot([r, r], [0.0, 0.05], color=C_NNIBP, lw=1.7, transform=trans,
                zorder=6, clip_on=False, solid_capstyle="butt")


def subdomain_bands(ax, lo: float, hi: float):
    """One label per p66 subdomain, above the top panel."""
    for name, a, b, row in SUBDOMAINS:
        if b < lo or a > hi:
            continue
        ax.annotate(name, xy=((max(a, lo) + min(b, hi)) / 2.0, 1.0),
                    xycoords=("data", "axes fraction"),
                    xytext=(0, 5 + 17 * row), textcoords="offset points",
                    ha="center", va="bottom", fontsize=15, color="0.3")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tables", type=Path,
                    default=Path("results/analysis/protein_level/tables"))
    ap.add_argument("--output-dir", type=Path,
                    default=Path("results/analysis/protein_level/plots"))
    args = ap.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plt = _style()
    rmsf = pd.read_csv(args.tables / "ca_rmsf_per_residue.csv")
    dr = wt_referenced(rmsf, "rmsf")

    fig, ax1 = plt.subplots(figsize=(13.0, 4.8))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.87, bottom=0.185)
    draw(ax1, dr, "Δ Cα RMSF (Å)")
    lo, hi = float(dr.index.min()), float(dr.index.max())
    shade_subdomains(ax1)
    highlight_regions(ax1)
    subdomain_bands(ax1, lo, hi)
    ax1.set_xlabel("RT p66 residues")

    handles = [plt.Line2D([], [], color=C_HIGHLIGHT, lw=2.2, label=HIGHLIGHT),
               plt.Line2D([], [], color=C_PANEL, lw=1.8, label="panel mean"),
               plt.Line2D([], [], color=C_OTHER, lw=0.9, alpha=0.6,
                          label="other genotypes"),
               plt.Line2D([], [], color=C_NNIBP, lw=2.6,
                          label="NNIBP lining"),
               plt.Rectangle((0, 0), 1, 1, color=C_GRIP, alpha=0.30,
                             label="primer grip 225\u2013235")]
    # Bottom left of the upper panel: the only corner no trace reaches.
    ax1.legend(handles=handles, frameon=False, fontsize=16, ncol=2,
               loc="lower left", handlelength=1.6, labelspacing=0.25,
               columnspacing=1.2)

    out = args.output_dir / "figure6_protein_level.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    print(f"wrote {out}  ({dr.shape[1]} genotypes, {dr.shape[0]} p66 residues)")
    print(f"{HIGHLIGHT} peak ΔRMSF {dr[HIGHLIGHT].max():.2f} A at residue "
          f"{int(dr[HIGHLIGHT].idxmax())}; panel mean peak "
          f"{dr.mean(axis=1).max():.2f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
