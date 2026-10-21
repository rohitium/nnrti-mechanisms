#!/usr/bin/env python3
"""Table-of-contents graphic for the manuscript.

ACS requires a TOC graphic on the last page, as .tif or .pdf, with no caption
and all text legible at the printed size. The guidelines put the usable area at
about 3.25 x 1.75 inches, so everything here is sized for that: the figure is
built at those dimensions rather than shrunk afterwards, which is what makes
labels illegible in practice.

The graphic has to say what the paper does in one glance. Left: doravirine in
the NNIBP with the panel's substituted positions around it, which is the
system. Right: the computed relative binding free energy against the measured
fold change in susceptibility, which is the result and its honest scatter.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_toc_graphic
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESID_OFFSET = -3
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
SITES = (98, 100, 103, 106, 181, 188, 190, 225, 227, 230, 234, 318)


def main(argv: list[str] | None = None) -> int:
    import MDAnalysis as mda
    from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", type=Path,
                   default=Path("results/md_runs/wt/rep_01/wt_minimized_rep01.pdb"))
    p.add_argument("--ddg", type=Path, default=Path("results/analysis/fep_pmx/panel_ddg.csv"))
    p.add_argument("--fold", type=Path,
                   default=Path("results/analysis/susceptibility/tables/dor_susceptibility_values.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("revisions"))
    p.add_argument("--dpi", type=int, default=600)
    args = p.parse_args(argv)

    u = mda.Universe(str(args.structure))
    seg = sorted({s.segid for s in u.segments})[0]
    p66 = u.select_atoms(f"protein and segid {seg}")
    dor = u.select_atoms("resname 2KW")
    pocket = p66.select_atoms("name CA and resid "
                             + " ".join(str(r + RESID_OFFSET) for r in NNIBP_P66))
    centre = pocket.positions.mean(axis=0)
    _, _, vt = np.linalg.svd(pocket.positions - centre)
    proj = lambda xyz: (np.asarray(xyz) - centre) @ vt[:2].T

    fig = plt.figure(figsize=(3.3, 1.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.08,
                          left=0.015, right=0.985, top=0.99, bottom=0.17)

    ax = fig.add_subplot(gs[0, 0])
    ca = p66.select_atoms("name CA")
    near = np.linalg.norm(ca.positions - centre, axis=1) < 19.0
    xy = proj(ca.positions)[near]
    ax.plot(xy[:, 0], xy[:, 1], "-", color="0.88", lw=0.5, zorder=1)
    d_xy = proj(dor.select_atoms("not name H*").positions)
    ax.scatter(d_xy[:, 0], d_xy[:, 1], s=13, color="#1a9850", zorder=4, linewidths=0)
    # One colour for every substituted position. Colouring them by phenotype
    # would be wrong: residues 106 and 190 appear in both susceptible and
    # resistant genotypes, so no single category belongs to the position.
    for site in SITES:
        sel = p66.select_atoms(f"resid {site + RESID_OFFSET} and name CA")
        if sel.n_atoms == 0:
            continue
        ax.scatter(*proj(sel.positions)[0], s=20, color="#f16913",
                   edgecolor="black", lw=0.3, zorder=5)

    ax.set_xlim(-16, 16); ax.set_ylim(-14, 14)
    ax.set_aspect("equal"); ax.axis("off")
    ax.text(0.5, 0.02, "doravirine in the NNIBP",
            transform=ax.transAxes, fontsize=5.0, va="bottom", ha="center", color="0.3")
    ax.text(0.5, 0.97, "18 clinical genotypes", transform=ax.transAxes,
            fontsize=5.0, va="top", ha="center", color="0.3")

    ax2 = fig.add_subplot(gs[0, 1])
    ddg = pd.read_csv(args.ddg)
    fold = pd.read_csv(args.fold)
    fcol = next((c for c in fold.columns if "fold" in c.lower()), None)
    gcol = next((c for c in fold.columns if "geno" in c.lower() or "mut" in c.lower()), None)
    m = ddg.merge(fold[[gcol, fcol]].rename(columns={gcol: "genotype", fcol: "fold"}),
                  on="genotype", how="left")
    m["fold"] = m["fold"].fillna(m["dor_fold_reduction"])
    m = m.dropna(subset=["fold", "ddg_bind_kcal"])
    m = m[m.fold > 0]
    colors = ["#d7301f" if g in RESISTANT else "#2c7fb8" if g in SUSCEPTIBLE else "#756bb1"
              for g in m.genotype]
    ax2.errorbar(np.log10(m.fold), m.ddg_bind_kcal, yerr=m.sem_kcal, fmt="none",
                 ecolor="0.6", elinewidth=0.5, capsize=0, zorder=2)
    ax2.scatter(np.log10(m.fold), m.ddg_bind_kcal, c=colors, s=16, zorder=3,
                edgecolor="black", linewidths=0.3)
    ax2.axhline(0.5, color="0.35", ls="--", lw=0.6)
    ax2.set_xlabel("log$_{10}$ fold resistance", fontsize=5.6, labelpad=1.5)
    ax2.set_ylabel("computed $\\Delta\\Delta G_{bind}$ (kcal/mol)", fontsize=5.6, labelpad=1.5)
    ax2.tick_params(labelsize=5.0, length=2, pad=1)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax2.spines[s].set_linewidth(0.6)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=3.2,
                          markeredgecolor="black", markeredgewidth=0.3, label=k)
               for k, c in (("resistant", "#d7301f"), ("susceptible", "#2c7fb8"),
                            ("uncertain", "#756bb1"))]
    ax2.legend(handles=handles, fontsize=4.4, loc="upper left", frameon=False,
               handletextpad=0.3, borderpad=0.1, labelspacing=0.25)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdf = args.output_dir / "TOC-Graphic.pdf"
    fig.savefig(pdf, dpi=args.dpi, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(args.output_dir / "TOC-Graphic.tif", dpi=args.dpi,
                bbox_inches="tight", pad_inches=0.02, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(args.output_dir / "TOC-Graphic.png", dpi=args.dpi,
                bbox_inches="tight", pad_inches=0.02)
    print(f"wrote {pdf} (+ .tif, .png) at {args.dpi} dpi, {m.shape[0]} genotypes plotted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
