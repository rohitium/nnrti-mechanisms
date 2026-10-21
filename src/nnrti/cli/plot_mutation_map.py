#!/usr/bin/env python3
"""Map every panel mutation onto the RT-DOR structure.

Reviewer 1 (comment 4) asked for a figure placing all the panel's mutations on
the structure, showing where they sit relative to DOR and to the NNIBP.

Two panels:

  A. The p66 backbone projected onto the plane that best spreads the
     NNIBP-lining residues, so the pocket is seen face-on rather than down an
     arbitrary crystal axis. DOR is drawn in place; each mutated position is a
     marker coloured by its phenotype category.

  B. The minimum heavy-atom distance from each mutated residue to DOR, which
     turns the spatial picture into the number that matters: whether the
     substitution can touch the ligand at all. Residues beyond about 5 A cannot
     contact DOR directly in the crystallographic pose, so any effect they have
     is transmitted through the protein.

Residue numbering follows the project convention: topology resid = canonical
HXB2 numbering - 3.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_mutation_map
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN

RESID_OFFSET = -3

#: NNIBP-lining residues, p66, canonical numbering (Cilento, Kirby & Sarafianos)
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)

#: every substituted position in the panel, canonical numbering
SITES = {98: "A98G", 100: "L100I", 103: "K103N", 106: "V106A/I/M", 181: "Y181C",
         188: "Y188L", 190: "G190A/E/S", 225: "P225H", 227: "F227C/L",
         230: "M230L", 234: "L234I", 318: "Y318F"}

CATEGORY_COLOR = {"Susceptible": "#2c7fb8", "Resistant": "#d7301f",
                  "Uncertain": "#756bb1", "Mixed": "#525252"}


def site_category(site: int) -> str:
    """Category of a position, from the genotypes in which it appears."""
    cats = set()
    for group, name in ((SUSCEPTIBLE, "Susceptible"), (RESISTANT, "Resistant"),
                        (UNCERTAIN, "Uncertain")):
        for geno in group:
            for part in geno.split("+"):
                digits = "".join(c for c in part if c.isdigit())
                if digits and int(digits) == site:
                    cats.add(name)
    if not cats:
        return "Mixed"
    if len(cats) == 1:
        return cats.pop()
    return "Mixed"


def main(argv: list[str] | None = None) -> int:
    import MDAnalysis as mda

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", type=Path,
                   default=Path("results/md_runs/wt/rep_01/wt_minimized_rep01.pdb"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/mechanisms/plots"))
    p.add_argument("--ligand", default="2KW")
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)

    u = mda.Universe(str(args.structure))
    segs = sorted({s.segid for s in u.segments})
    p66 = u.select_atoms(f"protein and segid {segs[0]}")
    dor = u.select_atoms(f"resname {args.ligand}")
    if dor.n_atoms == 0:
        raise SystemExit(f"no ligand {args.ligand} in {args.structure}")

    pocket_resids = [r + RESID_OFFSET for r in NNIBP_P66]
    pocket = p66.select_atoms("name CA and resid " + " ".join(map(str, pocket_resids)))

    # Plane that best spreads the pocket residues: face-on view of the NNIBP.
    centre = pocket.positions.mean(axis=0)
    _, _, vt = np.linalg.svd(pocket.positions - centre)
    axes = vt[:2]
    proj = lambda xyz: (np.asarray(xyz) - centre) @ axes.T

    ca = p66.select_atoms("name CA")
    ca_xy = proj(ca.positions)
    dor_xy = proj(dor.positions)
    pocket_xy = proj(pocket.positions)

    fig = plt.figure(figsize=(13.5, 6.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.25)
    ax = fig.add_subplot(gs[0, 0])

    # Only the neighbourhood of the pocket: the whole p66 trace projected onto
    # this plane is an unreadable tangle and carries no information here.
    near = np.linalg.norm(ca.positions - centre, axis=1) < 24.0
    seg = ca_xy[near]
    ax.plot(seg[:, 0], seg[:, 1], "-", color="0.86", lw=0.6, zorder=1)
    ax.scatter(seg[:, 0], seg[:, 1], s=4, color="0.72", zorder=1)
    ax.scatter(pocket_xy[:, 0], pocket_xy[:, 1], s=90, facecolor="none",
               edgecolor="0.45", lw=1.0, zorder=2, label="NNIBP-lining residues")
    ax.scatter(dor_xy[:, 0], dor_xy[:, 1], s=26, color="#1a9850", zorder=4,
               label="doravirine")

    rows, marks = [], []
    for site, label in SITES.items():
        sel = p66.select_atoms(f"resid {site + RESID_OFFSET}")
        if sel.n_atoms == 0:
            continue
        heavy = sel.select_atoms("not name H*")
        d = float(np.min(np.linalg.norm(
            heavy.positions[:, None, :] - dor.select_atoms("not name H*").positions[None, :, :],
            axis=-1)))
        cat = site_category(site)
        rows.append((site, label, d, cat))

        xy = proj(sel.select_atoms("name CA").positions)[0]
        ax.scatter(*xy, s=120, color=CATEGORY_COLOR[cat], edgecolor="black",
                   lw=0.6, zorder=5)
        marks.append((xy, label))

    # Fan the labels out around the pocket so they do not sit on each other:
    # each marker gets a leader line to a label placed on a ring by angle.
    if marks:
        pts = np.array([m[0] for m in marks])
        hub = pts.mean(axis=0)
        ang = np.arctan2(pts[:, 1] - hub[1], pts[:, 0] - hub[0])
        order = np.argsort(ang)
        radius = 1.45 * max(6.0, np.abs(pts - hub).max())
        slots = np.linspace(-np.pi, np.pi, len(marks), endpoint=False)
        for slot, idx in zip(slots, order):
            xy, label = marks[idx]
            lx = hub[0] + radius * np.cos(slot)
            ly = hub[1] + radius * np.sin(slot)
            ax.annotate(label, xy=xy, xytext=(lx, ly), fontsize=7.5, fontweight="bold",
                        ha="center" if abs(np.cos(slot)) < 0.4 else ("left" if np.cos(slot) > 0 else "right"),
                        va="center",
                        arrowprops=dict(arrowstyle="-", color="0.55", lw=0.6,
                                        shrinkA=0, shrinkB=4))
        ax.set_xlim(hub[0] - radius * 1.55, hub[0] + radius * 1.55)
        ax.set_ylim(hub[1] - radius * 1.45, hub[1] + radius * 1.45)

    ax.set_xlabel("projection onto NNIBP plane, axis 1 (A)", fontsize=9)
    ax.set_ylabel("axis 2 (A)", fontsize=9)
    ax.set_title("A.  Panel mutations on the RT (p66) structure", fontsize=10, loc="left")
    ax.set_aspect("equal")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, markeredgecolor="k",
                          markersize=8, label=k) for k, c in CATEGORY_COLOR.items()]
    handles += [plt.Line2D([], [], marker="o", ls="", color="#1a9850", markersize=6,
                           label="doravirine"),
                plt.Line2D([], [], marker="o", ls="", markerfacecolor="none",
                           markeredgecolor="0.45", markersize=9, label="NNIBP lining")]
    ax.legend(handles=handles, fontsize=7.5, loc="upper left", frameon=False)

    ax2 = fig.add_subplot(gs[0, 1])
    rows.sort(key=lambda r: r[2])
    ypos = np.arange(len(rows))
    ax2.barh(ypos, [r[2] for r in rows],
             color=[CATEGORY_COLOR[r[3]] for r in rows], edgecolor="black", lw=0.5)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([r[1] for r in rows], fontsize=8)
    ax2.axvline(4.0, color="0.3", ls="--", lw=1.0)
    ax2.text(4.1, -0.45, "4 A contact", fontsize=7.5, color="0.3")
    ax2.set_xlabel("minimum heavy-atom distance to DOR (A)", fontsize=9)
    ax2.set_title("B.  Reach of each substituted position", fontsize=10, loc="left")
    ax2.invert_yaxis()
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    fig.tight_layout()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(args.output_dir / f"mutation_map.{ext}", dpi=args.dpi, bbox_inches="tight")

    print(f"wrote {args.output_dir}/mutation_map.png|pdf\n")
    print(f"{'site':12s} {'dist to DOR (A)':>16s}  category")
    for site, label, d, cat in rows:
        print(f"{label:12s} {d:16.2f}  {cat}")
    n_far = sum(1 for r in rows if r[2] > 5.0)
    print(f"\n{n_far} of {len(rows)} substituted positions lie more than 5 A from DOR "
          "and cannot contact it directly in the crystallographic pose.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
