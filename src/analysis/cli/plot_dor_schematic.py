#!/usr/bin/env python3
"""Figure 1B -- annotated schematic of doravirine and its NNIBP environment.

Why this exists
---------------
The manuscript refers to DOR's moieties by name throughout (chlorocyanophenyl,
pyridinone, triazolinone) and to the residues that engage them, but never shows
the reader which part of the molecule is which. This draws the 2D structure with
each moiety shaded and labelled, and places the key NNIBP residues beside the
moiety they actually contact, with the contact counts measured from the
wild-type trajectories rather than assigned by eye.

Contact counts are atom pairs within 4.5 A, averaged over three 100 ns WT
replicates (see compute_dor_moiety_contacts.py for the same convention).

Usage
-----
    PYTHONPATH=. python -m src.analysis.cli.plot_dor_schematic
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# WT atom-pair contacts (<4.5 A) per moiety, from the equilibrium trajectories.
# residue -> (moiety, contacts, note)
# residue -> (moiety, contacts, note, label_xy)
# Label positions are hand-placed against the deterministic RDKit 2D layout
# (ring centroids: triazolinone (-5.00, 0.26), pyridinone (-0.38, -1.66),
# chlorocyanophenyl (3.67, 1.59)) so no leader line crosses the structure.
RESIDUE_CONTACTS = {
    "Tyr188": ("chlorocyanophenyl", 27.5, "aromatic stacking", (11.2, 2.4)),
    "Trp229": ("chlorocyanophenyl", 16.9, "", (10.0, 6.0)),
    "Phe227": ("chlorocyanophenyl", 9.0, "", (4.6, 7.6)),
    "Tyr181": ("chlorocyanophenyl", 0.4, "rotated away;\nno contact with DOR", (-1.2, 7.4)),
    "Tyr318": ("triazolinone", 18.6, "", (-11.6, 4.2)),
    "Lys103": ("triazolinone", 10.9, "main-chain C=O accepts\nN–H···O=C, 3.0 Å", (-11.8, -3.2)),
    "Val106": ("pyridinone", 8.9, "hydrophobic packing", (-4.2, -7.8)),
    "Gly190": ("pyridinone", 0.1, "no direct contact", (4.4, -7.2)),
}

# moiety -> label anchor
MOIETY_LABEL_XY = {
    "triazolinone": (-6.6, 2.3),
    "pyridinone": (-3.0, -4.0),
    "chlorocyanophenyl": (7.4, 0.0),
}

MOIETY_COLOR = {
    "chlorocyanophenyl": "#4C72B0",
    "pyridinone": "#C44E52",
    "triazolinone": "#55A868",
}
MOIETY_LABEL = {
    "chlorocyanophenyl": "chlorocyanophenyl",
    "pyridinone": "pyridinone\n(central)",
    "triazolinone": "triazolinone",
}


def classify_rings(mol) -> dict[str, list[int]]:
    """Assign each ring of DOR to a moiety by its substituents."""
    out: dict[str, list[int]] = {}
    for ring in mol.GetRingInfo().AtomRings():
        nbrs = {
            n.GetSymbol()
            for a in ring
            for n in mol.GetAtomWithIdx(a).GetNeighbors()
            if n.GetIdx() not in ring
        }
        if "Cl" in nbrs:
            out["chlorocyanophenyl"] = list(ring)
        elif len(ring) == 5:
            out["triazolinone"] = list(ring)
        else:
            out["pyridinone"] = list(ring)
    return out


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sdf", type=Path, default=repo / "data/ligands/dor.sdf")
    ap.add_argument("--output", type=Path,
                    default=repo / "results/plots/figure1B_dor_schematic.pdf")
    args = ap.parse_args()

    mol = Chem.MolFromMolFile(str(args.sdf))
    if mol is None:
        raise SystemExit(f"could not read {args.sdf}")
    mol = Chem.RemoveHs(mol)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    pos = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                    for i in range(mol.GetNumAtoms())])

    rings = classify_rings(mol)
    atom_moiety = {a: m for m, atoms in rings.items() for a in atoms}

    fig, ax = plt.subplots(figsize=(12.5, 8.2))

    # bonds
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        mi, mj = atom_moiety.get(i), atom_moiety.get(j)
        col = MOIETY_COLOR[mi] if (mi and mi == mj) else "#444444"
        lw = 2.6 if (mi and mi == mj) else 1.9
        p, q = pos[i], pos[j]
        if b.GetBondTypeAsDouble() == 2.0:
            d = q - p
            n = np.array([-d[1], d[0]])
            n = n / (np.linalg.norm(n) + 1e-9) * 0.055
            for s in (+1, -1):
                ax.plot(*zip(p + n * s, q + n * s), color=col, lw=lw, zorder=2,
                        solid_capstyle="round")
        else:
            ax.plot(*zip(p, q), color=col, lw=lw, zorder=2, solid_capstyle="round")

    # heteroatom labels
    for i, a in enumerate(mol.GetAtoms()):
        s = a.GetSymbol()
        if s == "C":
            continue
        ax.scatter(*pos[i], s=210, color="white", zorder=3, edgecolors="none")
        ax.text(*pos[i], s, ha="center", va="center", fontsize=10.5,
                fontweight="bold", zorder=4,
                color={"N": "#1f4e9c", "O": "#c0392b", "F": "#7d3c98",
                       "Cl": "#1e8449"}.get(s, "#333333"))

    # moiety shading + label
    centre = pos.mean(axis=0)
    ring_centroid = {m: pos[a].mean(axis=0) for m, a in rings.items()}
    for m, c in ring_centroid.items():
        ax.scatter(*c, s=6400, color=MOIETY_COLOR[m], alpha=0.11, zorder=0)
        lx, ly = MOIETY_LABEL_XY[m]
        ax.text(lx, ly, MOIETY_LABEL[m], ha="center", va="center",
                fontsize=12, fontweight="bold", color=MOIETY_COLOR[m], zorder=5)

    # residue callouts at hand-placed anchors
    for res, (moi, n, note, xy) in RESIDUE_CONTACTS.items():
        c = ring_centroid[moi]
        tx, ty = xy
        faint = n < 1.0
        col = "#9aa0a6" if faint else MOIETY_COLOR[moi]
        is_hbond = res == "Lys103"
        ax.annotate("", xy=c, xytext=(tx, ty),
                    arrowprops=dict(arrowstyle="-",
                                    color="#c0392b" if is_hbond else col,
                                    lw=1.7 if is_hbond else 1.3,
                                    linestyle=(0, (3, 2)) if is_hbond
                                    else (":" if faint else "-"),
                                    alpha=0.85, shrinkA=10, shrinkB=34),
                    zorder=1)
        head = res if faint else f"{res}  ·  {n:.1f} contacts"
        if note:
            head += f"\n{note}"
        ax.text(tx, ty, head, ha="center", va="center", fontsize=9.4,
                linespacing=1.4,
                color="#6b7075" if faint else "#1a1a1a", zorder=6,
                bbox=dict(boxstyle="round,pad=0.5",
                          facecolor="#f5f5f5" if faint else "white",
                          edgecolor="#c0392b" if is_hbond else col,
                          linewidth=1.6 if is_hbond else 1.2, alpha=0.98))

    # mark the donor nitrogen itself so the dashed Lys103 leader has an origin
    tri = rings["triazolinone"]
    nh = [k for k in tri if mol.GetAtomWithIdx(k).GetSymbol() == "N"
          and mol.GetAtomWithIdx(k).GetTotalNumHs() > 0]
    if nh:
        ax.scatter(*pos[nh[0]], s=250, facecolor="none", edgecolor="#c0392b",
                   lw=1.8, zorder=5)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-15.0, 16.0)
    ax.set_ylim(-10.5, 9.5)
    ax.set_title(
        "Doravirine moieties and their NNIBP contacts",
        fontsize=13, fontweight="bold", pad=12,
    )
    ax.text(0.5, -0.055,
            "Contacts are protein–ligand heavy-atom pairs within 4.5 Å, "
            "averaged over three 100 ns wild-type replicates.",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=8.4, color="#666666")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.output.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
