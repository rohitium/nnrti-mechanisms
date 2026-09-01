#!/usr/bin/env python3
"""Figure 1B -- doravirine's moieties and the NNIBP residues that engage them.

Why this exists
---------------
The manuscript names DOR's moieties throughout (chlorocyanophenyl, pyridinone,
triazolinone) and the residues that contact them, but never shows the reader
which part of the molecule is which. This draws DOR with each moiety shaded and
labelled, and places a schematic of each key residue beside the moiety it
engages, annotated with the contact count measured from the wild-type
trajectories rather than assigned by eye.

Residues are positioned next to what they touch, so leader lines are
unnecessary; the only connector drawn is the Lys103 main-chain hydrogen bond,
which is the one interaction between specific atoms rather than a packing
contact.

Contacts are protein-ligand heavy-atom pairs within 4.0 A, averaged over three
100 ns wild-type replicates. See manuscript/contact-cutoff-sensitivity.md for
why 4.0 A.

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
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

#: residue -> (side-chain SMILES, contacts with whole ligand, caption, centre, scale)
#: Fragments stand in for the side chain (p-cresol for Tyr, toluene for Phe,
#: 3-methylindole for Trp, isobutane for Val); Lys103 is drawn as a main-chain
#: amide because it is the backbone, not the side chain, that binds DOR.
RESIDUES = {
    "Tyr188": ("Cc1ccc(O)cc1", 20.3, "stacks chlorocyanophenyl", (8.6, 7.4), 0.80),
    "Trp229": ("Cc1c[nH]c2ccccc12", 8.4, "", (14.0, 0.4), 0.80),
    "Phe227": ("Cc1ccccc1", 6.5, "", (-1.4, 7.5), 0.80),
    "Tyr318": ("Cc1ccc(O)cc1", 9.5, "", (-9.6, 5.2), 0.80),
    "Lys103": ("CC(=O)NC", 6.7, "main-chain C=O", (-10.6, -3.4), 0.85),
    "Val106": ("CC(C)C", 9.5, "hydrophobic packing", (-1.6, -7.4), 0.85),
    "Tyr181": ("Cc1ccc(O)cc1", 4.3, "rotated away; no ring contact", (-8.2, -7.6), 0.72),
    "Gly190": (None, 3.0, "no side chain", (5.6, -7.2), 0.72),
}
FAINT = {"Tyr181", "Gly190"}

MOIETY_COLOR = {
    "chlorocyanophenyl": "#3B6EA8",
    "pyridinone": "#B5474B",
    "triazolinone": "#3F8C5E",
}
MOIETY_LABEL = {
    "chlorocyanophenyl": "chlorocyanophenyl",
    "pyridinone": "pyridinone\n(central)",
    "triazolinone": "triazolinone",
}
MOIETY_LABEL_XY = {
    "triazolinone": (-6.4, 2.4),
    "pyridinone": (-3.2, -4.1),
    "chlorocyanophenyl": (7.4, 3.1),
}
HETERO_COLOR = {"N": "#1f4e9c", "O": "#c0392b", "F": "#7d3c98", "Cl": "#1e8449"}


def classify_rings(mol) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for ring in mol.GetRingInfo().AtomRings():
        nbrs = {n.GetSymbol() for a in ring
                for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in ring}
        if "Cl" in nbrs:
            out["chlorocyanophenyl"] = list(ring)
        elif len(ring) == 5:
            out["triazolinone"] = list(ring)
        else:
            out["pyridinone"] = list(ring)
    return out


def draw_mol(ax, mol, pos, *, colour_fn, lw=2.4, hetero_size=9.5, zorder=2,
             label_hetero=True, dot=170, dbl_offset=0.11):
    """Draw a 2D molecule as line bonds with heteroatom labels."""
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        col = colour_fn(i, j)
        p, q = pos[i], pos[j]
        if b.GetBondTypeAsDouble() == 2.0:
            d = q - p
            n = np.array([-d[1], d[0]])
            n = n / (np.linalg.norm(n) + 1e-9) * dbl_offset
            for s in (+1, -1):
                ax.plot(*zip(p + n * s, q + n * s), color=col, lw=lw,
                        zorder=zorder, solid_capstyle="round")
        else:
            ax.plot(*zip(p, q), color=col, lw=lw, zorder=zorder,
                    solid_capstyle="round")
    if not label_hetero:
        return
    for i, a in enumerate(mol.GetAtoms()):
        s = a.GetSymbol()
        if s == "C":
            continue
        ax.scatter(*pos[i], s=dot, color="white", zorder=zorder + 1, edgecolors="none")
        ax.text(*pos[i], s, ha="center", va="center", fontsize=hetero_size,
                fontweight="bold", zorder=zorder + 2,
                color=HETERO_COLOR.get(s, "#333333"))


def fragment_pos(smiles: str, centre, scale: float):
    frag = Chem.MolFromSmiles(smiles)
    # Kekulize: aromatic bonds report GetBondTypeAsDouble() == 1.5, which the
    # renderer would draw as a plain single line, making Tyr/Phe/Trp look
    # saturated. Kekulizing gives explicit alternating single/double bonds.
    Chem.Kekulize(frag, clearAromaticFlags=True)
    AllChem.Compute2DCoords(frag)
    c = frag.GetConformer()
    p = np.array([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y]
                  for i in range(frag.GetNumAtoms())])
    p = (p - p.mean(axis=0)) * scale + np.asarray(centre, dtype=float)
    return frag, p


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sdf", type=Path, default=repo / "data/ligands/dor.sdf")
    ap.add_argument("--output", type=Path,
                    default=repo / "results/plots/figure1B_dor_schematic.pdf")
    args = ap.parse_args()

    mol = Chem.RemoveHs(Chem.MolFromMolFile(str(args.sdf)))
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    pos = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                    for i in range(mol.GetNumAtoms())])

    rings = classify_rings(mol)
    atom_moiety = {a: m for m, atoms in rings.items() for a in atoms}
    ring_centroid = {m: pos[a].mean(axis=0) for m, a in rings.items()}

    fig, ax = plt.subplots(figsize=(12.0, 8.4))

    # --- doravirine
    for m, c in ring_centroid.items():
        ax.scatter(*c, s=6400, color=MOIETY_COLOR[m], alpha=0.12, zorder=0)

    def dor_colour(i, j):
        mi, mj = atom_moiety.get(i), atom_moiety.get(j)
        return MOIETY_COLOR[mi] if (mi and mi == mj) else "#3d3d3d"

    draw_mol(ax, mol, pos, colour_fn=dor_colour, lw=2.7, hetero_size=10.5,
             zorder=2, dot=210, dbl_offset=0.13)

    for m in rings:
        lx, ly = MOIETY_LABEL_XY[m]
        ax.text(lx, ly, MOIETY_LABEL[m], ha="center", va="center", fontsize=12.5,
                fontweight="bold", color=MOIETY_COLOR[m], zorder=5)

    # --- residues
    for res, (smi, n, note, centre, scale) in RESIDUES.items():
        faint = res in FAINT
        col = "#9aa0a6" if faint else "#2f2f2f"
        cx, cy = centre
        if smi is not None:
            frag, fp = fragment_pos(smi, centre, scale)
            draw_mol(ax, frag, fp, colour_fn=lambda i, j, c=col: c, lw=1.9,
                     hetero_size=8.6, zorder=3, dot=150,
                     dbl_offset=0.115 * scale)
            ylo = fp[:, 1].min()
        else:  # glycine: no side chain to draw
            ax.text(cx, cy, "H", ha="center", va="center", fontsize=13,
                    fontweight="bold", color=col, zorder=4)
            ylo = cy - 0.35
        head = f"{res}   {n:.1f} contacts" if not faint else f"{res}   {n:.1f}"
        ax.text(cx, ylo - 0.72, head, ha="center", va="top", fontsize=11,
                fontweight="bold", color=col, zorder=6)
        if note:
            ax.text(cx, ylo - 1.42, note, ha="center", va="top", fontsize=8.8,
                    style="italic", color="#6b6b6b" if not faint else "#a5a5a5",
                    zorder=6)

    # --- the one connector: triazolinone N-H ... Lys103 main-chain C=O
    tri = rings["triazolinone"]
    nh = [k for k in tri if mol.GetAtomWithIdx(k).GetSymbol() == "N"
          and mol.GetAtomWithIdx(k).GetTotalNumHs() > 0]
    if nh:
        p = pos[nh[0]]
        lys_c = np.array(RESIDUES["Lys103"][3], dtype=float)
        v = lys_c - p
        v = v / np.linalg.norm(v)
        start, end = p + v * 0.45, lys_c - v * 1.55
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="-", color="#c0392b", lw=2.0,
                                    linestyle=(0, (3, 2.4))), zorder=4)
        ax.scatter(*p, s=250, facecolor="none", edgecolor="#c0392b", lw=1.8,
                   zorder=6)
        mid = (start + end) / 2
        ax.text(mid[0], mid[1] + 0.55, "N–H···O=C\n3.0 Å", ha="center",
                va="bottom", fontsize=9.2, color="#c0392b", fontweight="bold",
                zorder=7,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="none", alpha=0.9))

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-15.5, 19.0)
    ax.set_ylim(-11.0, 10.0)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
