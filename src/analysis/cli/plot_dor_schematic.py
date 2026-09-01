#!/usr/bin/env python3
"""Figure 1B -- doravirine's moieties and its three defining NNIBP interactions.

What is drawn
-------------
Tyr188 pi-stacking on the chlorocyanophenyl ring, Val106 packing against the
central pyridinone, and the Lys103 main-chain carbonyl accepting a hydrogen bond
from the triazolinone N-H. Each residue is drawn with its backbone (N-CA-C=O)
as well as its side chain.

How faithful it is
------------------
This is a **schematic**, not a projection. DOR is laid out as a butterfly -- the
pyridinone ring as the body, the chlorocyanophenyl and triazolinone rings as the
two wings -- and each residue is placed off the group it contacts, at a position
chosen for legibility rather than measured from the structure.

What is real: the connectivity and stereochemistry of DOR, each residue's own
atoms drawn face-on in its own plane with its backbone, which atom of each
residue makes the contact, and the identity of every interaction. The true
bearings and minimum contact distances are printed when the script runs and
belong in the caption.

What is not: inter-group distances and the 3D orientation of each residue
relative to the drug. A literal projection was tried first and cannot work --
DOR is twisted about its aryl ether, so no viewing direction shows all three
rings open at once, and Val106 lies far enough off the ligand plane that it
projects on top of the triazolinone. That version is in git history (54c3698).

Usage
-----
    PYTHONPATH=. python -m src.analysis.cli.plot_dor_schematic
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

LIG = "2KW"
OFFSET = -3

#: canonical residue -> (partner moiety, label, colour, placement direction,
#:                        stand-off distance)
#: DOR lays out as a butterfly: the pyridinone is the body, and the
#: chlorocyanophenyl and triazolinone rings are the two wings. Each residue is
#: placed off its own wing (or below the body), which keeps all three
#: interactions unobstructed and the figure symmetric. Directions are chosen for
#: legibility, not measured -- see the module docstring.
INTERACTIONS = {
    188: ("chlorocyanophenyl", "π-stacking", "#3B6EA8", (0.62, 0.79), 8.2, +1),
    103: ("triazolinone", "N–H···O=C", "#c0392b", (-0.52, 0.86), 8.2, +1),
    106: ("pyridinone", "hydrophobic", "#B5474B", (0.10, -1.0), 7.6, -1),
}
NAMES = {188: "Tyr188", 106: "Val106", 103: "Lys103"}

MOIETY_COLOR = {"chlorocyanophenyl": "#3B6EA8",
                "pyridinone": "#B5474B",
                "triazolinone": "#3F8C5E"}
#: Explicit label anchors in the RDKit drawing frame. Ring centroids there are
#: chlorocyanophenyl (3.67, 1.59), pyridinone (-0.38, -1.66) and triazolinone
#: (-5.00, 0.26); the interaction lines leave up-right, down and up-left
#: respectively, so labels go on the free side of each ring. Placing these by a
#: rule kept dropping them onto the lines.
MOIETY_LABEL_XY = {"chlorocyanophenyl": (5.85, -2.45),
                   "pyridinone": (-4.95, -4.35),
                   "triazolinone": (-5.60, -1.95)}
HETERO_COLOR = {"N": "#1f4e9c", "O": "#c0392b", "F": "#7d3c98",
                "CL": "#1e8449", "S": "#b7950b"}
RESIDUE_COLOR = "#4a4a4a"


def bonds_within(P, cutoff=1.95):
    return [(i, j) for i, j in itertools.combinations(range(len(P)), 2)
            if np.linalg.norm(P[i] - P[j]) < cutoff]


def ring_moieties(names, coords):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(names)))
    for i, j in bonds_within(coords, 1.85):
        G.add_edge(i, j)
    out = {}
    for ring in [r for r in nx.cycle_basis(G) if len(r) >= 5]:
        nb = {n for a in ring for n in G[a] if n not in ring}
        k = ("chlorocyanophenyl" if any(names[n].upper().startswith("CL") for n in nb)
             else ("triazolinone" if len(ring) == 5 else "pyridinone"))
        for a in ring:
            out[a] = k
    return out


def kabsch2d(A, B):
    """Rotation (with reflection allowed) taking centred A onto centred B."""
    H = A.T @ B
    U, _, Vt = np.linalg.svd(H)
    return U @ Vt


def face_on(P):
    """Project a group onto its own best plane, centred: a readable depiction."""
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c)
    return np.column_stack([(P - c) @ vt[0], (P - c) @ vt[1]])


def rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def draw_group(ax, P2, elements, colour, lw, *, hetero=9.0, dot=170, zorder=3,
               dbl=None, alpha=1.0):
    for i, j in bonds_within(P2 if dbl is None else P2, 1.95 if dbl is None else 1.95):
        ax.plot(*zip(P2[i], P2[j]), color=colour, lw=lw, zorder=zorder,
                solid_capstyle="round", alpha=alpha)
    for i, el in enumerate(elements):
        e = el.upper()
        if e == "C":
            continue
        ax.scatter(*P2[i], s=dot, color="white", zorder=zorder + 1,
                   edgecolors="none", alpha=alpha)
        ax.text(*P2[i], "Cl" if e == "CL" else e, ha="center", va="center",
                fontsize=hetero, fontweight="bold", zorder=zorder + 2,
                color=HETERO_COLOR.get(e, "#333"), alpha=alpha)


def main() -> int:
    import mdtraj as md

    repo = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--structure", type=Path,
                    default=repo / "results/md_runs/wt/rep_01/wt_minimized_rep01.pdb")
    ap.add_argument("--sdf", type=Path, default=repo / "data/ligands/dor.sdf")
    ap.add_argument("--output", type=Path,
                    default=repo / "results/plots/figure1B_dor_schematic.pdf")
    args = ap.parse_args()

    # ---------- clean 2D depiction of DOR (all three rings readable)
    mol = Chem.RemoveHs(Chem.MolFromMolFile(str(args.sdf)))
    Chem.Kekulize(mol, clearAromaticFlags=True)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    D2 = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                   for i in range(mol.GetNumAtoms())])
    d2_names = [a.GetSymbol() for a in mol.GetAtoms()]
    rings2 = {}
    for ring in mol.GetRingInfo().AtomRings():
        nb = {n.GetSymbol() for a in ring
              for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in ring}
        k = ("chlorocyanophenyl" if "Cl" in nb
             else ("triazolinone" if len(ring) == 5 else "pyridinone"))
        for a in ring:
            rings2[a] = k
    cent2 = {m: D2[[i for i, v in rings2.items() if v == m]].mean(axis=0)
             for m in set(rings2.values())}

    # ---------- real geometry: bearings measured in the plane of DOR
    t = md.load(str(args.structure))
    top = t.topology
    X = t.xyz[0] * 10.0
    la = [a for a in top.atoms if a.residue.name == LIG and a.element.symbol != "H"]
    L = X[[a.index for a in la]]
    lc = L.mean(axis=0)
    _, _, vt = np.linalg.svd(L - lc)
    pr = lambda P: np.column_stack([(P - lc) @ vt[0], (P - lc) @ vt[1]])
    moi3 = ring_moieties([a.name for a in la], L)
    lig_xy3 = pr(L)
    cent3 = {m: lig_xy3[[i for i, v in moi3.items() if v == m]].mean(axis=0)
             for m in set(moi3.values())}

    order = ["chlorocyanophenyl", "pyridinone", "triazolinone"]
    A = np.array([cent3[m] for m in order])
    B = np.array([cent2[m] for m in order])
    R = kabsch2d(A - A.mean(axis=0), B - B.mean(axis=0))

    fig, ax = plt.subplots(figsize=(13.0, 9.4))

    for m in set(rings2.values()):
        ax.scatter(*cent2[m], s=4700, color=MOIETY_COLOR[m], alpha=0.14, zorder=0)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        mi, mj = rings2.get(i), rings2.get(j)
        col = MOIETY_COLOR[mi] if (mi and mi == mj) else "#333333"
        p, q = D2[i], D2[j]
        if b.GetBondTypeAsDouble() == 2.0:
            d = q - p
            n = np.array([-d[1], d[0]])
            n = n / (np.linalg.norm(n) + 1e-9) * 0.13
            for sgn in (+1, -1):
                ax.plot(*zip(p + n * sgn, q + n * sgn), color=col, lw=2.8,
                        zorder=3, solid_capstyle="round")
        else:
            ax.plot(*zip(p, q), color=col, lw=2.8, zorder=3, solid_capstyle="round")
    for i, el in enumerate(d2_names):
        if el == "C":
            continue
        ax.scatter(*D2[i], s=230, color="white", zorder=4, edgecolors="none")
        ax.text(*D2[i], el, ha="center", va="center", fontsize=10.5,
                fontweight="bold", zorder=5,
                color=HETERO_COLOR.get(el.upper(), "#333"))
    for m in set(rings2.values()):
        lp = MOIETY_LABEL_XY[m]
        ax.text(lp[0], lp[1], m, ha="center", va="center", fontsize=12.5,
                fontweight="bold", color=MOIETY_COLOR[m], zorder=6)

    # ---------- residues
    print(f"{'residue':9s}{'partner':20s}{'true bearing':>14s}{'true min d':>13s}"
          "   (bearing is measured but NOT used for layout)")
    for can, (part, label, colour, want_dir, standoff, side) in INTERACTIONS.items():
        ai = [a.index for a in top.atoms
              if a.residue.resSeq == can + OFFSET and a.residue.name != LIG
              and a.element.symbol != "H" and a.residue.chain.index == 0]
        P = X[ai]
        elements = [top.atom(k).element.symbol for k in ai]
        m_idx = [i for i, v in moi3.items() if v == part]
        Dm = np.linalg.norm(P[:, None, :] - L[None, m_idx, :], axis=-1)
        ri, lj = np.unravel_index(np.argmin(Dm), Dm.shape)
        mind = float(Dm.min())

        bearing3 = pr(P).mean(axis=0) - cent3[part]  # kept for the printout
        bearing = np.asarray(want_dir, dtype=float)
        bearing /= np.linalg.norm(bearing)

        # readable depiction, rotated so the contacting atom faces the partner
        F = face_on(P)
        v = F[ri] - F.mean(axis=0)
        if np.linalg.norm(v) < 1e-6:
            v = np.array([1.0, 0.0])
        cur = np.arctan2(v[1], v[0])
        want = np.arctan2(-bearing[1], -bearing[0])
        F = F @ rot(want - cur).T
        pos = cent2[part] + bearing * standoff
        F = F + pos

        draw_group(ax, F, elements, RESIDUE_COLOR, 2.1, hetero=8.6, dot=150,
                   zorder=2)

        # interaction line: contacting atom to the partner ring, clear of both
        a_pt, b_pt = F[ri], cent2[part]
        u = (a_pt - b_pt) / np.linalg.norm(a_pt - b_pt)
        p0 = b_pt + u * 1.85
        p1 = a_pt - u * 1.05
        ax.plot(*zip(p0, p1), color=colour, lw=2.4, linestyle=(0, (2.5, 2.2)),
                zorder=5, solid_capstyle="round")
        perp = np.array([-u[1], u[0]])
        mid = p0 + (p1 - p0) * 0.58 + perp * 1.30 * side
        ang = np.degrees(np.arctan2(u[1], u[0]))
        ang = ang - 180 if ang > 90 else (ang + 180 if ang < -90 else ang)
        ax.text(mid[0], mid[1], label, ha="center", va="center", fontsize=11,
                fontweight="bold", color=colour, zorder=8, rotation=ang,
                rotation_mode="anchor",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="none", alpha=0.95))
        lab = F.mean(axis=0) + bearing * (float(np.max((F - F.mean(axis=0)) @ bearing)) + 1.5)
        ax.text(lab[0], lab[1], NAMES[can], ha="center", va="center",
                fontsize=13, fontweight="bold", color=RESIDUE_COLOR, zorder=7)
        print(f"{NAMES[can]:9s}{part:20s}"
              f"{np.degrees(np.arctan2(bearing3[1], bearing3[0])):11.0f}°"
              f"{mind:11.2f} Å")

    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.14)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
