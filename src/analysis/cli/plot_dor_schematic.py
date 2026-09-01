#!/usr/bin/env python3
"""Figure 1B -- doravirine and its NNIBP environment, drawn from real geometry.

Why this exists
---------------
The manuscript names DOR's moieties throughout (chlorocyanophenyl, pyridinone,
triazolinone) and the residues that engage them, but never shows which part of
the molecule is which, or where each residue actually sits.

Rather than a hand-placed cartoon, this projects the real wild-type structure
onto the best plane through DOR's heavy atoms (principal axes of the ligand), so
**every position and orientation in the figure is the true one** -- residues
appear where they are relative to the drug, at the angle they actually adopt.
The only liberty taken is a uniform radial expansion (--expand) to stop
contacting groups from overlapping on the page; this preserves every direction
exactly and is disclosed in the caption.

Out-of-plane depth is shown by opacity: groups lying near DOR's plane are drawn
solid, those further above or below are faded. Depth in angstroms is available
in the printed table.

Bonds are inferred from interatomic distance, so no residue templates are
needed and the backbone (N-CA-C=O) is drawn along with the side chain.

Usage
-----
    PYTHONPATH=. python -m src.analysis.cli.plot_dor_schematic
    PYTHONPATH=. python -m src.analysis.cli.plot_dor_schematic --expand 1.6
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LIG = "2KW"
OFFSET = -3  # topology resSeq = canonical + OFFSET

#: canonical residue -> (display name, label direction hint in projected plane)
RESIDUES = {
    188: "Tyr188",
    229: "Trp229",
    227: "Phe227",
    181: "Tyr181",
    103: "Lys103",
    106: "Val106",
    318: "Tyr318",
    190: "Gly190",
    # Ser105 is deliberately omitted: it makes no contact with DOR in the
    # wild-type pose (4.36 A minimum heavy-atom distance, 0.0 contacts at
    # 4.0 A). It matters only to the V106A displacement discussed later, and it
    # projects into the most crowded region of this figure. Add it back by
    # restoring the entry if that discussion needs it here.
}

MOIETY_COLOR = {
    "chlorocyanophenyl": "#3B6EA8",
    "pyridinone": "#B5474B",
    "triazolinone": "#3F8C5E",
}
MOIETY_LABEL = {
    "chlorocyanophenyl": "chlorocyanophenyl",
    "pyridinone": "pyridinone",
    "triazolinone": "triazolinone",
}
HETERO_COLOR = {"N": "#1f4e9c", "O": "#c0392b", "F": "#7d3c98",
                "CL": "#1e8449", "S": "#b7950b"}
RESIDUE_COLOR = "#4a4a4a"


def bonds_within(coords: np.ndarray, cutoff: float = 1.95):
    out = []
    for i, j in itertools.combinations(range(len(coords)), 2):
        if np.linalg.norm(coords[i] - coords[j]) < cutoff:
            out.append((i, j))
    return out


def moiety_of_ligand(names, coords) -> dict[int, str]:
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(len(names)))
    for i, j in bonds_within(coords, 1.85):
        G.add_edge(i, j)
    out: dict[int, str] = {}
    for ring in [c for c in nx.cycle_basis(G) if len(c) >= 5]:
        nb = {n for a in ring for n in G[a] if n not in ring}
        nm = ("chlorocyanophenyl" if any(names[n].upper().startswith("CL") for n in nb)
              else ("triazolinone" if len(ring) == 5 else "pyridinone"))
        for a in ring:
            out[a] = nm
    return out


def main() -> int:
    import mdtraj as md

    repo = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--structure", type=Path,
                    default=repo / "results/md_runs/wt/rep_01/wt_minimized_rep01.pdb")
    ap.add_argument("--expand", type=float, default=1.72,
                    help="uniform radial expansion of residues about the ligand "
                         "centroid; 1.0 keeps true distances but overlaps")
    ap.add_argument("--output", type=Path,
                    default=repo / "results/plots/figure1B_dor_schematic.pdf")
    args = ap.parse_args()

    t = md.load(str(args.structure))
    top = t.topology
    X = t.xyz[0] * 10.0

    lig_atoms = [a for a in top.atoms
                 if a.residue.name == LIG and a.element.symbol != "H"]
    lig_idx = [a.index for a in lig_atoms]
    lig_names = [a.name for a in lig_atoms]
    L = X[lig_idx]
    centre = L.mean(axis=0)

    # principal plane of the ligand
    _, sv, vt = np.linalg.svd(L - centre)
    e1, e2, e3 = vt[0], vt[1], vt[2]
    proj = lambda P: np.column_stack([(P - centre) @ e1, (P - centre) @ e2])
    depth = lambda P: (P - centre) @ e3

    lig_xy = proj(L)
    moiety = moiety_of_ligand(lig_names, L)

    fig, ax = plt.subplots(figsize=(13.0, 10.0))

    # ---- ligand
    for m in set(moiety.values()):
        pts = lig_xy[[i for i, v in moiety.items() if v == m]]
        ax.scatter(*pts.mean(axis=0), s=5200, color=MOIETY_COLOR[m],
                   alpha=0.13, zorder=0)
    for i, j in bonds_within(L, 1.85):
        mi, mj = moiety.get(i), moiety.get(j)
        col = MOIETY_COLOR[mi] if (mi and mi == mj) else "#333333"
        ax.plot(*zip(lig_xy[i], lig_xy[j]), color=col, lw=3.4, zorder=3,
                solid_capstyle="round")
    for i, nm in enumerate(lig_names):
        el = "".join(ch for ch in nm if ch.isalpha()).upper()
        el = "CL" if el.startswith("CL") else el[:1]
        if el == "C":
            continue
        ax.scatter(*lig_xy[i], s=250, color="white", zorder=4, edgecolors="none")
        ax.text(*lig_xy[i], "Cl" if el == "CL" else el, ha="center", va="center",
                fontsize=10.5, fontweight="bold", zorder=5,
                color=HETERO_COLOR.get(el, "#333"))

    # Moiety labels are pushed clear of the ligand along a fixed direction per
    # moiety: purely radial placement put them back on top of the rings, since
    # the three ring centroids are nearly collinear in this projection.
    LABEL_DIR = {"chlorocyanophenyl": np.array([-0.95, -0.55]),
                 "pyridinone": np.array([0.62, 1.0]),
                 "triazolinone": np.array([0.30, -1.0])}
    for m in set(moiety.values()):
        pts = lig_xy[[i for i, v in moiety.items() if v == m]]
        c = pts.mean(axis=0)
        d = LABEL_DIR[m] / np.linalg.norm(LABEL_DIR[m])
        reach = float(np.max((pts - c) @ d))
        lp = c + d * (reach + 1.8)
        ax.text(lp[0], lp[1], MOIETY_LABEL[m], ha="center", va="center",
                fontsize=13, fontweight="bold", color=MOIETY_COLOR[m], zorder=6)

    # ---- residues, at their true projected positions and orientations
    print(f"{'residue':9s}{'u':>7s}{'v':>7s}{'depth':>8s}{'min d':>8s}")
    rows = []
    shifted: dict[int, dict[int, np.ndarray]] = {}
    for can, name in RESIDUES.items():
        ai = [a.index for a in top.atoms
              if a.residue.resSeq == can + OFFSET
              and a.residue.name != LIG
              and a.element.symbol != "H"
              and a.residue.chain.index == 0]
        if not ai:
            continue
        P = X[ai]
        xy = proj(P)
        w = depth(P).mean()
        mind = float(np.linalg.norm(P[:, None, :] - L[None, :, :], axis=-1).min())
        # uniform radial expansion about the ligand centroid: preserves every
        # direction and the shape of each residue, only adds page separation
        rc = xy.mean(axis=0)
        shift = rc * (args.expand - 1.0)
        xy = xy + shift
        alpha = float(np.clip(1.0 - abs(w) / 11.0, 0.42, 1.0))
        for i, j in bonds_within(P, 1.95):
            ax.plot(*zip(xy[i], xy[j]), color=RESIDUE_COLOR, lw=2.0,
                    alpha=alpha, zorder=2, solid_capstyle="round")
        for k, a_i in enumerate(ai):
            el = top.atom(a_i).element.symbol.upper()
            if el == "C":
                continue
            ax.scatter(*xy[k], s=150, color="white", zorder=2.5,
                       edgecolors="none", alpha=alpha)
            ax.text(*xy[k], el.capitalize() if el == "CL" else el, ha="center",
                    va="center", fontsize=8.4, fontweight="bold", zorder=2.6,
                    color=HETERO_COLOR.get(el, "#333"), alpha=alpha)
        # Place the label just beyond the residue's own outer edge, along the
        # direction pointing away from the ligand. Using a fixed radius from the
        # origin instead put labels on top of neighbouring residues.
        lc = xy.mean(axis=0)
        d = lc / (np.linalg.norm(lc) + 1e-9)
        reach = float(np.max((xy - lc) @ d))
        lab = lc + d * (reach + 1.7)
        ha = "left" if d[0] > 0.35 else ("right" if d[0] < -0.35 else "center")
        va = "bottom" if d[1] > 0.35 else ("top" if d[1] < -0.35 else "center")
        ax.text(lab[0], lab[1], name, ha=ha, va=va, fontsize=12.5,
                fontweight="bold", color=RESIDUE_COLOR, alpha=max(alpha, 0.8),
                zorder=6)
        shifted[can] = dict(zip(ai, xy))
        rows.append((name, rc[0], rc[1], w, mind))
        print(f"{name:9s}{rc[0]:7.1f}{rc[1]:7.1f}{w:8.1f}{mind:8.2f}")

    # ---- the one interaction drawn explicitly: Lys103 main chain -> triazolinone
    tri = [i for i, v in moiety.items() if v == "triazolinone"]
    donor = None
    for i in tri:
        if lig_names[i].upper().startswith("N"):
            nb = [j for j in range(len(L))
                  if j != i and np.linalg.norm(L[i] - L[j]) < 1.85]
            if len(nb) == 2:  # the ring N bearing H, not the substituted one
                donor = i
    lys_o = [a.index for a in top.atoms
             if a.residue.resSeq == 103 + OFFSET and a.name == "O"
             and a.residue.chain.index == 0]
    if donor is not None and lys_o and 103 in shifted:
        p = lig_xy[donor]
        q = shifted[103].get(lys_o[0])
        if q is not None:
            v = q - p
            n = np.linalg.norm(v)
            u = v / n
            ax.plot(*zip(p + u * 0.42, q - u * 0.42), color="#c0392b", lw=2.2,
                    linestyle=(0, (2.6, 2.2)), zorder=5)
            ax.scatter(*p, s=300, facecolor="none", edgecolor="#c0392b", lw=1.9,
                       zorder=6)
            mid = (p + q) / 2
            perp = np.array([-u[1], u[0]])
            lp = mid + perp * 0.85
            ax.text(lp[0], lp[1], "N–H···O=C", ha="center", va="center",
                    fontsize=10.5, fontweight="bold", color="#c0392b", zorder=7,
                    rotation=np.degrees(np.arctan2(u[1], u[0])),
                    rotation_mode="anchor",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="none", alpha=0.9))

    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.13)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
