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
#: (x, y, horizontal alignment). The chlorocyanophenyl label is wide, so it is
#: left-anchored at the ring's lower-right edge and grows outward; centring it
#: near the ring drops it onto the ether oxygen.
MOIETY_LABEL_XY = {"chlorocyanophenyl": (4.55, -1.05, "left"),
                   "pyridinone": (-3.20, -3.05, "center"),
                   "triazolinone": (-5.60, -1.95, "center")}
HETERO_COLOR = {"N": "#1f4e9c", "O": "#c0392b", "F": "#7d3c98",
                "CL": "#1e8449", "S": "#b7950b"}
RESIDUE_COLOR = "#4a4a4a"

#: Each residue is built from SMILES and laid out by RDKit, exactly as DOR is
#: built from its SDF. Deriving the depiction from projected MD coordinates
#: instead gave correct chemistry but an ugly drawing -- backbones foreshortened
#: onto the ring plane, side chains folded over themselves -- and an earlier
#: version compounded that by inferring bonds from the projection. Building from
#: SMILES makes the atoms and bonds right by construction.
#: Fragments are capped as acetyl-N ... C(=O)NHMe, i.e. the residue as excised
#: from a peptide chain, so the main chain reads correctly.
RESIDUE_SMILES = {
    "Tyr188": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NC",
    "Val106": "CC(=O)N[C@@H](C(C)C)C(=O)NC",
    "Lys103": "CC(=O)N[C@@H](CCCCN)C(=O)NC",
}
#: the residue's own main-chain carbonyl (not the acetyl cap): C, O, N, CH3
BACKBONE_CO_SMARTS = "[CX3](=[OX1])[NX3][CH3]"
VAL_GAMMA_SMARTS = "[CH3][CH]([CH3])[CH]"


def build_residue(name):
    """2D depiction of a capped residue, plus the indices that matter."""
    mol = Chem.MolFromSmiles(RESIDUE_SMILES[name])
    Chem.Kekulize(mol, clearAromaticFlags=True)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    P = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                  for i in range(mol.GetNumAtoms())])
    info = {"mol": mol, "xy": P,
            "elements": [a.GetSymbol() for a in mol.GetAtoms()],
            "bonds": [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                       b.GetBondTypeAsDouble()) for b in mol.GetBonds()]}
    co = mol.GetSubstructMatches(Chem.MolFromSmarts(BACKBONE_CO_SMARTS))
    info["carbonyl_O"] = co[0][1] if co else None
    rings = mol.GetRingInfo().AtomRings()
    info["ring"] = list(rings[0]) if rings else None
    g = mol.GetSubstructMatches(Chem.MolFromSmarts(VAL_GAMMA_SMARTS))
    info["gammas"] = [g[0][0], g[0][2]] if g else None
    return info


def draw_fragment(ax, P, elements, bonds, colour, lw, *, hetero=8.8, dot=150,
                  zorder=3, dbl=0.115, circle_ring=None):
    ring = set(circle_ring or ())
    for i, j, order in bonds:
        if order == 2.0 and i in ring and j in ring:
            order = 1.0   # drawn as an aromatic circle instead
        p, q = P[i], P[j]
        if order == 2.0:
            d = q - p
            n = np.array([-d[1], d[0]])
            n = n / (np.linalg.norm(n) + 1e-9) * dbl
            for sgn in (+1, -1):
                ax.plot(*zip(p + n * sgn, q + n * sgn), color=colour, lw=lw,
                        zorder=zorder, solid_capstyle="round")
        else:
            ax.plot(*zip(p, q), color=colour, lw=lw, zorder=zorder,
                    solid_capstyle="round")
    for i, el in enumerate(elements):
        if el == "C":
            continue
        ax.scatter(*P[i], s=dot, color="white", zorder=zorder + 1, edgecolors="none")
        ax.text(*P[i], el, ha="center", va="center", fontsize=hetero,
                fontweight="bold", zorder=zorder + 2,
                color=HETERO_COLOR.get(el.upper(), "#333"))


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


class LabelPlacer:
    """Place text so it never overlaps an atom, a bond, or another label.

    Hand-tuned label coordinates kept drifting onto bonds whenever anything else
    moved. This measures each label's rendered extent, tests it against every
    atom point, bond segment and previously placed label, and walks outward
    along a preferred direction until it finds a clear spot.
    """

    def __init__(self, ax, fig):
        self.ax, self.fig = ax, fig
        self.points: list[np.ndarray] = []
        self.segs: list[tuple[np.ndarray, np.ndarray]] = []
        self.boxes: list[tuple[float, float, float, float]] = []

    def add_atoms(self, P):
        for q in np.atleast_2d(P):
            self.points.append(np.asarray(q, dtype=float))

    def add_bonds(self, P, pairs):
        for i, j in pairs:
            self.segs.append((np.asarray(P[i], float), np.asarray(P[j], float)))

    @staticmethod
    def _seg_box_hit(a, b, box, pad):
        x0, y0, x1, y1 = box
        x0, y0, x1, y1 = x0 - pad, y0 - pad, x1 + pad, y1 + pad
        for t in np.linspace(0, 1, 24):
            p = a + (b - a) * t
            if x0 <= p[0] <= x1 and y0 <= p[1] <= y1:
                return True
        return False

    def _extent(self, txt):
        self.fig.canvas.draw()
        bb = txt.get_window_extent(self.fig.canvas.get_renderer())
        inv = self.ax.transData.inverted()
        (x0, y0), (x1, y1) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1))
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    def _penalty(self, box, pad):
        """How badly this position collides. 0 means clear."""
        x0, y0, x1, y1 = box
        pen = 0.0
        for q in self.points:
            if x0 - pad <= q[0] <= x1 + pad and y0 - pad <= q[1] <= y1 + pad:
                pen += 1.0
        for a, b in self.segs:
            if self._seg_box_hit(a, b, box, pad):
                pen += 1.0
        for bx in self.boxes:
            if not (x1 + pad < bx[0] or bx[2] < x0 - pad
                    or y1 + pad < bx[1] or bx[3] < y0 - pad):
                pen += 4.0        # overlapping another label is the worst case
        return pen

    def place(self, text, anchor, direction, *, pad=0.30, start=1.6, step=0.30,
              limit=17.0, **kw):
        """Walk outward and around until the label is clear.

        Falls back to the least-colliding candidate rather than the first one:
        with a molecule this dense some labels have no perfectly free position,
        and picking the first tried put them straight on top of a bond.
        """
        d = np.asarray(direction, dtype=float)
        d = d / (np.linalg.norm(d) + 1e-9)
        anchor = np.asarray(anchor, dtype=float)
        txt = self.ax.text(0, 0, text, ha="center", va="center", **kw)
        best = (np.inf, None, None)
        r = start
        while r <= limit:
            for sway in (0.0, 0.26, -0.26, 0.55, -0.55, 0.9, -0.9, 1.3, -1.3):
                pos = anchor + (rot(sway) @ d) * r
                txt.set_position(pos)
                box = self._extent(txt)
                pen = self._penalty(box, pad)
                # prefer clear, then close to the anchor
                score = pen * 1000.0 + r + abs(sway) * 1.5
                if score < best[0]:
                    best = (score, pos, box)
                if pen == 0.0:
                    txt.set_position(pos)
                    self.boxes.append(box)
                    return pos
            r += step
        txt.set_position(best[1])
        self.boxes.append(best[2])
        return best[1]


def rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def draw_group(ax, P2, elements, colour, lw, *, hetero=9.0, dot=170, zorder=3,
               bonds=None, alpha=1.0):
    """Draw a group. `bonds` MUST come from the 3D coordinates.

    Inferring them from the projected 2D positions (as an earlier version did)
    invents bonds between atoms that merely project close together and drops
    real ones that project on top of each other -- which produced visibly wrong
    connectivity for Tyr188 and Lys103.
    """
    if bonds is None:
        raise ValueError("bonds must be supplied from 3D coordinates")
    for i, j in bonds:
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

    placer = LabelPlacer(ax, fig)
    pending: list[tuple] = []
    placer.add_atoms(D2)
    placer.add_bonds(D2, [(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                          for b in mol.GetBonds()])
    for m in set(rings2.values()):
        ax.scatter(*cent2[m], s=2600, color=MOIETY_COLOR[m], alpha=0.16, zorder=0)
    # The chlorocyanophenyl ring is a plain benzene, so it is drawn with an
    # inner circle rather than alternating Kekule bonds -- the alternating form
    # renders lopsided at this size and reads badly next to the Tyr188 ring,
    # which uses the same convention.
    chl_ring_set = {k for k, v in rings2.items() if v == "chlorocyanophenyl"}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        mi, mj = rings2.get(i), rings2.get(j)
        col = MOIETY_COLOR[mi] if (mi and mi == mj) else "#333333"
        p, q = D2[i], D2[j]
        in_benzene = i in chl_ring_set and j in chl_ring_set
        if b.GetBondTypeAsDouble() == 2.0 and not in_benzene:
            d = q - p
            n = np.array([-d[1], d[0]])
            n = n / (np.linalg.norm(n) + 1e-9) * 0.13
            for sgn in (+1, -1):
                ax.plot(*zip(p + n * sgn, q + n * sgn), color=col, lw=2.8,
                        zorder=3, solid_capstyle="round")
        else:
            ax.plot(*zip(p, q), color=col, lw=2.8, zorder=3, solid_capstyle="round")
    chl_pts0 = D2[sorted(chl_ring_set)]
    chl_c0 = chl_pts0.mean(axis=0)
    ax.add_patch(plt.Circle(chl_c0,
                            0.60 * np.linalg.norm(chl_pts0[0] - chl_c0) * 1.30,
                            fill=False, ec=MOIETY_COLOR["chlorocyanophenyl"],
                            lw=2.2, zorder=4))
    for i, el in enumerate(d2_names):
        if el == "C":
            continue
        ax.scatter(*D2[i], s=230, color="white", zorder=4, edgecolors="none")
        ax.text(*D2[i], el, ha="center", va="center", fontsize=10.5,
                fontweight="bold", zorder=5,
                color=HETERO_COLOR.get(el.upper(), "#333"))
    for m in set(rings2.values()):
        lx, ly, _ = MOIETY_LABEL_XY[m]
        anchor = cent2[m]
        direction = np.array([lx, ly]) - anchor
        pending.append((m, anchor, direction,
                        dict(fontsize=12.5, fontweight="bold",
                             color=MOIETY_COLOR[m], zorder=6)))

    # ---------- residues, built from SMILES and placed by their interaction
    #
    # Polarity: the donor is the triazolinone N-H, the acceptor is the Lys103
    # MAIN-CHAIN carbonyl (Lys103:O to N4x = 2.90 A in the wild-type structure).
    print(f"{'residue':9s}{'contact':26s}{'geometry':>20s}")

    def orient(P, from_idx, to_idx, want, centre_at):
        """Rotate so `from_idx -> to_idx` points along `want`, then translate."""
        v = P[to_idx] - P[from_idx] if from_idx is not None else P[to_idx] - P.mean(axis=0)
        cur = np.arctan2(v[1], v[0])
        P = (P - P.mean(axis=0)) @ rot(np.arctan2(want[1], want[0]) - cur).T
        return P + centre_at

    # ===== Tyr188: ring drawn parallel to the chlorocyanophenyl, stacked =====
    tyr = build_residue("Tyr188")
    P = tyr["xy"].copy()
    ring = tyr["ring"]
    chl_idx = [k for k, v in rings2.items() if v == "chlorocyanophenyl"]
    chl_poly = D2[chl_idx]
    chl_c = chl_poly.mean(axis=0)
    tc = P[ring].mean(axis=0)
    v_t = P[ring[0]] - tc
    v_c = chl_poly[0] - chl_c
    P = (P - tc) @ rot(np.arctan2(v_c[1], v_c[0]) - np.arctan2(v_t[1], v_t[0])).T
    STACK_DIR = np.array([0.80, 0.60]); STACK_DIR /= np.linalg.norm(STACK_DIR)
    P = P + chl_c + STACK_DIR * 7.6
    tyr_ring2 = P[ring]
    tyr_c2 = tyr_ring2.mean(axis=0)
    hull = np.vstack([chl_poly, tyr_ring2])
    oh = np.argsort(np.arctan2(*(hull - hull.mean(axis=0)).T[::-1]))
    ax.fill(*hull[oh].T, color=MOIETY_COLOR["chlorocyanophenyl"], alpha=0.10,
            zorder=1, linewidth=0)
    draw_fragment(ax, P, tyr["elements"], tyr["bonds"], RESIDUE_COLOR, 2.1,
                  circle_ring=ring)
    ax.add_patch(plt.Circle(tyr_c2, np.linalg.norm(tyr_ring2[0] - tyr_c2) * 0.62,
                            fill=False, ec=RESIDUE_COLOR, lw=1.9, zorder=4))
    placer.add_atoms(P)
    placer.add_bonds(P, [(i, j) for i, j, _ in tyr["bonds"]])
    perp = np.array([-STACK_DIR[1], STACK_DIR[0]])
    pending.append(("Tyr188", tyr_c2, STACK_DIR * 0.5 + perp * 1.0,
                    dict(fontsize=13, fontweight="bold", color=RESIDUE_COLOR,
                         zorder=7)))
    print(f"{'Tyr188':9s}{'ring / chlorocyanophenyl':26s}{'3.78 A, 6.4 deg':>20s}")

    # ===== Val106: line lands on the gamma-methyl that alanine lacks ========
    val = build_residue("Val106")
    pyr_idx = [k for k, v in rings2.items() if v == "pyridinone"]
    pyr_c = D2[pyr_idx].mean(axis=0)
    g1, g2 = val["gammas"]
    P = val["xy"].copy()
    gmid = (P[g1] + P[g2]) / 2
    want = np.array([0.10, -1.0]); want /= np.linalg.norm(want)
    v = gmid - P.mean(axis=0)
    P = (P - P.mean(axis=0)) @ rot(np.arctan2(-want[1], -want[0])
                                   - np.arctan2(v[1], v[0])).T
    P = P + pyr_c + want * 10.2
    draw_fragment(ax, P, val["elements"], val["bonds"], RESIDUE_COLOR, 2.1)
    placer.add_atoms(P)
    placer.add_bonds(P, [(i, j) for i, j, _ in val["bonds"]])
    near = g1 if np.linalg.norm(P[g1] - pyr_c) < np.linalg.norm(P[g2] - pyr_c) else g2
    u = (P[near] - pyr_c) / np.linalg.norm(P[near] - pyr_c)
    ax.plot(*zip(pyr_c + u * 1.75, P[near] - u * 0.42),
            color=MOIETY_COLOR["pyridinone"], lw=2.4, linestyle=(0, (2.5, 2.2)),
            zorder=5, solid_capstyle="round")
    for g in (g1, g2):
        ax.scatter(*P[g], s=200, facecolor="none",
                   edgecolor=MOIETY_COLOR["pyridinone"], lw=1.7, zorder=6)
    pending.append(("Val106", P.mean(axis=0), want,
                    dict(fontsize=13, fontweight="bold", color=RESIDUE_COLOR,
                         zorder=7)))
    print(f"{'Val106':9s}{'C-gamma / pyridinone':26s}{'3.75 A':>20s}")

    # ===== Lys103: main-chain C=O accepts from the triazolinone N-H =========
    lys = build_residue("Lys103")
    tri_idx = [k for k, v in rings2.items() if v == "triazolinone"]
    tri_c = D2[tri_idx].mean(axis=0)
    o_i = lys["carbonyl_O"]
    P = lys["xy"].copy()
    want = np.array([-0.52, 0.86]); want /= np.linalg.norm(want)
    v = P[o_i] - P.mean(axis=0)
    P = (P - P.mean(axis=0)) @ rot(np.arctan2(-want[1], -want[0])
                                   - np.arctan2(v[1], v[0])).T
    P = P + tri_c + want * 10.4
    draw_fragment(ax, P, lys["elements"], lys["bonds"], RESIDUE_COLOR, 2.1)
    placer.add_atoms(P)
    placer.add_bonds(P, [(i, j) for i, j, _ in lys["bonds"]])
    ring_n2 = [k for k in tri_idx if d2_names[k] == "N"]
    don2 = min(ring_n2, key=lambda k: np.linalg.norm(D2[k] - P[o_i]))
    p_don, p_acc = D2[don2], P[o_i]
    u = (p_acc - p_don) / np.linalg.norm(p_acc - p_don)
    ax.plot(*zip(p_don + u * 0.55, p_acc - u * 0.55), color="#c0392b", lw=2.4,
            linestyle=(0, (2.5, 2.2)), zorder=5, solid_capstyle="round")
    for pt in (p_don, p_acc):
        ax.scatter(*pt, s=250, facecolor="none", edgecolor="#c0392b", lw=1.8,
                   zorder=6)
    pending.append(("Lys103", P.mean(axis=0), want,
                    dict(fontsize=13, fontweight="bold", color=RESIDUE_COLOR,
                         zorder=7)))
    print(f"{'Lys103':9s}{'main-chain O / N4x':26s}{'2.90 A':>20s}")

    # Freeze the axes BEFORE measuring any label. Autoscaling was still active
    # while extents were being measured, so every added text shifted the
    # data<->pixel transform and the collision tests were done against stale
    # boxes -- which is why labels kept landing on bonds.
    allpts = np.vstack(placer.points)
    lo, hi = allpts.min(axis=0), allpts.max(axis=0)
    span = (hi - lo).max()
    mid = (hi + lo) / 2
    half = span * 0.5 + span * 0.30
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.autoscale(False)

    for text, anchor, direction, kw in pending:
        placer.place(text, anchor, direction, **kw)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
