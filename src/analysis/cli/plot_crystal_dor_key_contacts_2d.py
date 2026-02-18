#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
COL = {
    "polar": ("#bfe3ff", "#1f77b4", (0.12, 0.47, 0.71)),
    "hydrophobic": ("#ffd9bf", "#d95f02", (0.85, 0.37, 0.01)),
    "mixed": ("#e7ccff", "#7a2c91", (0.53, 0.12, 0.47)),
}


def _elem(sym: str) -> str:
    s = str(sym).strip().upper()
    if s in {"", "D", "T"}:
        return "H" if s in {"D", "T"} else "C"
    if s == "CL":
        return "Cl"
    if s == "BR":
        return "Br"
    return s if len(s) == 1 else s[0] + s[1:].lower()


def _ligand_from_cif(cif: Path) -> tuple[Chem.Mol, dict[str, int]]:
    d = MMCIF2Dict(str(cif))
    rw = Chem.RWMol()
    idx = {}
    for c, a, e in zip(d["_chem_comp_atom.comp_id"], d["_chem_comp_atom.atom_id"], d["_chem_comp_atom.type_symbol"]):
        if c == "2KW":
            idx[str(a).strip()] = rw.AddAtom(Chem.Atom(_elem(e)))
    bmap = {"sing": Chem.BondType.SINGLE, "doub": Chem.BondType.DOUBLE, "trip": Chem.BondType.TRIPLE, "arom": Chem.BondType.AROMATIC}
    arom = set()
    for c, a1, a2, bo in zip(d["_chem_comp_bond.comp_id"], d["_chem_comp_bond.atom_id_1"], d["_chem_comp_bond.atom_id_2"], d["_chem_comp_bond.value_order"]):
        if c != "2KW":
            continue
        i, j = idx.get(str(a1).strip()), idx.get(str(a2).strip())
        if i is None or j is None:
            continue
        bt = bmap.get(str(bo).lower().strip(), Chem.BondType.SINGLE)
        rw.AddBond(i, j, bt)
        if bt == Chem.BondType.AROMATIC:
            arom.update([i, j])
    m = rw.GetMol()
    for i in arom:
        m.GetAtomWithIdx(i).SetIsAromatic(True)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.AROMATIC:
            b.SetIsAromatic(True)
    Chem.SanitizeMol(m)
    rdDepictor.Compute2DCoords(m)
    return m, idx


def _cat(cats: set[str]) -> str:
    return "mixed" if {"polar", "hydrophobic"}.issubset(cats) else ("polar" if "polar" in cats else "hydrophobic")


def _seg_point_dist(a: np.ndarray, b: np.ndarray, p: np.ndarray) -> float:
    ab = b - a
    den = float(np.dot(ab, ab))
    if den < 1e-12:
        return float(np.hypot(*(p - a)))
    t = float(np.clip(np.dot(p - a, ab) / den, 0.0, 1.0))
    q = a + t * ab
    return float(np.hypot(*(p - q)))


def _orient(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _on_seg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> bool:
    return (
        min(a[0], b[0]) - 1e-9 <= c[0] <= max(a[0], b[0]) + 1e-9
        and min(a[1], b[1]) - 1e-9 <= c[1] <= max(a[1], b[1]) + 1e-9
    )


def _seg_intersects(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
    o1 = _orient(a, b, c)
    o2 = _orient(a, b, d)
    o3 = _orient(c, d, a)
    o4 = _orient(c, d, b)
    if (o1 > 0.0) != (o2 > 0.0) and (o3 > 0.0) != (o4 > 0.0):
        return True
    if abs(o1) < 1e-9 and _on_seg(a, b, c):
        return True
    if abs(o2) < 1e-9 and _on_seg(a, b, d):
        return True
    if abs(o3) < 1e-9 and _on_seg(c, d, a):
        return True
    if abs(o4) < 1e-9 and _on_seg(c, d, b):
        return True
    return False


def _seg_seg_dist(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    if _seg_intersects(a, b, c, d):
        return 0.0
    return min(
        _seg_point_dist(a, b, c),
        _seg_point_dist(a, b, d),
        _seg_point_dist(c, d, a),
        _seg_point_dist(c, d, b),
    )


def _path_has_hard_clash(
    path: list[np.ndarray],
    li: int,
    atom_geom: list[tuple[int, np.ndarray]],
    bond_geom: list[tuple[np.ndarray, np.ndarray, int, int]],
    used_segs: list[tuple[np.ndarray, np.ndarray]],
    atom_hard: float,
    bond_hard: float,
    conn_hard: float,
) -> bool:
    segs = list(zip(path[:-1], path[1:]))
    for p, q in segs:
        for ai, ap in atom_geom:
            if ai == li:
                continue
            if _seg_point_dist(p, q, ap) < atom_hard:
                return True
        for bp, bq, bi, bj in bond_geom:
            bthr = bond_hard * 0.58 if li in (bi, bj) else bond_hard
            if _seg_seg_dist(p, q, bp, bq) < bthr:
                return True
        for up, uq in used_segs:
            if _seg_seg_dist(p, q, up, uq) < conn_hard:
                return True
    return False


def _path_soft_score(
    path: list[np.ndarray],
    li: int,
    atom_geom: list[tuple[int, np.ndarray]],
    bond_geom: list[tuple[np.ndarray, np.ndarray, int, int]],
    used_segs: list[tuple[np.ndarray, np.ndarray]],
) -> float:
    segs = list(zip(path[:-1], path[1:]))
    plen = float(sum(np.hypot(*(q - p)) for p, q in segs))
    score = plen * 0.013 + (len(path) - 2) * 7.5
    for p, q in segs:
        for ai, ap in atom_geom:
            if ai == li:
                continue
            da = _seg_point_dist(p, q, ap)
            if da < 12.0:
                score += (12.0 - da) * 10.0
        for bp, bq, bi, bj in bond_geom:
            db = _seg_seg_dist(p, q, bp, bq)
            bsoft = 5.2 if li in (bi, bj) else 7.4
            if db < bsoft:
                score += (bsoft - db) * (10.0 if li in (bi, bj) else 16.0)
        for up, uq in used_segs:
            du = _seg_seg_dist(p, q, up, uq)
            if du < 5.8:
                score += (5.8 - du) * 12.0
    return score


def _refine_node_angles_for_connector_clearance(
    ent: list[dict[str, object]],
    ring_cx: float,
    ring_cy: float,
    ring_rx: float,
    ring_ry: float,
    br: float,
    xy: dict[int, object],
    all_xy: dict[int, object],
    ml: float,
    mt: float,
    mol: Chem.Mol,
) -> None:
    n = len(ent)
    if n == 0:
        return

    reff = max(1.0, min(ring_rx, ring_ry))
    dmin = min((2.0 * np.pi / n) * 0.9, (2.0 * br + 18.0) / reff)
    atom_geom = [(i, np.array([float(ml + all_xy[i].x), float(mt + all_xy[i].y)], float)) for i in range(mol.GetNumAtoms())]
    bond_geom = []
    for b in mol.GetBonds():
        bi, bj = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        bp = np.array([float(ml + all_xy[bi].x), float(mt + all_xy[bi].y)], float)
        bq = np.array([float(ml + all_xy[bj].x), float(mt + all_xy[bj].y)], float)
        bond_geom.append((bp, bq, bi, bj))

    theta = np.array(
        [
            np.arctan2(
                (float(e["sy"]) - ring_cy) / max(ring_ry, 1e-6),
                (float(e["sx"]) - ring_cx) / max(ring_rx, 1e-6),
            )
            for e in ent
        ],
        float,
    )
    base = theta.copy()
    dts = [0.0, -0.03, 0.03, -0.06, 0.06, -0.1, 0.1, -0.14, 0.14, -0.18, 0.18]

    def _score(i: int, ang: float) -> float:
        node = np.array([ring_cx + ring_rx * np.cos(ang), ring_cy + ring_ry * np.sin(ang)], float)
        score = abs(ang - base[i]) * 26.0
        for li in ent[i]["lig"]:
            li = int(li)
            lx, ly = float(ml + xy[li].x), float(mt + xy[li].y)
            v = np.array([lx - node[0], ly - node[1]], float)
            vn = float(np.hypot(v[0], v[1]))
            s = node if vn < 1e-6 else node + br * v / vn
            t = np.array([lx, ly], float) - (12.0 * v / max(vn, 1e-6))
            for ai, ap in atom_geom:
                if ai == li:
                    continue
                da = _seg_point_dist(s, t, ap)
                if da < 14.0:
                    score += (14.0 - da) ** 2 * 0.9
            for bp, bq, bi, bj in bond_geom:
                db = _seg_seg_dist(s, t, bp, bq)
                bsoft = 5.2 if li in (bi, bj) else 8.2
                if db < bsoft:
                    score += (bsoft - db) ** 2 * (5.6 if li in (bi, bj) else 9.2)
        return score

    for _ in range(4):
        for i in range(n):
            best_t = float(theta[i])
            best_s = _score(i, best_t)
            for dt in dts[1:]:
                cand = float(theta[i] + dt)
                ok = True
                for j in range(n):
                    if j == i:
                        continue
                    dd = abs(np.arctan2(np.sin(cand - theta[j]), np.cos(cand - theta[j])))
                    if dd < dmin * 0.84:
                        ok = False
                        break
                if not ok:
                    continue
                sc = _score(i, cand)
                if sc < best_s:
                    best_s, best_t = sc, cand
            theta[i] = best_t

    for i, e in enumerate(ent):
        t = float(theta[i])
        e["sx"] = float(ring_cx + ring_rx * np.cos(t))
        e["sy"] = float(ring_cy + ring_ry * np.sin(t))
        e["t"] = t


def _apply_manual_clock_overrides(
    ent: list[dict[str, object]],
    ring_cx: float,
    ring_cy: float,
    ring_rx: float,
    ring_ry: float,
) -> None:
    # Clock mapping on display:
    # 12 -> -pi/2, 3 -> 0, 6 -> pi/2, 9 -> pi.
    clock_angles = {
        "K103:O": (5.0 * np.pi) / 6.0,  # 8 o'clock
        "K103:N": (2.0 * np.pi) / 3.0,  # 7 o'clock
        "Y181:CD1": np.pi / 2.0,  # 6 o'clock
        "Y188:C": (5.0 * np.pi) / 12.0,  # 5:30
        "V189:C": np.pi / 3.0,  # 5 o'clock
        "V179:CB": np.pi / 4.0,  # 4:30
        "V179:CG1": np.pi / 6.0,  # 4 o'clock
        "Y188:CB": 0.0,  # 3 o'clock
    }
    for e in ent:
        label = str(e["label"])
        if label not in clock_angles:
            continue
        t = float(clock_angles[label])
        e["sx"] = float(ring_cx + ring_rx * np.cos(t))
        e["sy"] = float(ring_cy + ring_ry * np.sin(t))
        e["t"] = t


def _place_contact_nodes(
    ent: list[dict[str, object]],
    anchor_cx: float,
    anchor_cy: float,
    ring_cx: float,
    ring_cy: float,
    ring_rx: float,
    ring_ry: float,
    br: float,
    avoid_pts: np.ndarray | None = None,
) -> None:
    n = len(ent)
    if n == 0:
        return

    ang = np.array([np.arctan2(float(e["py"]) - anchor_cy, float(e["px"]) - anchor_cx) for e in ent], float)
    order = np.argsort(ang)
    theta = ang[order].copy()

    reff = max(1.0, min(ring_rx, ring_ry))
    dmin = min((2.0 * np.pi / n) * 0.94, (2.0 * br + 20.0) / reff)
    for _ in range(18):
        for i in range(1, n):
            gap = theta[i] - theta[i - 1]
            if gap < dmin:
                shift = 0.5 * (dmin - gap)
                theta[i - 1] -= shift
                theta[i] += shift
        wrap_gap = (theta[0] + 2.0 * np.pi) - theta[-1]
        if wrap_gap < dmin:
            shift = 0.5 * (dmin - wrap_gap)
            theta[0] += shift
            theta[-1] -= shift

    mean0, mean1 = float(np.mean(ang[order])), float(np.mean(theta))
    theta += (mean0 - mean1)
    out_theta = np.empty_like(theta)
    out_theta[order] = theta

    by_res = defaultdict(list)
    for i, e in enumerate(ent):
        by_res[str(e["label"]).split(":")[0]].append(i)
    offset_step = min(0.16, dmin * 0.42)
    for idxs in by_res.values():
        if len(idxs) < 2:
            continue
        idxs = sorted(idxs, key=lambda i: str(ent[i]["label"]))
        mid = 0.5 * (len(idxs) - 1)
        for k, idx in enumerate(idxs):
            out_theta[idx] += (k - mid) * offset_step

    theta2 = np.sort(out_theta.copy())
    for _ in range(10):
        for i in range(1, n):
            gap = theta2[i] - theta2[i - 1]
            if gap < dmin * 0.92:
                shift = 0.5 * (dmin * 0.92 - gap)
                theta2[i - 1] -= shift
                theta2[i] += shift
        wrap_gap = (theta2[0] + 2.0 * np.pi) - theta2[-1]
        if wrap_gap < dmin * 0.92:
            shift = 0.5 * (dmin * 0.92 - wrap_gap)
            theta2[0] += shift
            theta2[-1] -= shift
    idx_sorted = np.argsort(out_theta)
    out_theta[idx_sorted] = theta2

    if avoid_pts is not None and avoid_pts.size:
        base_theta = out_theta.copy()
        avoid_clear = br + 44.0
        node_clear = 2.0 * br + 14.0
        for _ in range(26):
            moved = False
            for i in range(n):
                t = float(out_theta[i])
                p = np.array([ring_cx + ring_rx * np.cos(t), ring_cy + ring_ry * np.sin(t)], float)
                d0 = float(np.min(np.hypot(p[0] - avoid_pts[:, 0], p[1] - avoid_pts[:, 1])))
                if d0 >= avoid_clear:
                    continue

                best_t = t
                best_score = d0 - (abs(t - base_theta[i]) * reff * 0.35)
                for dt in (-0.10, -0.07, -0.04, 0.04, 0.07, 0.10):
                    tc = t + dt
                    pc = np.array([ring_cx + ring_rx * np.cos(tc), ring_cy + ring_ry * np.sin(tc)], float)
                    dav = float(np.min(np.hypot(pc[0] - avoid_pts[:, 0], pc[1] - avoid_pts[:, 1])))
                    pen = 0.0
                    for j in range(n):
                        if j == i:
                            continue
                        tj = float(out_theta[j])
                        pj = np.array([ring_cx + ring_rx * np.cos(tj), ring_cy + ring_ry * np.sin(tj)], float)
                        dj = float(np.hypot(pc[0] - pj[0], pc[1] - pj[1]))
                        if dj < node_clear:
                            pen += (node_clear - dj) * 1.25
                    score = dav - pen - (abs(tc - base_theta[i]) * reff * 0.35)
                    if score > best_score + 0.2:
                        best_score = score
                        best_t = tc
                if best_t != t:
                    out_theta[i] = best_t
                    moved = True

            idx = np.argsort(out_theta)
            ts = out_theta[idx].copy()
            dkeep = dmin * 0.9
            for _ in range(6):
                for i in range(1, n):
                    gap = ts[i] - ts[i - 1]
                    if gap < dkeep:
                        shift = 0.5 * (dkeep - gap)
                        ts[i - 1] -= shift
                        ts[i] += shift
                wrap_gap = (ts[0] + 2.0 * np.pi) - ts[-1]
                if wrap_gap < dkeep:
                    shift = 0.5 * (dkeep - wrap_gap)
                    ts[0] += shift
                    ts[-1] -= shift
            out_theta[idx] = ts
            if not moved:
                break

    for i, e in enumerate(ent):
        t = float(out_theta[i])
        e["sx"] = float(ring_cx + ring_rx * np.cos(t))
        e["sy"] = float(ring_cy + ring_ry * np.sin(t))
        e["t"] = t


def plot(cif: Path, defs_csv: Path, out: Path) -> None:
    df = pd.read_csv(defs_csv)
    if df.empty:
        raise ValueError("Empty contact definition table")
    df["category"] = df["category"].astype(str).str.lower()
    df["protein_resid_auth"] = pd.to_numeric(df["protein_resid_auth"], errors="coerce")
    df = df.dropna(subset=["ligand_atom", "protein_resid_auth", "category"]).copy()
    # Omit Y188:CB:F aromatic contact connector for clarity in this figure.
    df = df[
        ~(
            (df["protein_resid_auth"] == 188)
            & (df["protein_atom"].astype(str).str.strip().str.upper() == "CB")
            & (df["ligand_atom"].astype(str).str.strip().str.upper().isin(["F"]))
        )
    ].copy()
    df["protein_label"] = df.apply(lambda r: f"{AA1.get(str(r['protein_resname']).upper(), str(r['protein_resname'])[:1])}{int(r['protein_resid_auth'])}:{r['protein_atom']}", axis=1)

    mol, atom_map = _ligand_from_cif(cif)
    w, h = 1800, 1200
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    dopts = d.drawOptions()
    dopts.bondLineWidth = 2.2
    dopts.clearBackground = False

    lig_rows = defaultdict(list)
    for r in df.itertuples(index=False):
        lig_rows[str(r.ligand_atom)].append(r)

    h_atoms, h_cols, h_rad = [], {}, {}
    for lig, rows in lig_rows.items():
        i = atom_map.get(lig)
        if i is None:
            continue
        cats = {str(r.category) for r in rows}
        c = _cat(cats)
        h_atoms.append(i)
        h_cols[i] = COL[c][2]
        h_rad[i] = 0.48

    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=h_atoms, highlightAtomColors=h_cols, highlightAtomRadii=h_rad)
    xy = {i: d.GetDrawCoords(i) for i in h_atoms}
    all_xy = {i: d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())}
    d.FinishDrawing()

    img = plt.imread(BytesIO(d.GetDrawingText()), format="png")
    ml, mr, mt, mb = 320, 320, 220, 240
    cw, ch = w + ml + mr, h + mt + mb
    fig, ax = plt.subplots(figsize=(13.5, 9.2))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.imshow(img, extent=(ml, ml + w, mt + h, mt), zorder=2)
    ax.set_xlim(0, cw)
    ax.set_ylim(ch, 0)
    ax.axis("off")

    x = np.array([ml + xy[i].x for i in h_atoms], float)
    y = np.array([mt + xy[i].y for i in h_atoms], float)
    cx, cy = float(x.mean()), float(y.mean())
    sx, sy = max(280.0, np.ptp(x) * 1.2), max(220.0, np.ptp(y) * 1.35)
    gcx, gcy = cw / 2.0, ch / 2.0
    pocket_w, pocket_h = sx * 1.48, sy * 1.54
    ax.add_patch(
        Ellipse(
            (gcx, gcy),
            pocket_w,
            pocket_h,
            facecolor=(0.53, 0.83, 0.58, 0.16),
            edgecolor=(0.12, 0.68, 0.35, 0.95),
            linewidth=2.2,
            zorder=1,
        )
    )

    groups = {}
    for r in df.itertuples(index=False):
        lig = atom_map.get(str(r.ligand_atom))
        if lig is None or lig not in xy:
            continue
        g = groups.setdefault(str(r.protein_label), {"cats": set(), "lig": set()})
        g["cats"].add(str(r.category).lower())
        g["lig"].add(int(lig))

    ent = []
    for label, g in groups.items():
        pts = np.array([(ml + xy[i].x, mt + xy[i].y) for i in sorted(g["lig"])], float)
        ent.append({"label": label, "cat": _cat(g["cats"]), "lig": sorted(g["lig"]), "px": float(pts[:, 0].mean()), "py": float(pts[:, 1].mean())})

    ent.sort(key=lambda e: np.arctan2(float(e["py"]) - cy, float(e["px"]) - cx))
    n = len(ent)
    br = 46.0
    if n:
        ring_rx = max(260.0, pocket_w * 0.5 - (br + 7.0))
        ring_ry = max(220.0, pocket_h * 0.5 - (br + 7.0))
        avoid_idx = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6)]
        avoid_pts = np.array([(ml + all_xy[i].x, mt + all_xy[i].y) for i in avoid_idx], float)
        _place_contact_nodes(ent, cx, cy, gcx, gcy, ring_rx, ring_ry, br, avoid_pts=avoid_pts)
        _refine_node_angles_for_connector_clearance(ent, gcx, gcy, ring_rx, ring_ry, br, xy, all_xy, ml, mt, mol)
        _apply_manual_clock_overrides(ent, gcx, gcy, ring_rx, ring_ry)

    atom_geom = [(i, np.array([float(ml + all_xy[i].x), float(mt + all_xy[i].y)], float)) for i in range(mol.GetNumAtoms())]
    bond_geom = []
    for b in mol.GetBonds():
        bi, bj = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        bp = np.array([float(ml + all_xy[bi].x), float(mt + all_xy[bi].y)], float)
        bq = np.array([float(ml + all_xy[bj].x), float(mt + all_xy[bj].y)], float)
        bond_geom.append((bp, bq, bi, bj))

    center = np.array([gcx, gcy], float)
    con_specs = []
    for e in ent:
        _, edge, _ = COL[e["cat"]]
        node = np.array([float(e["sx"]), float(e["sy"])], float)
        rv = center - node
        rn = float(np.hypot(rv[0], rv[1]))
        rin = np.array([0.0, 1.0], float) if rn < 1e-6 else rv / rn
        tan = np.array([-rin[1], rin[0]], float)
        force_straight = str(e["label"]) in {"K103:O", "K103:N"}
        for li in e["lig"]:
            lx, ly = float(ml + xy[li].x), float(mt + xy[li].y)
            v = np.array([lx - node[0], ly - node[1]], float)
            vn = float(np.hypot(v[0], v[1]))
            s = node if vn < 1e-6 else node + br * v / vn
            t = np.array([lx, ly], float) - (12.0 * v / max(vn, 1e-6))
            con_specs.append(
                {"li": int(li), "color": edge, "s": s, "t": t, "rin": rin, "tan": tan, "force_straight": force_straight}
            )

    con_specs.sort(key=lambda c: -float(np.hypot(*(c["t"] - c["s"]))))
    used_segs: list[tuple[np.ndarray, np.ndarray]] = []
    for c in con_specs:
        s = c["s"]
        t = c["t"]
        rin = c["rin"]
        tan = c["tan"]
        li = int(c["li"])
        force_straight = bool(c.get("force_straight", False))
        if force_straight:
            best = [s, t]
            ax.plot([pt[0] for pt in best], [pt[1] for pt in best], color=c["color"], linewidth=1.35, alpha=0.88, zorder=3)
            for p, q in zip(best[:-1], best[1:]):
                used_segs.append((np.array(p, float), np.array(q, float)))
            continue
        candidates = [[s, t]]
        for off_tan, off_r in ((24.0, 14.0), (38.0, 20.0), (56.0, 28.0)):
            candidates.append([s, s + tan * off_tan + rin * off_r, t])
            candidates.append([s, s - tan * off_tan + rin * off_r, t])
            candidates.append([s, s + tan * off_tan + rin * off_r, t + tan * (0.55 * off_tan), t])
            candidates.append([s, s - tan * off_tan + rin * off_r, t - tan * (0.55 * off_tan), t])
            # Outward doglegs can avoid bond overlays near the ligand core.
            candidates.append([s, s + tan * off_tan - rin * (0.55 * off_r), t])
            candidates.append([s, s - tan * off_tan - rin * (0.55 * off_r), t])
            candidates.append([s, s + tan * off_tan - rin * (0.55 * off_r), t + tan * (0.4 * off_tan), t])
            candidates.append([s, s - tan * off_tan - rin * (0.55 * off_r), t - tan * (0.4 * off_tan), t])

        atom_hard = 8.0
        bond_hard = 6.5
        conn_hard = 2.6

        straight = candidates[0]
        if not _path_has_hard_clash(straight, li, atom_geom, bond_geom, used_segs, atom_hard, bond_hard, conn_hard):
            best = straight
        else:
            valid = [
                path
                for path in candidates[1:]
                if not _path_has_hard_clash(path, li, atom_geom, bond_geom, used_segs, atom_hard, bond_hard, conn_hard)
            ]
            if valid:
                best = min(
                    valid,
                    key=lambda path: _path_soft_score(path, li, atom_geom, bond_geom, used_segs),
                )
            else:
                # Fallback: keep rendering even if every candidate collides; choose least-bad.
                best = min(
                    candidates,
                    key=lambda path: _path_soft_score(path, li, atom_geom, bond_geom, used_segs),
                )

        ax.plot([pt[0] for pt in best], [pt[1] for pt in best], color=c["color"], linewidth=1.35, alpha=0.88, zorder=3)
        for p, q in zip(best[:-1], best[1:]):
            used_segs.append((np.array(p, float), np.array(q, float)))

    for e in ent:
        fill, edge, _ = COL[e["cat"]]
        ax.add_patch(Circle((e["sx"], e["sy"]), br + 1.8, facecolor="white", edgecolor="none", alpha=0.36, zorder=3.7))
        ax.add_patch(Circle((e["sx"], e["sy"]), br, facecolor=fill, edgecolor=edge, linewidth=1.9, zorder=4))
        res, atm = (e["label"].split(":") + [""])[:2]
        ax.text(e["sx"], e["sy"] - 8, res, ha="center", va="center", fontsize=9.8, fontweight="bold", color="#1f2328", zorder=5)
        ax.text(e["sx"], e["sy"] + 12, atm, ha="center", va="center", fontsize=8.3, color="#2f3842", zorder=5)

    inv = {v: k for k, v in atom_map.items()}
    for i in h_atoms:
        lx, ly = float(ml + xy[i].x), float(mt + xy[i].y)
        ax.add_patch(Circle((lx, ly), 10.8, facecolor=h_cols[i], edgecolor="#1f2933", linewidth=1.1, zorder=7))
        ax.text(lx, ly, inv.get(i, ""), ha="center", va="center", fontsize=8.3, color="white", fontweight="bold", zorder=8)

    legend_keys = [k for k in ("polar", "hydrophobic", "mixed") if any(str(e["cat"]) == k for e in ent)]
    if legend_keys:
        ly = 66.0
        entry_w = 170.0
        box_w = 56.0 + entry_w * len(legend_keys)
        lx0 = (cw - box_w) / 2.0
        ax.add_patch(
            FancyBboxPatch(
                (lx0, ly - 28.0),
                box_w,
                58.0,
                boxstyle="round,pad=0.35,rounding_size=16",
                facecolor=(1.0, 1.0, 1.0, 0.82),
                edgecolor=(0.72, 0.78, 0.84, 0.95),
                linewidth=1.1,
                zorder=10,
            )
        )
        for j, k in enumerate(legend_keys):
            fill, edge, _ = COL[k]
            x0 = lx0 + 30.0 + j * entry_w
            ax.add_patch(Circle((x0, ly), 10.6, facecolor=fill, edgecolor=edge, linewidth=1.45, zorder=11))
            ax.text(x0 + 17.0, ly, k.capitalize(), ha="left", va="center", fontsize=10.1, color="#1f2933", zorder=11)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.99)
    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print(f"Saved: {out}")


def main() -> int:
    p = argparse.ArgumentParser(description="Plot Doravirine 2D crystal-contact diagram.")
    p.add_argument("--cif", type=Path, default=Path("data/structures/4NCG.cif"))
    p.add_argument("--contact-defs", type=Path, default=Path("results/dor_key_contact_definitions_4ncg.csv"))
    p.add_argument("--output", type=Path, default=Path("results/plots/crystal_dor_key_contacts_2d.png"))
    a = p.parse_args()
    if not a.cif.exists():
        raise FileNotFoundError(f"Missing CIF: {a.cif}")
    if not a.contact_defs.exists():
        raise FileNotFoundError(f"Missing contact definition table: {a.contact_defs}")
    plot(a.cif, a.contact_defs, a.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
