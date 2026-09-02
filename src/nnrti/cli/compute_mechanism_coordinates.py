#!/usr/bin/env python3
"""Trajectory coordinates for the molecular-mechanism section.

Computes a battery of candidate order parameters in a single pass per
replicate, so the most contrastive coordinate can be chosen afterwards:

  Y188L   aromatic stacking between residue 188 and the DOR chlorocyanophenyl
          ring (centroid distance, interplanar angle, minimum contact)
  V106A   dislocation of DOR toward Ser105
  K103N   polar contact between the residue 103 side chain and DOR, plus the
          conserved backbone hydrogen bond to the triazolinone
  G190E   Val179 packing, residue 190 side-chain contact, and Glu190 salt
          bridges to Lys101/Lys103

DOR ring systems are identified from bond connectivity rather than atom names
(the MD topology renames the crystallographic ligand atoms).

Residue numbering follows the repo convention: topology resid = canonical - 3.

Output: results/analysis/mechanisms/mechanism_coordinates.csv (per-frame rows)
"""
from __future__ import annotations

import argparse
import itertools
import json
import logging
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

LIG_RESNAME = "2KW"
P66_CHAIN = 0  # chain 0 is p66; p51 reuses the same residue numbering
OFFSET = -3  # topology resid = canonical + OFFSET

# canonical residue numbers
R188, R106, R105, R103, R101, R179, R190 = 188, 106, 105, 103, 101, 179, 190
R229 = 229  # Trp229: conserved pocket aromatic, never mutated in this panel
R227 = 227
CONTACT_CUT = 4.0  # A  (switched from 4.5 on 2026-08-31; see paper/contact-cutoff-sensitivity.md)

TYR_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
TRP_RING = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def _remap(candidate: Path, repo_root: Path) -> Path:
    if candidate.exists():
        return candidate
    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker in text:
        mapped = repo_root / text.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    return candidate


def _replicate_inputs(row: pd.Series, repo_root: Path) -> tuple[Path, Path, Path]:
    js = _remap(Path(str(row["output_json"])), repo_root)
    data = json.loads(js.read_text())
    topo = _remap(Path(str(data.get("analysis_topology_pdb") or "").strip()), repo_root)
    dcd = _remap(Path(str(data.get("analysis_dcd") or "").strip()), repo_root)
    if not topo.exists() or not dcd.exists():
        raise FileNotFoundError(f"missing files for {row['mutation']} rep{row['replicate']}")
    return topo, dcd, js


def _total_ns(js: Path) -> float:
    try:
        j = json.loads(js.read_text())
        steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
        if steps > 0:
            return steps * 2.0 / 1_000_000.0
    except Exception:
        pass
    m = re.match(r"^(.+)_rep(\d{2})\.json$", js.name)
    if m:
        csv = js.parent / f"{m.group(1)}_rep{m.group(2)}_md_state.csv"
        if csv.exists():
            try:
                sdf = pd.read_csv(csv)
                for col in ('#"Step"', "Step"):
                    if col in sdf.columns:
                        s = pd.to_numeric(sdf[col], errors="coerce").dropna()
                        if not s.empty:
                            return float(s.max()) * 2.0 / 1_000_000.0
            except Exception:
                pass
    return 100.0


def _ligand_rings(traj) -> dict:
    """Identify DOR ring systems from connectivity in the first frame."""
    top = traj.topology
    lig = [a for a in top.atoms if a.residue.name == LIG_RESNAME and a.element.symbol != "H"]
    if not lig:
        raise ValueError("no ligand heavy atoms found")
    idx = {a.name: a.index for a in lig}
    xyz = traj.xyz[0] * 10.0  # nm -> A
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(idx)
    for i, j in itertools.combinations(idx, 2):
        if np.linalg.norm(xyz[idx[i]] - xyz[idx[j]]) < 1.85:
            G.add_edge(i, j)
    rings = [sorted(c) for c in nx.cycle_basis(G) if len(c) >= 5]

    out: dict = {"all_heavy": [idx[n] for n in idx]}
    for r in rings:
        nbrs = {n for a in r for n in G[a] if n not in r}
        if any(n.startswith("Cl") for n in nbrs):
            out["ring_chlorophenyl"] = [idx[a] for a in r]
        elif len(r) == 5:
            out["ring_triazolinone"] = [idx[a] for a in r]
        else:
            out["ring_pyridinone"] = [idx[a] for a in r]
    out["polar"] = [idx[n] for n in idx if n[0] in ("N", "O")]
    # aryl ether torsion: pyridinone C - O(ether) - phenyl C - phenyl C
    ether = [n for n in idx if n.startswith("O") and len(list(G[n])) == 2]
    for e in ether:
        nb = list(G[e])
        a_in_ph = [x for x in nb if idx[x] in out.get("ring_chlorophenyl", [])]
        a_in_py = [x for x in nb if idx[x] in out.get("ring_pyridinone", [])]
        if a_in_ph and a_in_py:
            ph = a_in_ph[0]
            ph_nb = [x for x in G[ph] if idx[x] in out.get("ring_chlorophenyl", [])][0]
            out["tau2_atoms"] = [idx[a_in_py[0]], idx[e], idx[ph], idx[ph_nb]]
            break
    out["triaz_N"] = [
        idx[n] for n in idx if n.startswith("N") and idx[n] in out.get("ring_triazolinone", [])
    ]
    return out


def _res_atoms(top, canonical: int, sidechain_only: bool = False, names=None):
    """Atoms of a canonical p66 residue.

    p66 and p51 share the same residue numbering, so the chain must be pinned
    or every selection silently picks up the p51 copy as well.
    """
    resid = canonical + OFFSET
    sel = []
    for a in top.atoms:
        if a.residue.name == LIG_RESNAME:
            continue
        if a.residue.chain.index != P66_CHAIN:
            continue
        if a.residue.resSeq != resid:
            continue
        if a.element.symbol == "H":
            continue
        if names is not None and a.name not in names:
            continue
        if sidechain_only and a.name in BACKBONE:
            continue
        sel.append(a.index)
    return sel


def _centroid(xyz, idx):
    return xyz[idx].mean(axis=0)


def _normal(xyz, idx):
    pts = xyz[idx]
    return np.linalg.svd(pts - pts.mean(axis=0))[2][2]


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def _ncontacts(xyz, a, b, cut=CONTACT_CUT):
    if not a or not b:
        return np.nan
    d = np.linalg.norm(xyz[a][:, None, :] - xyz[b][None, :, :], axis=-1)
    return int((d < cut).sum())


def _mindist(xyz, a, b):
    if not a or not b:
        return np.nan
    d = np.linalg.norm(xyz[a][:, None, :] - xyz[b][None, :, :], axis=-1)
    return float(d.min())


def _process(row: pd.Series, repo_root: Path, stride: int) -> list[dict]:
    import mdtraj as md

    topo, dcd, js = _replicate_inputs(row, repo_root)
    traj = md.load(str(dcd), top=str(topo), stride=max(1, stride))
    traj.make_molecules_whole(inplace=True)
    top = traj.topology

    rings = _ligand_rings(traj)
    chl = rings.get("ring_chlorophenyl", [])
    triazN = rings.get("triaz_N", [])
    lig_all = rings["all_heavy"]
    lig_polar = rings["polar"]
    tau2 = rings.get("tau2_atoms")

    r188_sc = _res_atoms(top, R188, sidechain_only=True)
    r188_ring = _res_atoms(top, R188, names=TYR_RING)
    is_tyr188 = len(r188_ring) == 6

    r106_sc = _res_atoms(top, R106, sidechain_only=True)
    r105_og = _res_atoms(top, R105, names={"OG"}) or _res_atoms(top, R105, names={"CB"})
    r105_ca = _res_atoms(top, R105, names={"CA"})
    r105_sc = _res_atoms(top, R105, sidechain_only=True)

    r103_sc_polar = _res_atoms(top, R103, names={"NZ", "OD1", "ND2"})
    r103_sc = _res_atoms(top, R103, sidechain_only=True)
    r103_bb = _res_atoms(top, R103, names={"N", "O"})

    r229_ring = _res_atoms(top, R229, names=set(TRP_RING))
    r229_sc = _res_atoms(top, R229, sidechain_only=True)
    r227_sc = _res_atoms(top, R227, sidechain_only=True)
    prot_heavy = [
        a.index for a in top.atoms
        if a.residue.name != LIG_RESNAME and a.element.symbol != "H"
    ]
    r101_nz = _res_atoms(top, R101, names={"NZ"})
    r179_sc = _res_atoms(top, R179, sidechain_only=True)
    r190_sc = _res_atoms(top, R190, sidechain_only=True)
    r190_ca = _res_atoms(top, R190, names={"CA"})
    r190_carboxyl = _res_atoms(top, R190, names={"OE1", "OE2"})

    total = _total_ns(js)
    n = traj.n_frames
    out = []
    for f in range(n):
        xyz = traj.xyz[f] * 10.0
        rec = {
            "mutation": str(row["mutation"]),
            "safe_label": str(row.get("safe_label", "")),
            "replicate": int(row["replicate"]),
            "frame": f,
            "time_ns": (f / max(1, n - 1)) * total,
        }
        # --- Y188L : stacking ---
        if chl and r188_sc:
            rec["y188_cent_dist"] = float(
                np.linalg.norm(_centroid(xyz, r188_sc) - _centroid(xyz, chl))
            )
            rec["y188_mindist"] = _mindist(xyz, r188_sc, chl)
        if is_tyr188 and chl:
            nrm_a, nrm_b = _normal(xyz, r188_ring), _normal(xyz, chl)
            rec["y188_interplanar_deg"] = float(
                np.degrees(np.arccos(np.clip(abs(np.dot(nrm_a, nrm_b)), 0, 1)))
            )
            rec["y188_ringcent_dist"] = float(
                np.linalg.norm(_centroid(xyz, r188_ring) - _centroid(xyz, chl))
            )
        if chl and r229_ring:
            nA, nB = _normal(xyz, chl), _normal(xyz, r229_ring)
            rec["chl_trp229_interplanar_deg"] = float(
                np.degrees(np.arccos(np.clip(abs(np.dot(nA, nB)), 0, 1)))
            )
            rec["chl_trp229_cent_dist"] = float(
                np.linalg.norm(_centroid(xyz, chl) - _centroid(xyz, r229_ring))
            )
        if tau2:
            rec["dor_tau2_deg"] = abs(
                _dihedral(xyz[tau2[0]], xyz[tau2[1]], xyz[tau2[2]], xyz[tau2[3]])
            )
        if r188_sc:
            rec["n_contacts_188_dor"] = _ncontacts(xyz, r188_sc, lig_all)
            rec["n_contacts_188_chl"] = _ncontacts(xyz, r188_sc, chl)
        if chl:
            rec["chl_ring_burial"] = _ncontacts(xyz, chl, prot_heavy)
        rec["res227_to_dor_mindist"] = _mindist(xyz, r227_sc, lig_all)
        rec["res229_to_dor_mindist"] = _mindist(xyz, r229_sc, lig_all)
        # --- V106A : dislocation toward Ser105 ---
        if r105_og:
            rec["dor_cent_to_S105OG"] = float(
                np.linalg.norm(_centroid(xyz, lig_all) - _centroid(xyz, r105_og))
            )
        if r105_ca:
            rec["dor_cent_to_S105CA"] = float(
                np.linalg.norm(_centroid(xyz, lig_all) - _centroid(xyz, r105_ca))
            )
        rec["dor_to_S105_mindist"] = _mindist(xyz, lig_all, r105_sc)
        if r106_sc:
            rec["dor_to_res106_mindist"] = _mindist(xyz, lig_all, r106_sc)
            rec["dor_cent_to_res106"] = float(
                np.linalg.norm(_centroid(xyz, lig_all) - _centroid(xyz, r106_sc))
            )
        # --- K103N : polar contact ---
        rec["res103_sc_to_dor_polar"] = _mindist(xyz, r103_sc_polar, lig_polar)
        rec["res103_sc_to_dor"] = _mindist(xyz, r103_sc, lig_all)
        rec["res103_bb_to_triazN"] = _mindist(xyz, r103_bb, triazN)
        # --- G190E ---
        rec["res179_to_dor_mindist"] = _mindist(xyz, r179_sc, lig_all)
        rec["res190_to_dor_mindist"] = _mindist(xyz, r190_sc, lig_all)
        rec["res190CA_to_dor_mindist"] = _mindist(xyz, r190_ca, lig_all)
        if r190_carboxyl and r101_nz:
            rec["res190_to_K101NZ"] = _mindist(xyz, r190_carboxyl, r101_nz)
        out.append(rec)
    return out


def main() -> int:
    warnings.filterwarnings("ignore")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/analysis/mechanisms/mechanism_coordinates.csv"),
    )
    args = p.parse_args()

    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(args.manifest)
    if args.mutations:
        mf = mf[mf["mutation"].isin(args.mutations)].copy()
    if mf.empty:
        logging.error("no manifest rows selected")
        return 1

    rows: list[dict] = []
    for _, r in mf.iterrows():
        try:
            got = _process(r, repo_root, args.stride)
            rows.extend(got)
            logging.info(f"{r['mutation']} rep{int(r['replicate'])}: {len(got)} frames")
        except Exception as exc:
            logging.error(f"FAILED {r['mutation']} rep{int(r['replicate'])}: {exc}")

    if not rows:
        return 1
    df = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    logging.info(f"wrote {args.output_csv} ({len(df)} rows, {df.shape[1]} cols)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
