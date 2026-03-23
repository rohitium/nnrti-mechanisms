#!/usr/bin/env python3
"""Compute doravirine torsion angles from a crystal mmCIF (default: 4NCG).

This script reports the same four torsions used by the MD trajectory analysis
pipeline, but for a single crystal pose.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser


REPO = Path(__file__).resolve().parents[3]

# MD-analysis torsion definitions (x-suffixed atom names)
TORSIONS_X: list[tuple[str, tuple[str, str, str, str]]] = [
    ("tau1", ("C12x", "C2x", "O1x", "C9x")),
    ("tau2", ("C2x", "O1x", "C9x", "C10x")),
    ("tau3", ("C4x", "N2x", "C15x", "C14x")),
    ("tau4", ("N2x", "C15x", "C14x", "N5x")),
]

# Mapping from MD x-suffixed names -> 4NCG ligand atom IDs (2KW).
# This mapping was derived by graph isomorphism (element + bond topology)
# between the analysis topology ligand and the 4NCG chem_comp definition.
X_TO_CIF = {
    "C12x": "C12",
    "C2x": "C7",
    "O1x": "O",
    "C9x": "C2",
    "C10x": "C1",
    "C4x": "C10",
    "N2x": "N11",
    "C15x": "C17",
    "C14x": "C18",
    "N5x": "N22",
}


def _calc_dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Return dihedral angle in degrees in the range (-180, 180]."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(y, x)))


def _select_atom(residue, atom_name: str):
    """Pick atom by name, resolving disordered altloc to A/blank/max-occupancy."""
    if atom_name not in residue:
        raise KeyError(f"Atom {atom_name} not found in residue {residue.get_resname()} {residue.id}")
    atom = residue[atom_name]
    if not getattr(atom, "is_disordered", lambda: 0)():
        return atom
    alts = list(atom.disordered_get_list())
    if not alts:
        return atom
    for preferred in ("A", " "):
        for a in alts:
            if str(a.get_altloc()).strip() == preferred.strip():
                return a
    alts_sorted = sorted(alts, key=lambda a: float(a.get_occupancy() or 0.0), reverse=True)
    return alts_sorted[0]


def _pick_ligand_residue(structure, ligand_resname: str, chain_id: str, resid: int) -> tuple[int, str, int, object]:
    hits: list[tuple[int, str, int, object]] = []
    want_chain = str(chain_id).strip()
    want_resid = int(resid) if resid is not None else None

    for model in structure:
        for chain in model:
            if want_chain and str(chain.id).strip() != want_chain:
                continue
            for residue in chain:
                if str(residue.get_resname()).strip().upper() != ligand_resname.upper():
                    continue
                seqnum = int(residue.id[1])
                if want_resid is not None and seqnum != want_resid:
                    continue
                hits.append((int(model.id), str(chain.id), seqnum, residue))

    if not hits:
        raise ValueError(
            f"No residue found for ligand={ligand_resname}, chain={want_chain or '*'}, resid={want_resid if want_resid is not None else '*'}"
        )
    return hits[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute DOR torsion angles from crystal 4NCG.cif")
    parser.add_argument("--cif", type=Path, default=Path("data/structures/4NCG.cif"))
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--chain-id",
        type=str,
        default="",
        help="Optional chain ID for the ligand residue (default: first match).",
    )
    parser.add_argument(
        "--resid",
        type=int,
        default=None,
        help="Optional residue number for the ligand residue (default: first match).",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("results/tables/holo/dor_torsions_4ncg.csv"),
        help="Output CSV path.",
    )
    args = parser.parse_args()

    cif_path = args.cif
    if not cif_path.exists():
        raise FileNotFoundError(f"Missing CIF: {cif_path}")

    parser_obj = MMCIFParser(QUIET=True)
    structure = parser_obj.get_structure("crystal", str(cif_path))
    model_id, chain_id, ligand_resid, residue = _pick_ligand_residue(
        structure=structure,
        ligand_resname=str(args.ligand_resname),
        chain_id=str(args.chain_id),
        resid=args.resid,
    )

    values: dict[str, float] = {}
    atom_map_text: list[str] = []
    for tname, atoms_x in TORSIONS_X:
        atoms_cif = [X_TO_CIF[a] for a in atoms_x]
        coords = [np.asarray(_select_atom(residue, a).coord, dtype=float) for a in atoms_cif]
        values[tname] = _calc_dihedral(*coords)
        atom_map_text.append(f"{tname}: {'-'.join(atoms_x)} => {'-'.join(atoms_cif)}")

    row = {
        "structure": cif_path.stem,
        "cif_path": str(cif_path if cif_path.is_absolute() else (REPO / cif_path)),
        "model_id": int(model_id),
        "chain_id": str(chain_id),
        "ligand_resname": str(args.ligand_resname).upper(),
        "ligand_resid": int(ligand_resid),
        **values,
    }

    out_csv = args.out_csv
    if not out_csv.is_absolute():
        out_csv = REPO / out_csv
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out_csv, index=False)

    print(f"Saved {out_csv}")
    for k in ("tau1", "tau2", "tau3", "tau4"):
        print(f"{k}: {values[k]: .3f} deg")
    print("Atom mapping:")
    for s in atom_map_text:
        print(f"  {s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
