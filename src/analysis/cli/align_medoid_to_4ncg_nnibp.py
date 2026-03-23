#!/usr/bin/env python3
"""Align a medoid structure to 4NCG using NNIBP pocket C-alpha atoms.

This creates a transformed copy of the medoid where the NNIBP pocket is
best-fit aligned to the reference crystal pocket.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


# Canonical HIV-1 RT p66 NNIBP-lining residues used across this repo.
NNIBP_P66_CANONICAL = [100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318]


def _parse_resseq_csv(text: str) -> list[int]:
    out: list[int] = []
    for tok in str(text).replace(";", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        out.append(int(tok))
    if not out:
        raise ValueError("No residue numbers parsed.")
    return out


def _largest_protein_chain_index(topology) -> int:
    best_idx = -1
    best_count = -1
    for c in topology.chains:
        count = sum(1 for r in c.residues if r.is_protein)
        if count > best_count:
            best_count = count
            best_idx = int(c.index)
    if best_idx < 0:
        raise ValueError("No protein chain found.")
    return best_idx


def _pick_atoms(topology, chain_idx: int, resseqs: list[int], atom_name: str) -> tuple[np.ndarray, list[tuple[int, str]]]:
    idx: list[int] = []
    missing: list[tuple[int, str]] = []
    for rs in resseqs:
        sel = topology.select(f"chainid {int(chain_idx)} and protein and resSeq {int(rs)} and name {atom_name}")
        if sel.size != 1:
            missing.append((int(rs), f"found={int(sel.size)}"))
            continue
        idx.append(int(sel[0]))
    return np.asarray(idx, dtype=int), missing


def _paired_rmsd_nm(a_xyz: np.ndarray, b_xyz: np.ndarray) -> float:
    d2 = np.sum((a_xyz - b_xyz) ** 2, axis=1)
    return float(np.sqrt(np.mean(d2)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Align medoid pocket to 4NCG NNIBP pocket.")
    parser.add_argument("--medoid-pdb", type=Path, default=Path("results/average_structures/WT_medoid_structure.pdb"))
    parser.add_argument("--reference-cif", type=Path, default=Path("data/structures/4NCG.cif"))
    parser.add_argument(
        "--canonical-pocket-resseq",
        type=str,
        default=",".join(str(x) for x in NNIBP_P66_CANONICAL),
        help="Comma-separated canonical p66 residue numbers for pocket alignment.",
    )
    parser.add_argument(
        "--model-resseq-offset",
        type=int,
        default=-3,
        help="Offset applied to canonical residue numbering in model (medoid) structure.",
    )
    parser.add_argument("--atom-name", type=str, default="CA")
    parser.add_argument(
        "--model-chain-index",
        type=int,
        default=-1,
        help="Protein chain index in medoid (default: auto-detect largest protein chain).",
    )
    parser.add_argument(
        "--reference-chain-index",
        type=int,
        default=-1,
        help="Protein chain index in reference (default: auto-detect largest protein chain).",
    )
    parser.add_argument(
        "--output-pdb",
        type=Path,
        default=Path("__AUTO__"),
        help="Output aligned PDB. Default: <medoid_stem>_aligned_to_4NCG_nnibp.pdb",
    )
    args = parser.parse_args()

    import mdtraj as md

    medoid_pdb = args.medoid_pdb
    ref_cif = args.reference_cif
    if not medoid_pdb.exists():
        raise FileNotFoundError(f"Missing medoid PDB: {medoid_pdb}")
    if not ref_cif.exists():
        raise FileNotFoundError(f"Missing reference CIF: {ref_cif}")

    out_pdb = args.output_pdb
    if str(out_pdb) == "__AUTO__":
        out_pdb = medoid_pdb.with_name(f"{medoid_pdb.stem}_aligned_to_4NCG_nnibp.pdb")
    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    med = md.load(str(medoid_pdb))
    ref = md.load(str(ref_cif))

    canon_res = _parse_resseq_csv(args.canonical_pocket_resseq)
    model_res = [int(r + int(args.model_resseq_offset)) for r in canon_res]

    model_chain = int(args.model_chain_index)
    if model_chain < 0:
        model_chain = _largest_protein_chain_index(med.topology)
    ref_chain = int(args.reference_chain_index)
    if ref_chain < 0:
        ref_chain = _largest_protein_chain_index(ref.topology)

    model_idx, model_missing = _pick_atoms(med.topology, model_chain, model_res, args.atom_name)
    ref_idx, ref_missing = _pick_atoms(ref.topology, ref_chain, canon_res, args.atom_name)
    if model_missing or ref_missing:
        raise ValueError(
            "Could not build complete pocket atom mapping.\n"
            f"model_missing={model_missing}\n"
            f"ref_missing={ref_missing}"
        )
    if model_idx.size != ref_idx.size or model_idx.size < 3:
        raise ValueError(
            f"Pocket atom count mismatch/too small: model={model_idx.size}, ref={ref_idx.size}"
        )

    rmsd_before_nm = _paired_rmsd_nm(med.xyz[0, model_idx], ref.xyz[0, ref_idx])
    med_aligned = med[:]
    med_aligned.superpose(ref, atom_indices=model_idx, ref_atom_indices=ref_idx)
    rmsd_after_nm = _paired_rmsd_nm(med_aligned.xyz[0, model_idx], ref.xyz[0, ref_idx])
    med_aligned.save_pdb(str(out_pdb))

    mapping = []
    for c, m_res, mi, ri in zip(canon_res, model_res, model_idx.tolist(), ref_idx.tolist()):
        ma = med.topology.atom(int(mi))
        ra = ref.topology.atom(int(ri))
        mapping.append(
            {
                "canonical_resSeq": int(c),
                "model_resSeq": int(m_res),
                "model_resname": str(ma.residue.name),
                "reference_resname": str(ra.residue.name),
                "atom_name": str(args.atom_name),
            }
        )

    summary = {
        "medoid_pdb": str(medoid_pdb),
        "reference_cif": str(ref_cif),
        "output_pdb": str(out_pdb),
        "atom_name": str(args.atom_name),
        "model_chain_index": int(model_chain),
        "reference_chain_index": int(ref_chain),
        "model_resseq_offset": int(args.model_resseq_offset),
        "n_alignment_atoms": int(model_idx.size),
        "rmsd_before_angstrom": float(rmsd_before_nm * 10.0),
        "rmsd_after_angstrom": float(rmsd_after_nm * 10.0),
        "pocket_mapping": mapping,
    }
    summary_json = out_pdb.with_suffix(".summary.json")
    summary_json.write_text(json.dumps(summary, indent=2))

    print(f"Saved {out_pdb}")
    print(f"Saved {summary_json}")
    print(f"Pocket RMSD before: {rmsd_before_nm * 10.0:.3f} A")
    print(f"Pocket RMSD after:  {rmsd_after_nm * 10.0:.3f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
