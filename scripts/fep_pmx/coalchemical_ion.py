#!/usr/bin/env python3
"""Co-alchemical ion for charge-changing pmx NEQ FEP legs.

A mutation that changes the protein's net charge (e.g. K103N: Lys+->Asn0, or
G190E: Gly0->Glu-, both Delta_q = -1) leaves the periodic box non-neutral at
lambda=1, which contaminates DeltaG through the Ewald/finite-size term. The
co-alchemical-ion fix (King et al., JMB 2019): pick one counter-ion in bulk
solvent and make it a dual-state particle whose charge changes by -Delta_q over
the same lambda, so the total box charge is invariant at every lambda.

For our legs Delta_q = -1, so we need +1 of compensation -> decouple one Cl-
(charge -1 -> 0, LJ -> off via the switch mdp's gapsys soft-core). The chosen
Cl- is the most bulk-solvated one (max of its minimum distance to any protein
atom), so its decoupling free energy is bulk-like and cancels between the holo
and apo legs in DeltaDeltaG = DeltaG_holo - DeltaG_apo.

This edits an already-neutralized GROMACS system (system.gro + system.top from
`build_solvated_system`): it moves the chosen ion to the end of the coordinate
file and rewrites one `CL`/`NA` entry in `[ molecules ]` into a dual-state
`*_coalch` moleculetype (plus a zero-LJ dummy atomtype).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Non-protein residue names to exclude when locating "bulk" (protein-distant) ions.
_NON_PROTEIN = {"SOL", "WAT", "HOH", "CL", "NA", "K", "MG", "CA", "RB", "IB", "2KW"}

# GROMACS ion species we know how to decouple: resname -> (atomtype, mass, sigma, eps, charge)
_ION_SPEC = {
    "CL": ("Cl", 35.45, "4.47766e-01", "1.48913e-01", -1.0),
    "NA": ("Na", 22.99, "2.43928e-01", "3.65846e-02", +1.0),
}


def _read_gro(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[1])
    return lines[0], lines[2 : 2 + n], lines[2 + n]


def _gro_resname(line: str) -> str:
    return line[5:10].strip()


def _gro_atomname(line: str) -> str:
    return line[10:15].strip()


def _gro_xyz(line: str) -> tuple[float, float, float]:
    return float(line[20:28]), float(line[28:36]), float(line[36:44])


def _renumber_gro(atoms: list[str]) -> list[str]:
    """Rewrite the 5-col atom-serial field (cols 15:20) sequentially from 1."""
    out = []
    for i, line in enumerate(atoms, start=1):
        out.append(line[:15] + f"{i % 100000:>5}" + line[20:])
    return out


def _pick_bulk_ion(atoms: list[str], ion_resname: str) -> int:
    """Index (into `atoms`) of the ion of `ion_resname` most distant from protein."""
    protein_xyz, ion_idx, ion_xyz = [], [], []
    for i, line in enumerate(atoms):
        rn = _gro_resname(line)
        if rn == ion_resname and _gro_atomname(line) == ion_resname:
            ion_idx.append(i)
            ion_xyz.append(_gro_xyz(line))
        elif rn not in _NON_PROTEIN:
            protein_xyz.append(_gro_xyz(line))
    if not ion_idx:
        raise ValueError(f"No {ion_resname} ions found in coordinates")
    if not protein_xyz:
        raise ValueError("No protein atoms found in coordinates")
    P = np.asarray(protein_xyz)
    I = np.asarray(ion_xyz)
    # min distance from each ion to any protein atom, then take the farthest ion.
    d2 = ((I[:, None, :] - P[None, :, :]) ** 2).sum(-1)
    min_d = np.sqrt(d2.min(axis=1))
    return ion_idx[int(min_d.argmax())]


def _split_molecules_line(top_lines: list[str], ion_resname: str, coalch_name: str) -> list[str]:
    """In [ molecules ], decrement `ion_resname` by 1 and append `coalch_name 1`."""
    out, in_mols, done = [], False, False
    for line in top_lines:
        s = line.strip()
        if s.startswith("[") and "molecules" in s.lower():
            in_mols = True
            out.append(line)
            continue
        if in_mols and not done:
            parts = s.split()
            if len(parts) == 2 and parts[0] == ion_resname:
                count = int(parts[1])
                if count < 1:
                    raise ValueError(f"{ion_resname} count < 1 in [ molecules ]")
                out.append(f"{ion_resname:<15} {count - 1}")
                out.append(f"{coalch_name:<15} 1")
                done = True
                continue
        out.append(line)
    if not done:
        raise ValueError(f"Could not find '{ion_resname}' in [ molecules ]")
    return out


def _insert_after_forcefield_include(top_lines: list[str], block: list[str]) -> list[str]:
    """Insert `block` right after the first #include of a forcefield.itp."""
    for i, line in enumerate(top_lines):
        if line.strip().startswith("#include") and "forcefield.itp" in line:
            return top_lines[: i + 1] + [""] + block + top_lines[i + 1 :]
    # Fallback: prepend (still before any [ moleculetype ]).
    return block + [""] + top_lines


def _insert_before_system(top_lines: list[str], block: list[str]) -> list[str]:
    """Insert `block` immediately before the [ system ] section."""
    for i, line in enumerate(top_lines):
        if line.strip().lower().startswith("[ system"):
            return top_lines[:i] + block + [""] + top_lines[i:]
    raise ValueError("No [ system ] section found in topology")


def add_coalchemical_ion(top_path: Path, gro_path: Path, *, delta_q: int) -> dict:
    """Convert one bulk counter-ion into a dual-state (real->dummy) co-alchemical ion.

    delta_q: the protein's net-charge change (A->B). For delta_q = -1 we decouple
    a Cl- (compensation +1); for delta_q = +1 we decouple a Na+ (compensation -1).
    Only |delta_q| == 1 is supported here.
    """
    if delta_q not in (-1, +1):
        raise ValueError(f"Only delta_q == +/-1 supported, got {delta_q}")
    ion_resname = "CL" if delta_q == -1 else "NA"
    coalch_name = f"{ion_resname}_coalch"
    atomtype, mass, sigma, eps, chg = _ION_SPEC[ion_resname]
    dummy_type = f"DUM_{ion_resname}"

    title, atoms, box = _read_gro(gro_path)
    top_lines = top_path.read_text().splitlines()
    if coalch_name in "\n".join(top_lines):
        return {"status": "already-present", "coalch_name": coalch_name}

    # 1) pick the most bulk-solvated ion and move it to the end of the coord list.
    idx = _pick_bulk_ion(atoms, ion_resname)
    chosen = atoms[idx]
    reordered = atoms[:idx] + atoms[idx + 1 :] + [chosen]
    reordered = _renumber_gro(reordered)
    gro_path.write_text("\n".join([title, str(len(reordered)), *reordered, box]) + "\n")

    # 2) dummy atomtype (zero LJ, zero charge) + dual-state moleculetype.
    atomtypes_block = [
        "[ atomtypes ]",
        f"; co-alchemical dummy (vanished {ion_resname} at lambda=1)",
        f"{dummy_type}   0   {mass}   0.00000   A   0.00000e+00   0.00000e+00",
    ]
    moltype_block = [
        "[ moleculetype ]",
        f"; co-alchemical ion: {ion_resname} (state A) -> neutral dummy (state B)",
        f"{coalch_name}    1",
        "[ atoms ]",
        ";  nr   type  resnr  residue  atom  cgnr    charge     mass    typeB    chargeB   massB",
        f"    1   {atomtype:<4}  1      {ion_resname:<3}      {ion_resname:<3}   1"
        f"    {chg:>9.5f}  {mass}   {dummy_type}   0.00000   {mass}",
    ]

    top_lines = _insert_after_forcefield_include(top_lines, atomtypes_block)
    top_lines = _insert_before_system(top_lines, moltype_block)
    top_lines = _split_molecules_line(top_lines, ion_resname, coalch_name)
    top_path.write_text("\n".join(top_lines) + "\n")

    ix, iy, iz = _gro_xyz(chosen)
    return {
        "status": "ok",
        "ion_resname": ion_resname,
        "coalch_name": coalch_name,
        "delta_q": delta_q,
        "chosen_ion_gro_index": idx,
        "chosen_ion_xyz_nm": [ix, iy, iz],
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Insert a co-alchemical ion into a neutral GROMACS system.")
    p.add_argument("--top", type=Path, required=True)
    p.add_argument("--gro", type=Path, required=True)
    p.add_argument("--delta-q", type=int, required=True, help="Protein net-charge change A->B (+/-1)")
    args = p.parse_args(argv)
    try:
        info = add_coalchemical_ion(args.top, args.gro, delta_q=args.delta_q)
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(info)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
