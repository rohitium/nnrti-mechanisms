from __future__ import annotations

from pathlib import Path
from typing import Iterable


def _require(module_name: str):
    try:
        import importlib

        return importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover - dependency gate
        raise RuntimeError(
            f"Missing dependency '{module_name}'. Install required packages and retry."
        ) from exc


def _three_letter(res_name: str) -> str:
    res_name = res_name.strip().upper()
    if len(res_name) == 3:
        return res_name
    one_to_three = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR",
    }
    if res_name not in one_to_three:
        raise ValueError(f"Unknown residue code: {res_name}")
    return one_to_three[res_name]


def _mutation_strings(old_res: str, res_id: str, new_res: str) -> Iterable[str]:
    old_res = _three_letter(old_res)
    new_res = _three_letter(new_res)
    return (
        f"{old_res}-{res_id}-{new_res}",
        f"{old_res}{res_id}{new_res}",
    )


def _residue_name_in_chain(fixer, chain_id: str, residue_id: str) -> str:
    for chain in fixer.topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            if residue.id == residue_id:
                return residue.name
    raise ValueError(f"Residue {residue_id} not found in chain {chain_id}.")


def apply_mutation(
    cif_path: Path,
    chain_id: str,
    residue_id: str,
    new_residue: str,
    output_path: Path,
) -> Path:
    pdbfixer = _require("pdbfixer")
    app = _require("openmm.app")

    with open(cif_path, "r") as handle:
        fixer_ref = pdbfixer.PDBFixer(pdbxfile=handle)
    for mutation_str in _mutation_strings(
        old_res=_residue_name_in_chain(
            fixer_ref,
            chain_id=chain_id,
            residue_id=residue_id,
        ),
        res_id=residue_id,
        new_res=new_residue,
    ):
        with open(cif_path, "r") as handle:
            fixer = pdbfixer.PDBFixer(pdbxfile=handle)
        try:
            fixer.applyMutations([mutation_str], chain_id)
        except Exception:
            continue
        # Verify mutation applied.
        mutated_name = _residue_name_in_chain(fixer, chain_id, residue_id)
        if _three_letter(mutated_name) == _three_letter(new_residue):
            with open(output_path, "w") as handle:
                app.PDBxFile.writeFile(fixer.topology, fixer.positions, handle)
            return output_path

    raise RuntimeError(
        f"Failed to apply mutation for {chain_id}:{residue_id} -> {new_residue}."
    )


def apply_mutations(
    cif_path: Path,
    mutations: Iterable[tuple[str, str, str]],
    output_path: Path,
) -> Path:
    mutations = list(mutations)
    if not mutations:
        raise ValueError("No mutations provided.")
    if len(mutations) == 1:
        chain_id, residue_id, new_residue = mutations[0]
        return apply_mutation(
            cif_path=cif_path,
            chain_id=chain_id,
            residue_id=residue_id,
            new_residue=new_residue,
            output_path=output_path,
        )

    current_path = cif_path
    for idx, (chain_id, residue_id, new_residue) in enumerate(mutations):
        is_last = idx == len(mutations) - 1
        step_path = (
            output_path
            if is_last
            else output_path.with_name(f"{output_path.stem}_step{idx + 1}.cif")
        )
        apply_mutation(
            cif_path=current_path,
            chain_id=chain_id,
            residue_id=residue_id,
            new_residue=new_residue,
            output_path=step_path,
        )
        current_path = step_path
    return output_path
