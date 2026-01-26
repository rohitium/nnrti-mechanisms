from __future__ import annotations

from pathlib import Path
from typing import Iterable

from .helpers import (
    mutation_strings,
    require_module,
    residue_name_in_chain,
    three_letter,
)


def apply_mutation(
    cif_path: Path,
    chain_id: str,
    residue_id: str,
    new_residue: str,
    output_path: Path,
) -> Path:
    pdbfixer = require_module("pdbfixer")
    app = require_module("openmm.app")

    with open(cif_path, "r") as handle:
        fixer_ref = pdbfixer.PDBFixer(pdbxfile=handle)
    for mutation_str in mutation_strings(
        old_res=residue_name_in_chain(
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
        mutated_name = residue_name_in_chain(fixer, chain_id, residue_id)
        if three_letter(mutated_name) == three_letter(new_residue):
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

    pdbfixer = require_module("pdbfixer")
    app = require_module("openmm.app")

    with open(cif_path, "r") as handle:
        fixer = pdbfixer.PDBFixer(pdbxfile=handle)

    by_chain = {}
    for chain_id, residue_id, new_residue in mutations:
        old_res = residue_name_in_chain(fixer, chain_id, residue_id)
        by_chain.setdefault(chain_id, []).append(
            f"{three_letter(old_res)}-{residue_id}-{three_letter(new_residue)}"
        )

    for chain_id, mutation_strs in by_chain.items():
        fixer.applyMutations(mutation_strs, chain_id)
    with open(output_path, "w") as handle:
        app.PDBxFile.writeFile(fixer.topology, fixer.positions, handle)
    return output_path
