from __future__ import annotations

from pathlib import Path

from .helpers import require_module, three_letter


def verify_mutations(
    base_cif_path: Path,
    mut_cif_path: Path,
    targets: list[tuple[str, str, str]],
) -> None:
    pdbfixer = require_module("pdbfixer")
    with open(base_cif_path, "r") as handle:
        base = pdbfixer.PDBFixer(pdbxfile=handle)
    with open(mut_cif_path, "r") as handle:
        mutated = pdbfixer.PDBFixer(pdbxfile=handle)

    by_chain = {}
    for chain_id, old_res, new_res in targets:
        by_chain.setdefault(chain_id, []).append((old_res, new_res))

    for chain_id, expected_pairs in by_chain.items():
        base_res = [res.name for res in _chain_residues(base, chain_id)]
        mut_res = [res.name for res in _chain_residues(mutated, chain_id)]
        if len(base_res) != len(mut_res):
            raise RuntimeError(
                f"Verification failed for chain {chain_id}: residue count changed "
                f"({len(base_res)} -> {len(mut_res)})"
            )
        diffs = [
            (i, base_res[i], mut_res[i])
            for i in range(len(base_res))
            if base_res[i] != mut_res[i]
        ]
        if len(diffs) != len(expected_pairs):
            raise RuntimeError(
                f"Verification failed for chain {chain_id}: expected {len(expected_pairs)} "
                f"changes, found {len(diffs)}"
            )
        expected_set = {(three_letter(o), three_letter(n)) for o, n in expected_pairs}
        found_set = {(three_letter(o), three_letter(n)) for _, o, n in diffs}
        if expected_set != found_set:
            raise RuntimeError(
                f"Verification failed for chain {chain_id}: expected {expected_set}, "
                f"found {found_set}"
            )


def _chain_residues(fixer, chain_id: str):
    for chain in fixer.topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            yield residue
