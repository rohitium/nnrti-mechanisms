from __future__ import annotations

from typing import Iterable


def require_module(module_name: str):
    try:
        import importlib

        return importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            f"Missing dependency '{module_name}'. Install required packages and retry."
        ) from exc


def three_letter(res_name: str) -> str:
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


def mutation_strings(old_res: str, res_id: str, new_res: str) -> Iterable[str]:
    old_res = three_letter(old_res)
    new_res = three_letter(new_res)
    return (
        f"{old_res}-{res_id}-{new_res}",
        f"{old_res}{res_id}{new_res}",
    )


def residue_name_in_chain(fixer, chain_id: str, residue_id: str) -> str:
    for chain in fixer.topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            if residue.id == residue_id:
                return residue.name
    raise ValueError(f"Residue {residue_id} not found in chain {chain_id}.")
