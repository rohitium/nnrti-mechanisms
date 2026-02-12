from __future__ import annotations

from pathlib import Path

from ...utils import load_residue_mappings


def detect_numbering_scheme(cif_path: Path, chain_map: dict[str, str]) -> dict[str, str]:
    import pdbfixer

    residue_maps = load_residue_mappings(cif_path)
    numbering_scheme: dict[str, str] = {}
    with open(cif_path, "r") as handle:
        fixer = pdbfixer.PDBFixer(pdbxfile=handle)
    for chain in fixer.topology.chains():
        chain_id = chain.id
        if chain_id not in chain_map:
            continue
        maps = residue_maps.get(chain_id, {})
        auth_ids = set(maps.get("auth_map", {}).keys())
        label_ids = set(maps.get("label_map", {}).keys())
        auth_hits = 0
        label_hits = 0
        for res in chain.residues():
            res_id = res.id
            res_name = res.name
            if res_id in auth_ids and maps["auth_map"].get(res_id) == res_name:
                auth_hits += 1
            if res_id in label_ids and maps["label_map"].get(res_id) == res_name:
                label_hits += 1
        if auth_hits == 0 and label_hits == 0:
            raise ValueError(
                f"Unable to determine numbering for chain {chain_id} in {cif_path}"
            )
        numbering_scheme[chain_id] = "label" if label_hits > auth_hits else "auth"
    return numbering_scheme
