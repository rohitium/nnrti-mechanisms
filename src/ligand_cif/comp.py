from __future__ import annotations

from .block import bond_order, normalize_element
from .types import LigandBond


def chem_comp_atoms(block, comp_id: str) -> dict[str, tuple[str, bool]]:
    cat = block.get_mmcif_category("_chem_comp_atom.")
    if not cat:
        raise ValueError("Category _chem_comp_atom not found in CIF.")
    rows = zip(
        cat["comp_id"],
        cat["atom_id"],
        cat["type_symbol"],
        cat["pdbx_aromatic_flag"],
    )
    atoms = {}
    for comp, atom_id, element, aromatic_flag in rows:
        if comp != comp_id:
            continue
        atoms[atom_id] = (normalize_element(element), aromatic_flag.upper() == "Y")
    if not atoms:
        raise ValueError(f"No chem_comp atoms found for {comp_id}.")
    return atoms


def chem_comp_bonds(block, comp_id: str) -> list[LigandBond]:
    cat = block.get_mmcif_category("_chem_comp_bond.")
    if not cat:
        raise ValueError("Category _chem_comp_bond not found in CIF.")
    rows = zip(
        cat["comp_id"],
        cat["atom_id_1"],
        cat["atom_id_2"],
        cat["value_order"],
        cat["pdbx_aromatic_flag"],
    )
    bonds = []
    for comp, atom1, atom2, order, aromatic_flag in rows:
        if comp != comp_id:
            continue
        bonds.append(
            LigandBond(
                atom1=atom1,
                atom2=atom2,
                order=bond_order(order),
                aromatic=aromatic_flag.upper() == "Y",
            )
        )
    return bonds
