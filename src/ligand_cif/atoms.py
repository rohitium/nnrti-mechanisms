from __future__ import annotations

from .types import LigandAtom


def atom_site_category(block) -> dict:
    cat = block.get_mmcif_category("_atom_site.")
    if not cat:
        raise ValueError("Category _atom_site not found in CIF.")
    return cat


def ligand_atoms(
    block, comp_id: str, chain_id: str, atom_elements: dict[str, tuple[str, bool]]
) -> list[LigandAtom]:
    cat = atom_site_category(block)
    keys = set(cat.keys())
    comp_id_key = "label_comp_id" if "label_comp_id" in keys else "auth_comp_id"
    asym_id_key = "label_asym_id" if "label_asym_id" in keys else "auth_asym_id"
    atom_id_key = "label_atom_id" if "label_atom_id" in keys else "auth_atom_id"
    atoms = []
    for idx in range(len(cat[comp_id_key])):
        if cat[comp_id_key][idx] != comp_id:
            continue
        if cat[asym_id_key][idx] != chain_id:
            continue
        atom_name = cat[atom_id_key][idx]
        atom_info = atom_elements.get(atom_name)
        if atom_info is None:
            continue
        element, aromatic = atom_info
        atoms.append(
            LigandAtom(
                name=atom_name,
                element=element,
                x=float(cat["Cartn_x"][idx]),
                y=float(cat["Cartn_y"][idx]),
                z=float(cat["Cartn_z"][idx]),
                aromatic=aromatic,
            )
        )
    if not atoms:
        raise ValueError(f"No ligand atoms found for {comp_id} chain {chain_id}.")
    return atoms
