from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .config import dor_spec, rpv_spec
from .utils import ensure_dirs, project_paths


def _require(module_name: str):
    try:
        import importlib

        return importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover - dependency gate
        raise RuntimeError(
            f"Missing dependency '{module_name}'. Install required packages and retry."
        ) from exc


@dataclass(frozen=True)
class LigandAtom:
    name: str
    element: str
    x: float
    y: float
    z: float
    aromatic: bool


@dataclass(frozen=True)
class LigandBond:
    atom1: str
    atom2: str
    order: int
    aromatic: bool


def _bond_order(value_order: str) -> int:
    mapping = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 4}
    return mapping.get(value_order.upper(), 1)


def _normalize_element(element: str) -> str:
    element = element.strip()
    if not element:
        return element
    return element[0].upper() + element[1:].lower()


def _load_block(cif_path: Path):
    gemmi = _require("gemmi")
    doc = gemmi.cif.read_file(str(cif_path))
    return doc.sole_block()


def _chem_comp_atoms(block, comp_id: str) -> dict[str, tuple[str, bool]]:
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
        atoms[atom_id] = (_normalize_element(element), aromatic_flag.upper() == "Y")
    if not atoms:
        raise ValueError(f"No chem_comp atoms found for {comp_id}.")
    return atoms


def _chem_comp_bonds(block, comp_id: str) -> list[LigandBond]:
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
                order=_bond_order(order),
                aromatic=aromatic_flag.upper() == "Y",
            )
        )
    return bonds


def _atom_site_category(block) -> dict:
    cat = block.get_mmcif_category("_atom_site.")
    if not cat:
        raise ValueError("Category _atom_site not found in CIF.")
    return cat


def _ligand_atoms(
    block, comp_id: str, chain_id: str, atom_elements: dict[str, tuple[str, bool]]
) -> list[LigandAtom]:
    cat = _atom_site_category(block)
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


def _build_rdkit_mol(atoms: list[LigandAtom], bonds: list[LigandBond]):
    _require("rdkit")
    from rdkit import Chem

    mol = Chem.RWMol()
    atom_index = {}
    for atom in atoms:
        rd_atom = Chem.Atom(atom.element)
        rd_atom.SetIsAromatic(atom.aromatic)
        idx = mol.AddAtom(rd_atom)
        atom_index[atom.name] = idx

    for bond in bonds:
        bond_type = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.AROMATIC,
        }[bond.order]
        mol.AddBond(atom_index[bond.atom1], atom_index[bond.atom2], bond_type)
        if bond.aromatic:
            bond_obj = mol.GetBondBetweenAtoms(
                atom_index[bond.atom1], atom_index[bond.atom2]
            )
            bond_obj.SetIsAromatic(True)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom in atoms:
        idx = atom_index[atom.name]
        conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
    mol.AddConformer(conf, assignId=True)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def generate_ligand_sdf(cif_path: Path, comp_id: str, chain_id: str, out_path: Path) -> Path:
    _require("rdkit")
    block = _load_block(cif_path)
    atom_elements = _chem_comp_atoms(block, comp_id)
    bonds = _chem_comp_bonds(block, comp_id)
    atoms = _ligand_atoms(block, comp_id, chain_id, atom_elements)
    atom_names = {atom.name for atom in atoms}
    bonds = [b for b in bonds if b.atom1 in atom_names and b.atom2 in atom_names]
    mol = _build_rdkit_mol(atoms, bonds)
    mol = _add_explicit_hydrogens(mol)
    _write_sdf_file(mol, out_path)
    return out_path


def _add_explicit_hydrogens(mol):
    from rdkit import Chem

    mol_h = Chem.AddHs(mol, addCoords=True)
    if mol_h is None:
        raise RuntimeError("RDKit failed to add hydrogens.")
    return mol_h


def _write_sdf_file(mol, out_path: Path) -> None:
    from rdkit import Chem

    writer = Chem.SDWriter(str(out_path))
    if writer is None:
        raise RuntimeError(f"Failed to open SDF writer for {out_path}")
    writer.write(mol)
    writer.flush()
    writer.close()
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError(f"Failed to write hydrogenated SDF: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate ligand SDF from CIF metadata.")
    parser.add_argument("--cif", type=Path, help="Path to CIF file.")
    parser.add_argument("--resname", type=str, help="Ligand residue name (comp_id).")
    parser.add_argument("--chain", type=str, help="Ligand chain id.")
    parser.add_argument("--out", type=Path, help="Output SDF path.")
    args = parser.parse_args()

    if args.cif and args.resname and args.chain and args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        generate_ligand_sdf(args.cif, args.resname, args.chain, args.out)
        return

    root = Path(__file__).resolve().parents[1]
    paths = project_paths(root)
    ensure_dirs([paths.ligands])
    for spec in [rpv_spec(root), dor_spec(root)]:
        generate_ligand_sdf(
            cif_path=spec.structure.cif_path,
            comp_id=spec.structure.ligand_resname,
            chain_id=spec.structure.ligand_chain,
            out_path=spec.structure.ligand_sdf,
        )


if __name__ == "__main__":
    main()
