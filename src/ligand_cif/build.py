from __future__ import annotations

from pathlib import Path

from .types import LigandAtom, LigandBond


def build_rdkit_mol(atoms: list[LigandAtom], bonds: list[LigandBond]):
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


def add_explicit_hydrogens(mol):
    from rdkit import Chem

    mol_h = Chem.AddHs(mol, addCoords=True)
    if mol_h is None:
        raise RuntimeError("RDKit failed to add hydrogens.")
    return mol_h


def write_sdf_file(mol, out_path: Path) -> None:
    from rdkit import Chem

    writer = Chem.SDWriter(str(out_path))
    if writer is None:
        raise RuntimeError(f"Failed to open SDF writer for {out_path}")
    writer.write(mol)
    writer.flush()
    writer.close()
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError(f"Failed to write hydrogenated SDF: {out_path}")
