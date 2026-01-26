from __future__ import annotations

import argparse
from pathlib import Path

from ..config import dor_spec, rpv_spec
from .atoms import ligand_atoms
from .block import load_block
from .build import add_explicit_hydrogens, build_rdkit_mol, write_sdf_file
from .comp import chem_comp_atoms, chem_comp_bonds
from ..utils import ensure_dirs, project_paths


def generate_ligand_sdf(
    cif_path: Path, comp_id: str, chain_id: str, out_path: Path
) -> Path:
    block = load_block(cif_path)
    atom_elements = chem_comp_atoms(block, comp_id)
    bonds = chem_comp_bonds(block, comp_id)
    atoms = ligand_atoms(block, comp_id, chain_id, atom_elements)
    atom_names = {atom.name for atom in atoms}
    bonds = [b for b in bonds if b.atom1 in atom_names and b.atom2 in atom_names]
    mol = build_rdkit_mol(atoms, bonds)
    mol = add_explicit_hydrogens(mol)
    write_sdf_file(mol, out_path)
    return out_path


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
