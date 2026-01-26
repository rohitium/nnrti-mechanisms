from __future__ import annotations

from pathlib import Path

from .ligand import build_forcefield, load_ligand_molecule
from .minimizer import minimize_system
from .require import require_module
from .restraints import heavy_atom_indices, restrained_indices


def minimize_with_restraints(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius_angstrom: float,
    restraint_k_kj_mol_nm2: float,
    output_path: Path,
) -> tuple:
    app = require_module("openmm.app")
    pdbfixer = require_module("pdbfixer")

    if not ligand_sdf.exists():
        raise FileNotFoundError(
            f"Ligand SDF not found: {ligand_sdf}. Provide a matching SDF."
        )

    ligand_mol = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand_mol])
    with open(cif_path, "r") as handle:
        fixer = pdbfixer.PDBFixer(pdbxfile=handle)
    for residue in fixer.topology.residues():
        if residue.name == "OMC":
            residue.name = "DC"
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)

    omc_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if atom.residue.name == "DC" and atom.name in {"O2'", "CM2"}
    ]
    if omc_atoms:
        modeller.delete(omc_atoms)
    hydrogens = [
        atom
        for atom in modeller.topology.atoms()
        if atom.element and atom.element.symbol == "H" and atom.residue.name != ligand_resname
    ]
    if hydrogens:
        modeller.delete(hydrogens)

    ligand_residues = [res for res in modeller.topology.residues() if res.name == ligand_resname]
    if ligand_residues:
        modeller.delete(ligand_residues)
        off_unit = require_module("openff.units").unit
        omm_unit = require_module("openmm.unit")
        ligand_topology = ligand_mol.to_topology().to_openmm()
        for residue in ligand_topology.residues():
            residue.name = ligand_resname
        ligand_positions = ligand_mol.conformers[0].to(off_unit.nanometer).magnitude
        modeller.add(ligand_topology, ligand_positions * omm_unit.nanometer)
    modeller.addHydrogens(forcefield)

    ligand_indices = [
        atom.index for atom in modeller.topology.atoms() if atom.residue.name == ligand_resname
    ]
    heavy_indices = heavy_atom_indices(modeller.topology, ligand_resname)
    restraint_indices = restrained_indices(
        modeller.positions, ligand_indices, heavy_indices, restraint_radius_angstrom
    )

    _, positions = minimize_system(
        modeller.topology,
        modeller.positions,
        forcefield,
        restraint_indices,
        restraint_k_kj_mol_nm2,
    )

    with open(output_path, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, positions, handle)
    return modeller.topology, positions, forcefield
