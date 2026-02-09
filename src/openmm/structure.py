from __future__ import annotations

from pathlib import Path

import numpy as np

from .ligand import build_forcefield, load_ligand_molecule
from .minimizer import minimize_system
from .require import require_module
from .restraints import heavy_atom_indices, restrained_indices


def _jitter_positions(positions, seed: int | None, jitter_angstrom: float):
    if jitter_angstrom <= 0.0:
        return positions
    unit = require_module("openmm.unit")
    openmm = require_module("openmm")
    rng = np.random.default_rng(seed)
    pos_nm = np.array([p.value_in_unit(unit.nanometer) for p in positions])
    noise = rng.normal(scale=jitter_angstrom / 10.0, size=pos_nm.shape)
    pos_nm = pos_nm + noise
    return [openmm.Vec3(*xyz) * unit.nanometer for xyz in pos_nm]


def minimize_with_restraints(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius_angstrom: float,
    restraint_k_kj_mol_nm2: float,
    output_path: Path,
    jitter_seed: int | None = None,
    jitter_angstrom: float = 0.0,
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
    # Keep missing-residue handling consistent across WT and mutated CIFs.
    # Mutated CIFs lose missing-residue annotations, so we avoid adding them here too.
    fixer.missingResidues = {}
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

    # Rebuild ligand from the SDF template so ForceField templates match, but preserve
    # the bound pose by rigidly aligning template coordinates to the input CIF ligand.
    ligand_atoms = [
        atom for atom in modeller.topology.atoms() if atom.residue.name == ligand_resname
    ]
    if not ligand_atoms:
        raise ValueError(
            f"Ligand residue '{ligand_resname}' not found in {cif_path}. "
            "Cannot prepare bound complex."
        )
    omm_unit = require_module("openmm.unit")
    off_unit = require_module("openff.units").unit
    def _is_hydrogen(atom) -> bool:
        if getattr(atom, "element", None) is not None:
            return atom.element.symbol.upper() == "H"
        return atom.name.upper().startswith("H")

    original_ligand_xyz = np.array(
        [modeller.positions[a.index].value_in_unit(omm_unit.nanometer) for a in ligand_atoms],
        dtype=float,
    )
    original_heavy_xyz = np.array(
        [
            modeller.positions[a.index].value_in_unit(omm_unit.nanometer)
            for a in ligand_atoms
            if not _is_hydrogen(a)
        ],
        dtype=float,
    )
    ligand_residues = [
        res for res in modeller.topology.residues() if res.name == ligand_resname
    ]
    modeller.delete(ligand_residues)
    ligand_topology = ligand_mol.to_topology().to_openmm()
    for residue in ligand_topology.residues():
        residue.name = ligand_resname
    ligand_template_xyz = np.asarray(
        ligand_mol.conformers[0].to(off_unit.nanometer).magnitude,
        dtype=float,
    )
    template_heavy_mask = np.array(
        [a.atomic_number != 1 for a in ligand_mol.atoms],
        dtype=bool,
    )
    template_heavy_xyz = ligand_template_xyz[template_heavy_mask]
    if template_heavy_xyz.shape[0] != original_heavy_xyz.shape[0]:
        raise ValueError(
            "Ligand heavy-atom count mismatch between CIF and SDF; cannot align bound pose "
            f"({original_heavy_xyz.shape[0]} vs {template_heavy_xyz.shape[0]})."
        )

    # Kabsch rigid alignment on heavy atoms: template (mobile) -> CIF ligand (target).
    mobile_centroid = template_heavy_xyz.mean(axis=0)
    target_centroid = original_heavy_xyz.mean(axis=0)
    m = template_heavy_xyz - mobile_centroid
    t = original_heavy_xyz - target_centroid
    cov = m.T @ t
    u, _, vt = np.linalg.svd(cov)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    aligned_xyz = ((ligand_template_xyz - mobile_centroid) @ r) + target_centroid
    modeller.add(ligand_topology, aligned_xyz * omm_unit.nanometer)
    modeller.addHydrogens(forcefield)
    modeller.positions = _jitter_positions(
        modeller.positions, seed=jitter_seed, jitter_angstrom=jitter_angstrom
    )

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
