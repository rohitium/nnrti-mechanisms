from __future__ import annotations

from dataclasses import dataclass
import logging
import os
import time
from pathlib import Path
from typing import Sequence

import numpy as np


def _require(module_name: str):
    try:
        import importlib

        return importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover - dependency gate
        raise RuntimeError(
            f"Missing dependency '{module_name}'. Install required packages and retry."
        ) from exc


@dataclass(frozen=True)
class EnergyResult:
    complex_kj_mol: float
    receptor_kj_mol: float
    ligand_kj_mol: float
    binding_proxy_kj_mol: float


def load_ligand_molecule(ligand_sdf: Path):
    openff = _require("openff.toolkit")
    mols = openff.topology.Molecule.from_file(
        str(ligand_sdf), allow_undefined_stereo=True
    )
    if isinstance(mols, list):
        if not mols:
            raise ValueError(f"No molecules found in {ligand_sdf}")
        molecules = mols
    else:
        molecules = [mols]

    for mol in molecules:
        mol.assign_partial_charges(partial_charge_method="gasteiger")

    return molecules[0]


def build_forcefield(ligand_molecules) -> "openmm.app.ForceField":
    app = _require("openmm.app")
    generators = _require("openmmforcefields.generators")

    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/DNA.bsc1.xml")
    generator = generators.SMIRNOFFTemplateGenerator(
        molecules=ligand_molecules,
        template_generator_kwargs={"partial_charge_method": "gasteiger"},
    )
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield


def _minimize(
    topology,
    positions,
    forcefield,
    restraint_indices: Sequence[int],
    restraint_k_kj_mol_nm2: float,
):
    app = _require("openmm.app")
    openmm = _require("openmm")
    unit = _require("openmm.unit")
    platform, properties = _get_platform()

    system = forcefield.createSystem(
        topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds
    )

    if restraint_indices:
        restraint = openmm.CustomExternalForce("k*(x-x0)^2 + k*(y-y0)^2 + k*(z-z0)^2")
        restraint.addGlobalParameter("k", restraint_k_kj_mol_nm2)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        for idx in restraint_indices:
            pos = positions[idx]
            restraint.addParticle(idx, pos.value_in_unit(unit.nanometer))
        system.addForce(restraint)

    integrator = openmm.LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    start = time.perf_counter()
    simulation.minimizeEnergy()
    logging.info("OpenMM minimization completed in %.2fs", time.perf_counter() - start)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    return system, state.getPositions(), state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )


def _modeller_without(topology, positions, residues_to_delete):
    app = _require("openmm.app")
    modeller = app.Modeller(topology, positions)
    modeller.delete(residues_to_delete)
    return modeller


def _compute_energy(topology, positions, forcefield) -> float:
    app = _require("openmm.app")
    openmm = _require("openmm")
    unit = _require("openmm.unit")
    platform, properties = _get_platform()

    system = forcefield.createSystem(
        topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds
    )
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def _heavy_atom_indices(topology, exclude_resname: str) -> list[int]:
    indices = []
    for atom in topology.atoms():
        if atom.element is None:
            continue
        if atom.residue.name == exclude_resname:
            continue
        if atom.element.symbol != "H":
            indices.append(atom.index)
    return indices


def _restrained_indices(
    positions,
    ligand_indices: Sequence[int],
    candidate_indices: Sequence[int],
    radius_angstrom: float,
) -> list[int]:
    unit = _require("openmm.unit")
    radius_nm = radius_angstrom / 10.0
    ligand_pos = np.array(
        [positions[i].value_in_unit(unit.nanometer) for i in ligand_indices]
    )
    restrained = []
    for idx in candidate_indices:
        pos = positions[idx].value_in_unit(unit.nanometer)
        d = np.min(np.linalg.norm(ligand_pos - pos, axis=1))
        if d > radius_nm:
            restrained.append(idx)
    return restrained


def minimize_with_restraints(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius_angstrom: float,
    restraint_k_kj_mol_nm2: float,
    output_path: Path,
) -> tuple:
    app = _require("openmm.app")
    pdbfixer = _require("pdbfixer")
    logging.info("Loading structure: %s", cif_path)

    if not ligand_sdf.exists():
        raise FileNotFoundError(
            f"Ligand SDF not found: {ligand_sdf}. Provide a matching SDF."
        )

    ligand_mol = load_ligand_molecule(ligand_sdf)
    logging.info("Ligand loaded and charged from %s", ligand_sdf)
    forcefield = build_forcefield([ligand_mol])
    logging.info("Forcefield initialized")
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
    omc_atom_names = {"O2'", "CM2"}
    omc_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if atom.residue.name == "DC" and atom.name in omc_atom_names
    ]
    if omc_atoms:
        modeller.delete(omc_atoms)
    hydrogens = [
        atom
        for atom in modeller.topology.atoms()
        if atom.element
        and atom.element.symbol == "H"
        and atom.residue.name != ligand_resname
    ]
    if hydrogens:
        modeller.delete(hydrogens)

    # Replace ligand with the hydrogenated SDF version to ensure template matching.
    ligand_residues = [res for res in modeller.topology.residues() if res.name == ligand_resname]
    if ligand_residues:
        modeller.delete(ligand_residues)
        off_unit = _require("openff.units").unit
        omm_unit = _require("openmm.unit")
        ligand_topology = ligand_mol.to_topology().to_openmm()
        for residue in ligand_topology.residues():
            residue.name = ligand_resname
        ligand_positions = ligand_mol.conformers[0].to(off_unit.nanometer).magnitude
        modeller.add(ligand_topology, ligand_positions * omm_unit.nanometer)
    modeller.addHydrogens(forcefield)
    logging.info("Hydrogens added for protein/DNA")

    ligand_indices = [
        atom.index
        for atom in modeller.topology.atoms()
        if atom.residue.name == ligand_resname
    ]
    heavy_indices = _heavy_atom_indices(modeller.topology, ligand_resname)
    restraint_indices = _restrained_indices(
        modeller.positions, ligand_indices, heavy_indices, restraint_radius_angstrom
    )
    logging.info("Restraints applied to %d atoms", len(restraint_indices))

    system, positions, _ = _minimize(
        modeller.topology,
        modeller.positions,
        forcefield,
        restraint_indices,
        restraint_k_kj_mol_nm2,
    )

    with open(output_path, "w") as handle:
        app.PDBxFile.writeFile(modeller.topology, positions, handle)
    pdb_path = output_path.with_suffix(".pdb")
    with open(pdb_path, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, positions, handle)
    logging.info("Wrote minimized structures: %s and %s", output_path, pdb_path)
    return modeller.topology, positions, forcefield


def _get_platform():
    openmm = _require("openmm")
    platform_name = os.environ.get("OPENMM_PLATFORM", "").strip()
    if not platform_name:
        platform_name = (
            "Metal"
            if any(
                openmm.Platform.getPlatform(i).getName() == "Metal"
                for i in range(openmm.Platform.getNumPlatforms())
            )
            else "CPU"
        )
    platform = openmm.Platform.getPlatformByName(platform_name)
    properties = {}
    if platform_name == "CPU":
        threads = os.environ.get("OPENMM_CPU_THREADS")
        if threads:
            properties["Threads"] = threads
    if platform_name in {"OpenCL", "Metal", "CUDA"}:
        device_index = os.environ.get("OPENMM_DEVICE_INDEX")
        if device_index:
            properties["DeviceIndex"] = device_index
    logging.info("OpenMM platform: %s (%s)", platform_name, properties or "default")
    return platform, properties


def compute_binding_proxy(
    topology, positions, forcefield, ligand_resname: str
) -> EnergyResult:
    ligand_residues = [
        res for res in topology.residues() if res.name == ligand_resname
    ]
    if not ligand_residues:
        raise ValueError(f"Ligand residue '{ligand_resname}' not found.")

    complex_energy = _compute_energy(topology, positions, forcefield)

    modeller_receptor = _modeller_without(topology, positions, ligand_residues)
    receptor_energy = _compute_energy(
        modeller_receptor.topology, modeller_receptor.positions, forcefield
    )

    modeller_ligand = _modeller_without(
        topology,
        positions,
        [res for res in topology.residues() if res.name != ligand_resname],
    )
    ligand_energy = _compute_energy(
        modeller_ligand.topology, modeller_ligand.positions, forcefield
    )

    return EnergyResult(
        complex_kj_mol=complex_energy,
        receptor_kj_mol=receptor_energy,
        ligand_kj_mol=ligand_energy,
        binding_proxy_kj_mol=complex_energy - receptor_energy - ligand_energy,
    )
