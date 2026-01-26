from __future__ import annotations

from dataclasses import dataclass

from .platform import get_platform
from .require import require_module


@dataclass(frozen=True)
class EnergyResult:
    complex_kj_mol: float
    receptor_kj_mol: float
    ligand_kj_mol: float
    binding_proxy_kj_mol: float


def _modeller_without(topology, positions, residues_to_delete):
    app = require_module("openmm.app")
    modeller = app.Modeller(topology, positions)
    modeller.delete(residues_to_delete)
    return modeller


def _compute_energy(topology, positions, forcefield) -> float:
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()

    system = forcefield.createSystem(
        topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds
    )
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def compute_binding_proxy(topology, positions, forcefield, ligand_resname: str) -> EnergyResult:
    ligand_residues = [res for res in topology.residues() if res.name == ligand_resname]
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
