from __future__ import annotations

import logging
import time
from typing import Sequence

from .platform import get_platform
from .require import require_module


def minimize_system(
    topology,
    positions,
    forcefield,
    restraint_indices: Sequence[int],
    restraint_k_kj_mol_nm2: float,
):
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()

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
    return system, state.getPositions()


def run_implicit_md(
    topology,
    positions,
    forcefield,
    restraint_indices: Sequence[int],
    restraint_k_kj_mol_nm2: float,
    temperature_k: float,
    friction_per_ps: float,
    timestep_ps: float,
    steps: int,
    report_interval: int,
    output_dcd,
    output_pdb,
    output_cif,
):
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    # Add an implicit solvent model manually because the SMIRNOFF ligand
    # templates do not provide GBSA parameters for ForceField.createSystem().
    gbsa = openmm.GBSAOBCForce()
    gbsa.setSolventDielectric(78.5)
    gbsa.setSoluteDielectric(1.0)
    element_radii_nm = {
        "H": 0.12,
        "C": 0.17,
        "N": 0.155,
        "O": 0.152,
        "F": 0.147,
        "P": 0.18,
        "S": 0.18,
        "CL": 0.175,
        "BR": 0.185,
        "I": 0.198,
    }
    nonbonded = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nonbonded = force
            break
    if nonbonded is None:
        raise RuntimeError("NonbondedForce not found; cannot assign GBSA parameters.")

    for atom in topology.atoms():
        charge, sigma, _ = nonbonded.getParticleParameters(atom.index)
        symbol = atom.element.symbol.upper() if atom.element else ""
        radius = element_radii_nm.get(symbol, 0.17)
        gbsa.addParticle(charge, radius, 1.0)
    system.addForce(gbsa)
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
        temperature_k * unit.kelvin,
        friction_per_ps / unit.picosecond,
        timestep_ps * unit.picoseconds,
    )
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
    if output_dcd is not None:
        simulation.reporters.append(app.DCDReporter(str(output_dcd), report_interval))

    start = time.perf_counter()
    simulation.step(steps)
    logging.info("Implicit MD completed in %.2fs", time.perf_counter() - start)

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    md_positions = state.getPositions()
    if output_pdb is not None:
        with open(output_pdb, "w") as handle:
            app.PDBFile.writeFile(topology, md_positions, handle)
    if output_cif is not None:
        with open(output_cif, "w") as handle:
            app.PDBxFile.writeFile(topology, md_positions, handle)
    return md_positions
