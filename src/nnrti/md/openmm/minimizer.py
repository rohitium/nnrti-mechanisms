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
