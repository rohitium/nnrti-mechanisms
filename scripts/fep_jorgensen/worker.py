from __future__ import annotations

"""Minimal Sherlock worker: standard library + OpenMM only."""

import argparse
import csv
import json
from pathlib import Path


LAMBDA_PARAMETERS = (
    "lambda_sterics_core",
    "lambda_electrostatics_core",
    "lambda_sterics_insert",
    "lambda_sterics_delete",
    "lambda_electrostatics_insert",
    "lambda_electrostatics_delete",
    "lambda_bonds",
    "lambda_angles",
    "lambda_torsions",
)


def _set_lambda(context, value: float) -> None:
    """Set Perses default linear protocol on a serialized hybrid system."""
    available = set(context.getParameters())
    for name in LAMBDA_PARAMETERS:
        if name in available:
            context.setParameter(name, value)


def run_window(
    phase_dir: Path,
    output_dir: Path,
    state_index: int,
    temperature_k: float,
    timestep_fs: float,
    collision_rate_per_ps: float,
    equilibration_steps: int,
    production_steps: int,
    energy_interval: int,
    checkpoint_interval: int,
    platform_name: str,
) -> Path:
    from openmm import LangevinMiddleIntegrator, Platform, XmlSerializer, unit
    from openmm import app

    schedule = json.loads((phase_dir / "schedule.json").read_text())
    lambdas = [float(x) for x in schedule["lambda_values"]]
    if not 0 <= state_index < len(lambdas):
        raise IndexError(f"state-index {state_index} outside 0..{len(lambdas)-1}")
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = output_dir / f"state_{state_index:02d}.chk"
    csv_path = output_dir / f"state_{state_index:02d}_energies.csv"

    pdb = app.PDBFile(str(phase_dir / "hybrid_topology.pdb"))
    system = XmlSerializer.deserialize((phase_dir / "hybrid_system.xml").read_text())
    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        collision_rate_per_ps / unit.picosecond,
        timestep_fs * unit.femtosecond,
    )
    platform = Platform.getPlatformByName(platform_name)
    properties = {"Precision": "mixed"} if platform_name in {"CUDA", "OpenCL"} else {}
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    _set_lambda(simulation.context, lambdas[state_index])

    if checkpoint.exists():
        simulation.loadCheckpoint(str(checkpoint))
    else:
        simulation.minimizeEnergy(maxIterations=500)
        simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
        simulation.step(equilibration_steps)

    beta = 1.0 / (
        unit.MOLAR_GAS_CONSTANT_R
        * temperature_k
        * unit.kelvin
    )
    write_header = not csv_path.exists()
    with csv_path.open("a", newline="") as handle:
        writer = csv.writer(handle)
        if write_header:
            writer.writerow(["sample", "origin_state", *[f"u_{i}" for i in range(len(lambdas))]])
        samples = production_steps // energy_interval
        for sample in range(samples):
            simulation.step(energy_interval)
            reduced = []
            for value in lambdas:
                _set_lambda(simulation.context, value)
                energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
                reduced.append(float((beta * energy).value_in_unit(unit.dimensionless)))
            _set_lambda(simulation.context, lambdas[state_index])
            writer.writerow([sample, state_index, *reduced])
            handle.flush()
            if (sample + 1) * energy_interval % checkpoint_interval == 0:
                simulation.saveCheckpoint(str(checkpoint))
    simulation.saveCheckpoint(str(checkpoint))
    return csv_path


def main() -> int:
    p = argparse.ArgumentParser(description="Run one OpenMM FEP lambda window")
    p.add_argument("--phase-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--state-index", type=int, required=True)
    p.add_argument("--temperature-k", type=float, default=300.0)
    p.add_argument("--timestep-fs", type=float, default=2.0)
    p.add_argument("--collision-rate-per-ps", type=float, default=1.0)
    p.add_argument("--equilibration-steps", type=int, default=250_000)
    p.add_argument("--production-steps", type=int, default=2_500_000)
    p.add_argument("--energy-interval", type=int, default=2_500)
    p.add_argument("--checkpoint-interval", type=int, default=25_000)
    p.add_argument("--platform", default="CUDA")
    args = p.parse_args()
    run_window(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
