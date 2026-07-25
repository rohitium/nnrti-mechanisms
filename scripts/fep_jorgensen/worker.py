"""Minimal Sherlock worker: standard library + OpenMM only."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def _perses_default_parameters(value: float) -> dict[str, float]:
    """Evaluate Perses' default LambdaProtocol without importing Perses."""
    first_half = min(2.0 * value, 1.0)
    second_half = max(2.0 * (value - 0.5), 0.0)
    return {
        "lambda_sterics_core": value,
        "lambda_electrostatics_core": value,
        "lambda_sterics_insert": first_half,
        "lambda_sterics_delete": second_half,
        "lambda_electrostatics_insert": second_half,
        "lambda_electrostatics_delete": first_half,
        "lambda_bonds": value,
        "lambda_angles": value,
        "lambda_torsions": value,
    }


def _find_nonbonded_force(system):
    from openmm import NonbondedForce

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            return force
    raise ValueError("System does not contain an OpenMM NonbondedForce.")


def _capture_nonbonded_parameters(nonbonded):
    particles = [nonbonded.getParticleParameters(i) for i in range(nonbonded.getNumParticles())]
    exceptions = [nonbonded.getExceptionParameters(i) for i in range(nonbonded.getNumExceptions())]
    return particles, exceptions


def _set_nonbonded_strength(nonbonded, original, alchemical_atoms: set[int], strength: float) -> None:
    particles, exceptions = original
    charge_scale = math.sqrt(strength)
    lj_scale = strength
    for atom_index, (charge, sigma, epsilon) in enumerate(particles):
        if atom_index in alchemical_atoms:
            nonbonded.setParticleParameters(atom_index, charge * charge_scale, sigma, epsilon * lj_scale)
        else:
            nonbonded.setParticleParameters(atom_index, charge, sigma, epsilon)
    for exception_index, (i, j, charge_prod, sigma, epsilon) in enumerate(exceptions):
        if int(i) in alchemical_atoms or int(j) in alchemical_atoms:
            nonbonded.setExceptionParameters(
                exception_index,
                i,
                j,
                charge_prod * strength,
                sigma,
                epsilon * lj_scale,
            )
        else:
            nonbonded.setExceptionParameters(exception_index, i, j, charge_prod, sigma, epsilon)


def _set_lambda(
    context,
    schedule: dict,
    value: float,
    *,
    original_nonbonded,
) -> None:
    protocol = schedule.get("lambda_parameter_functions", "nonbonded-scaling")
    if protocol == "perses-default":
        available = set(context.getParameters())
        for name, parameter_value in _perses_default_parameters(value).items():
            if name in available:
                context.setParameter(name, parameter_value)
        return
    if protocol != "nonbonded-scaling":
        raise ValueError(f"Unsupported lambda parameter functions: {protocol}")
    nonbonded = _find_nonbonded_force(context.getSystem())
    alchemical_atoms = set(schedule["alchemical_plan"]["alchemical_atom_indices"])
    _set_nonbonded_strength(nonbonded, original_nonbonded, alchemical_atoms, float(value))
    nonbonded.updateParametersInContext(context)


def _completed_samples(csv_path: Path) -> int:
    if not csv_path.exists():
        return 0
    with csv_path.open(newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _truncate_samples(csv_path: Path, keep: int) -> None:
    with csv_path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    if rows:
        with csv_path.open("w", newline="") as handle:
            csv.writer(handle).writerows(rows[: keep + 1])


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
    protocol = schedule.get("lambda_parameter_functions", "nonbonded-scaling")
    if protocol not in {"perses-default", "nonbonded-scaling"}:
        raise ValueError(f"Unsupported lambda parameter functions: {protocol}")
    if not 0 <= state_index < len(lambdas):
        raise IndexError(f"state-index {state_index} outside 0..{len(lambdas)-1}")
    if production_steps % energy_interval:
        raise ValueError("production_steps must be divisible by energy_interval")
    if checkpoint_interval % energy_interval:
        raise ValueError("checkpoint_interval must be divisible by energy_interval")
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
    original_nonbonded = None
    if protocol == "nonbonded-scaling":
        original_nonbonded = _capture_nonbonded_parameters(_find_nonbonded_force(system))
    _set_lambda(
        simulation.context,
        schedule,
        lambdas[state_index],
        original_nonbonded=original_nonbonded,
    )

    target_samples = production_steps // energy_interval
    existing_samples = _completed_samples(csv_path)
    if checkpoint.exists() and existing_samples > target_samples:
        _truncate_samples(csv_path, target_samples)
        return csv_path
    if checkpoint.exists() and existing_samples == target_samples:
        return csv_path
    if checkpoint.exists():
        samples_per_checkpoint = checkpoint_interval // energy_interval
        resumable_samples = (
            existing_samples // samples_per_checkpoint
        ) * samples_per_checkpoint
        if resumable_samples == 0:
            raise RuntimeError(
                f"Checkpoint {checkpoint} exists without a complete checkpoint block in {csv_path}"
            )
        if resumable_samples != existing_samples:
            _truncate_samples(csv_path, resumable_samples)
        completed_samples = resumable_samples
        simulation.loadCheckpoint(str(checkpoint))
    else:
        completed_samples = 0
        if csv_path.exists():
            csv_path.unlink()
        simulation.minimizeEnergy(maxIterations=500)
        simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
        simulation.step(equilibration_steps)

    beta_per_kj_mol = 1.0 / (
        unit.MOLAR_GAS_CONSTANT_R * temperature_k * unit.kelvin
    ).value_in_unit(unit.kilojoule_per_mole)
    write_header = not csv_path.exists()
    with csv_path.open("a", newline="") as handle:
        writer = csv.writer(handle)
        if write_header:
            writer.writerow(["sample", "origin_state", *[f"u_{i}" for i in range(len(lambdas))]])
        for sample in range(completed_samples, target_samples):
            simulation.step(energy_interval)
            reduced = []
            for value in lambdas:
                _set_lambda(
                    simulation.context,
                    schedule,
                    value,
                    original_nonbonded=original_nonbonded,
                )
                energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
                energy_kj_mol = energy.value_in_unit(unit.kilojoule_per_mole)
                reduced.append(float(beta_per_kj_mol * energy_kj_mol))
            _set_lambda(
                simulation.context,
                schedule,
                lambdas[state_index],
                original_nonbonded=original_nonbonded,
            )
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
    run_window(
        phase_dir=args.phase_dir,
        output_dir=args.output_dir,
        state_index=args.state_index,
        temperature_k=args.temperature_k,
        timestep_fs=args.timestep_fs,
        collision_rate_per_ps=args.collision_rate_per_ps,
        equilibration_steps=args.equilibration_steps,
        production_steps=args.production_steps,
        energy_interval=args.energy_interval,
        checkpoint_interval=args.checkpoint_interval,
        platform_name=args.platform,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
