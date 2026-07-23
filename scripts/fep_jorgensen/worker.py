"""Minimal Sherlock worker: standard library + OpenMM only."""

from __future__ import annotations

import argparse
import csv
import json
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


def _set_lambda(context, value: float) -> None:
    """Set Perses' default staged protocol on a serialized hybrid system."""
    available = set(context.getParameters())
    for name, parameter_value in _perses_default_parameters(value).items():
        if name in available:
            context.setParameter(name, parameter_value)


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
    protocol = schedule.get("lambda_parameter_functions", "perses-default")
    if protocol != "perses-default":
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
    _set_lambda(simulation.context, lambdas[state_index])

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
                _set_lambda(simulation.context, value)
                energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
                energy_kj_mol = energy.value_in_unit(unit.kilojoule_per_mole)
                reduced.append(float(beta_per_kj_mol * energy_kj_mol))
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
