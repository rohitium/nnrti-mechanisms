"""Jorgensen-inspired OpenMM equilibration for FEP input structures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .approx_protocol import ApproxJorgensenProtocol
from .config import FEPConfig
from src.md.openmm.ligand import build_forcefield, load_ligand_molecule
from src.md.openmm.platform import get_platform


def equilibrate_complex(
    input_pdb: Path,
    output_pdb: Path,
    ligand_sdf: Path,
    protocol: ApproxJorgensenProtocol | None = None,
    platform_name: str | None = None,
) -> Path:
    """Relax a bound complex before alchemical mutation setup."""
    from openmm import LangevinMiddleIntegrator, Platform, unit
    from openmm import app

    settings = protocol or ApproxJorgensenProtocol()
    settings.validate()
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])
    with input_pdb.open() as handle:
        pdb = app.PDBFile(handle)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    system.addForce(
        app.MonteCarloBarostat(
            1.0 * unit.atmosphere,
            settings.md_equilibration_temperature_k * unit.kelvin,
            25,
        )
    )
    if platform_name:
        platform = Platform.getPlatformByName(platform_name)
        properties = {"Precision": "mixed"} if platform_name in {"CUDA", "OpenCL"} else {}
    else:
        platform, properties = get_platform()
    integrator = LangevinMiddleIntegrator(
        settings.md_initial_temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        settings.md_timestep_fs * unit.femtosecond,
    )
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)

    simulation.minimizeEnergy(maxIterations=settings.md_minimization_iterations)
    simulation.context.setVelocitiesToTemperature(settings.md_initial_temperature_k * unit.kelvin)
    simulation.step(settings.md_initial_equilibration_steps)

    integrator.setTemperature(settings.md_equilibration_temperature_k * unit.kelvin)
    simulation.context.setVelocitiesToTemperature(settings.md_equilibration_temperature_k * unit.kelvin)
    simulation.step(settings.md_equilibration_steps)

    if settings.quench_blocks > 0:
        block_temps = [
            settings.quench_start_k
            + (settings.quench_end_k - settings.quench_start_k)
            * block
            / settings.quench_blocks
            for block in range(1, settings.quench_blocks + 1)
        ]
        for temperature_k in block_temps:
            integrator.setTemperature(temperature_k * unit.kelvin)
            simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
            simulation.step(settings.quench_block_steps)

    integrator.setTemperature(settings.md_equilibration_temperature_k * unit.kelvin)
    simulation.context.setVelocitiesToTemperature(settings.md_equilibration_temperature_k * unit.kelvin)
    with output_pdb.open("w") as handle:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            handle,
            keepIds=True,
        )
    return output_pdb


def equilibrate_leg(config: FEPConfig, platform_name: str | None = None) -> Path:
    config.validate(require_inputs=True)
    settings = config.approx_protocol
    input_pdb = config.wt_complex_pdb
    output_pdb = config.run_dir / "inputs" / "equilibrated_complex.pdb"
    if output_pdb.exists():
        return output_pdb
    equilibrate_complex(
        input_pdb,
        output_pdb,
        config.ligand_sdf,
        settings,
        platform_name=platform_name,
    )
    manifest = {
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb),
        "protocol": settings.to_dict(),
    }
    (config.run_dir / "inputs" / "equilibration.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return output_pdb


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Equilibrate one manuscript complex with Jorgensen-inspired OpenMM MD"
    )
    parser.add_argument("--mutation", default="V106A")
    parser.add_argument("--start-label", default="WT")
    parser.add_argument("--end-label")
    parser.add_argument("--input-complex-pdb", "--wt-complex-pdb", dest="input_complex_pdb", type=Path)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    parser.add_argument("--platform")
    args = parser.parse_args()
    from .mutations import MutationLeg

    end_label = args.end_label or args.mutation
    leg = MutationLeg(args.start_label, end_label, args.mutation)
    config = FEPConfig.for_leg(
        leg,
        wt_complex_pdb=args.input_complex_pdb or leg.input_complex_pdb(args.replicate),
        output_dir=args.output_dir,
    )
    output = equilibrate_leg(config, platform_name=args.platform)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
