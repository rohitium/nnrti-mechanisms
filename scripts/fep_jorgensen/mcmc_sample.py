"""Integrated openmmtools MCMC sampling for one Perses-prepared holo leg."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from .config import FEPConfig
from .mutations import MutationLeg
from .perses_hybrid import perses_available, prepare_holo_hybrid


def _run_holo_mcmc(htf, run_dir: Path, config: FEPConfig, args) -> Path:
    from openmm import unit
    from openmmtools import cache, mcmc, utils
    from openmmtools.multistate import MultiStateReporter
    from perses.annihilation.lambda_protocol import LambdaProtocol
    from perses.dispersed.utils import configure_platform
    from perses.samplers.multistate import HybridRepexSampler

    reporter_path = run_dir / "holo" / "multistate.nc"
    reporter = MultiStateReporter(
        str(reporter_path),
        checkpoint_interval=args.checkpoint_interval,
    )
    move = mcmc.LangevinSplittingDynamicsMove(
        timestep=args.timestep_fs * unit.femtoseconds,
        collision_rate=args.collision_rate / unit.picosecond,
        n_steps=args.steps_per_cycle,
        reassign_velocities=False,
        n_restart_attempts=20,
        splitting="V R R R O R R R V",
        constraint_tolerance=1e-6,
    )
    sampler = HybridRepexSampler(mcmc_moves=move, hybrid_factory=htf)
    sampler.setup(
        n_states=len(config.lambda_schedule.values),
        temperature=config.temperature_k * unit.kelvin,
        storage_file=reporter,
        lambda_protocol=LambdaProtocol(functions="default"),
        endstates=False,
    )
    platform = configure_platform(args.platform or utils.get_fastest_platform().getName())
    sampler.energy_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.sampler_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.extend(args.n_cycles)
    reporter.close()
    return reporter_path


def _load_hybrid_factory(config: FEPConfig):
    from openmm import MonteCarloBarostat, app, unit
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

    input_dir = config.run_dir / "inputs"
    return PointMutationExecutor(
        protein_filename=str(input_dir / "wt_receptor_no_ligand.pdb"),
        mutation_chain_id=config.chain_id,
        mutation_residue_id=config.residue_id,
        proposed_residue=config.mutant_residue,
        old_residue=config.wt_residue,
        ligand_input=str(input_dir / "dor_bound_pose.sdf"),
        ligand_index=0,
        allow_undefined_stereo_sdf=True,
        is_solvated=False,
        forcefield_files=[
            config.approx_protocol.force_field_protein,
            config.approx_protocol.force_field_dna,
            config.approx_protocol.force_field_water,
        ],
        small_molecule_forcefields=config.approx_protocol.small_molecule_forcefield,
        barostat=MonteCarloBarostat(
            config.pressure_atm * unit.atmosphere,
            config.temperature_k * unit.kelvin,
            25,
        ),
        forcefield_kwargs={
            "removeCMMotion": False,
            "constraints": app.HBonds,
            "hydrogenMass": 1.0 * unit.amu,
        },
        periodic_forcefield_kwargs={
            "nonbondedMethod": app.PME,
            "nonbondedCutoff": 1.0 * unit.nanometer,
            "ewaldErrorTolerance": 1e-4,
        },
        conduct_endstate_validation=False,
        generate_unmodified_hybrid_topology_factory=True,
        generate_rest_capable_hybrid_topology_factory=False,
    ).get_complex_htf()


def sample_holo_leg(config: FEPConfig, args) -> Path:
    holo_dir = config.run_dir / "holo"
    input_dir = config.run_dir / "inputs"
    if not (holo_dir / "hybrid_system.xml").is_file():
        if not perses_available():
            raise ImportError("Perses/openmmtools required for mcmc_sample")
        prepare_holo_hybrid(config)
    elif not input_dir.is_dir():
        prepare_holo_hybrid(config)
    htf = _load_hybrid_factory(config)
    return _run_holo_mcmc(htf, config.run_dir, config, args)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run openmmtools multistate MCMC on one Perses-prepared holo leg"
    )
    parser.add_argument("--mutation", default="V106A")
    parser.add_argument("--start-label", default="WT")
    parser.add_argument("--end-label")
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    parser.add_argument("--n-cycles", type=int, default=5000)
    parser.add_argument("--steps-per-cycle", type=int, default=250)
    parser.add_argument("--timestep-fs", type=float, default=4.0)
    parser.add_argument("--collision-rate", type=float, default=5.0)
    parser.add_argument("--checkpoint-interval", type=int, default=10)
    parser.add_argument("--platform", default="CUDA")
    args = parser.parse_args()
    end_label = args.end_label or args.mutation
    leg = MutationLeg(args.start_label, end_label, args.mutation)
    config = FEPConfig.for_leg(leg, output_dir=args.output_dir, prepare_backend="perses")
    start = time.time()
    reporter = sample_holo_leg(config, args)
    summary = {
        "leg_id": config.leg.leg_id,
        "mutation": config.mutation,
        "reporter": str(reporter),
        "n_cycles": args.n_cycles,
        "steps_per_cycle": args.steps_per_cycle,
        "elapsed_seconds": time.time() - start,
        "sampler": "openmmtools HybridRepexSampler",
    }
    (config.run_dir / "mcmc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
