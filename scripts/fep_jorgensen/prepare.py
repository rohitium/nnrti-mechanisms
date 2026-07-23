from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import FEPConfig
from .equilibrate import equilibrate_leg
from .mutations import MutationLeg
from src.analysis.cli.run_perses_point_mutation_cycle import extract_protein_and_ligand


def _write_phase(htf, phase: str, config: FEPConfig) -> None:
    """Serialize Perses hybrid products for OpenMM-only workers."""
    from openmm import XmlSerializer, app

    phase_dir = config.run_dir / phase
    phase_dir.mkdir(parents=True, exist_ok=True)
    (phase_dir / "hybrid_system.xml").write_text(
        XmlSerializer.serialize(htf.hybrid_system)
    )
    with (phase_dir / "hybrid_topology.pdb").open("w") as handle:
        app.PDBFile.writeFile(
            htf.hybrid_topology, htf.hybrid_positions, handle, keepIds=True
        )
    schedule = {
        "phase": phase,
        "mutation": config.mutation,
        "start_label": config.start_label,
        "end_label": config.end_label,
        "leg_id": config.leg.leg_id,
        "lambda_values": list(config.lambda_schedule.values),
        "lambda_parameter_functions": "perses-default",
        "thermodynamic_cycle": config.approx_protocol.thermodynamic_cycle,
    }
    (phase_dir / "schedule.json").write_text(json.dumps(schedule, indent=2) + "\n")


def prepare(config: FEPConfig, platform_name: str = "CPU") -> None:
    """Local alchemical setup for holo-only mutation FEP.

    Flow:
    1. Jorgensen-inspired OpenMM equilibration of the bound complex.
    2. Perses hybrid construction for the protein-side-chain mutation leg.
    3. Serialization of the holo hybrid for Sherlock/OpenMM sampling workers.
    """
    from openmm import MonteCarloBarostat, app, unit
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

    config.validate(require_inputs=True)
    if not config.skip_equilibration:
        equilibrate_leg(config, platform_name=platform_name)
    protein_pdb, bound_ligand_sdf = extract_protein_and_ligand(
        config.preparation_complex_pdb,
        config.ligand_sdf,
        config.run_dir / "inputs",
        ligand_resname=config.ligand_resname,
    )
    executor = PointMutationExecutor(
        protein_filename=str(protein_pdb),
        mutation_chain_id=config.chain_id,
        mutation_residue_id=config.residue_id,
        proposed_residue=config.mutant_residue,
        old_residue=config.wt_residue,
        ligand_input=str(bound_ligand_sdf),
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
        conduct_endstate_validation=True,
        generate_unmodified_hybrid_topology_factory=True,
        generate_rest_capable_hybrid_topology_factory=False,
    )
    _write_phase(executor.get_complex_htf(), "holo", config)
    config.write(config.run_dir / "config.json")
    config.approx_protocol.write(config.run_dir / "approx_protocol.json")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Equilibrate and prepare one holo mutation FEP leg locally"
    )
    parser.add_argument("--mutation", default="V106A", help="Single substitution made in this leg")
    parser.add_argument("--start-label", default="WT")
    parser.add_argument("--end-label", help="Resulting single or compound mutant label")
    parser.add_argument("--input-complex-pdb", "--wt-complex-pdb", dest="input_complex_pdb", type=Path)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    parser.add_argument("--skip-equilibration", action="store_true")
    parser.add_argument("--platform", default="CPU")
    args = parser.parse_args()
    end_label = args.end_label or args.mutation
    leg = MutationLeg(args.start_label, end_label, args.mutation)
    config = FEPConfig.for_leg(
        leg,
        wt_complex_pdb=args.input_complex_pdb or leg.input_complex_pdb(args.replicate),
        output_dir=args.output_dir,
        skip_equilibration=args.skip_equilibration,
    )
    prepare(config, platform_name=args.platform)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
