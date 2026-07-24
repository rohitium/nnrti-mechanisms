"""Perses hybrid-topology preparation for holo and apo mutation FEP legs."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from src.analysis.cli.run_perses_point_mutation_cycle import extract_protein_and_ligand

from .openeye_shim import install_openeye_shim
from .perses_patches import patch_perses_proline_support
from .structure import extract_protein_only

from .alchemical import resolve_mutation_site
from .config import FEPConfig
from .mutations import Mutation


def _configure_perses_runtime() -> None:
    """Ensure AMBER GAFF is discoverable before openmoltools import."""
    if os.environ.get("AMBERHOME"):
        return
    candidates: list[Path] = []
    for value in (sys.prefix, os.environ.get("CONDA_PREFIX")):
        if value:
            prefix = Path(value)
            if prefix not in candidates:
                candidates.append(prefix)
    for prefix in candidates:
        gaff_dat = prefix / "dat" / "leap" / "parm" / "gaff.dat"
        if not gaff_dat.is_file():
            continue
        os.environ["AMBERHOME"] = str(prefix)
        bindir = prefix / "bin"
        if bindir.is_dir():
            path = os.environ.get("PATH", "")
            bindir_str = str(bindir)
            if bindir_str not in path.split(os.pathsep):
                os.environ["PATH"] = bindir_str + os.pathsep + path
        return


def _remove_scaling_artifacts(run_dir: Path) -> None:
    """Drop MD-asset scaling outputs if a leg is re-prepared with Perses."""
    scaling_plan = run_dir / "alchemical_plan.json"
    if scaling_plan.is_file():
        scaling_plan.unlink()


def _ensure_ligand_sdf(config: FEPConfig, input_dir: Path) -> Path:
    ligand_sdf = input_dir / "dor_bound_pose.sdf"
    if ligand_sdf.is_file():
        return ligand_sdf
    _, ligand_sdf = extract_protein_and_ligand(
        config.wt_complex_pdb,
        config.ligand_sdf,
        input_dir,
        ligand_resname=config.ligand_resname,
    )
    return ligand_sdf


def _build_point_mutation_executor(
    config: FEPConfig,
    *,
    protein_pdb: Path,
    ligand_sdf: Path,
    pdb_residue_id: str,
):
    from openmm import MonteCarloBarostat, app, unit
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

    return PointMutationExecutor(
        protein_filename=str(protein_pdb),
        mutation_chain_id=config.chain_id,
        mutation_residue_id=pdb_residue_id,
        proposed_residue=config.mutant_residue,
        old_residue=config.wt_residue,
        ligand_input=str(ligand_sdf),
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
    )


def _write_phase(htf, config: FEPConfig, phase: str, source_pdb: Path) -> None:
    from openmm import XmlSerializer, app

    if phase == "holo":
        _remove_scaling_artifacts(config.run_dir)
    phase_dir = config.run_dir / phase
    phase_dir.mkdir(parents=True, exist_ok=True)
    (phase_dir / "hybrid_system.xml").write_text(XmlSerializer.serialize(htf.hybrid_system))
    schedule = {
        "phase": phase,
        "mutation": config.mutation,
        "start_label": config.start_label,
        "end_label": config.end_label,
        "leg_id": config.leg.leg_id,
        "lambda_values": list(config.lambda_schedule.values),
        "lambda_parameter_functions": "perses-default",
        "thermodynamic_cycle": config.approx_protocol.thermodynamic_cycle,
        "prepare_backend": "perses",
        "source_pdb": str(source_pdb),
    }
    (phase_dir / "schedule.json").write_text(json.dumps(schedule, indent=2) + "\n")
    pdb_path = phase_dir / "hybrid_topology.pdb"
    try:
        topology = htf.omm_hybrid_topology
        topology.setPeriodicBoxVectors(htf.hybrid_system.getDefaultPeriodicBoxVectors())
        with pdb_path.open("w") as handle:
            app.PDBFile.writeFile(topology, htf.hybrid_positions, handle, keepIds=True)
    except Exception:
        import mdtraj as md
        import numpy as np
        from openmm import unit

        xyz = np.array(htf.hybrid_positions.value_in_unit(unit.nanometer))[None, ...]
        md.Trajectory(xyz, topology=htf.hybrid_topology).save_pdb(str(pdb_path))


def _write_prepare_metadata(
    config: FEPConfig,
    replicate: int,
    *,
    protein_pdb: Path,
    ligand_sdf: Path | None,
    source_pdb: Path,
    pdb_residue_id: str,
    phases: tuple[str, ...],
) -> None:
    (config.run_dir / "prepare_backend.json").write_text(
        json.dumps(
            {
                "backend": "perses",
                "replicate": replicate,
                "phases": list(phases),
                "protein_pdb": str(protein_pdb),
                "ligand_sdf": str(ligand_sdf) if ligand_sdf is not None else None,
                "source_pdb": str(source_pdb),
                "auth_residue_id": config.residue_id,
                "pdb_residue_id": pdb_residue_id,
            },
            indent=2,
        )
        + "\n"
    )


def prepare_holo_hybrid(config: FEPConfig, replicate: int = 1, *, write_metadata: bool = True) -> None:
    """Build a Perses hybrid for mutation in the inhibitor-bound complex."""
    _configure_perses_runtime()
    install_openeye_shim()
    patch_perses_proline_support()

    config.validate(require_inputs=True)
    mutation = Mutation.parse(config.mutation)
    site = resolve_mutation_site(
        config.leg.input_complex_pdb(replicate=replicate),
        config.leg.endpoint_complex_pdb(replicate=replicate),
        mutation,
        chain_id=config.chain_id,
    )
    input_dir = config.run_dir / "inputs"
    protein_pdb, bound_ligand_sdf = extract_protein_and_ligand(
        config.preparation_complex_pdb,
        config.ligand_sdf,
        input_dir,
        ligand_resname=config.ligand_resname,
    )
    executor = _build_point_mutation_executor(
        config,
        protein_pdb=protein_pdb,
        ligand_sdf=bound_ligand_sdf,
        pdb_residue_id=site.pdb_residue_id,
    )
    _write_phase(executor.get_complex_htf(), config, "holo", config.preparation_complex_pdb)
    if write_metadata:
        config.write(config.run_dir / "config.json")
        config.approx_protocol.write(config.run_dir / "approx_protocol.json")
        _write_prepare_metadata(
            config,
            replicate,
            protein_pdb=protein_pdb,
            ligand_sdf=bound_ligand_sdf,
            source_pdb=config.preparation_complex_pdb,
            pdb_residue_id=site.pdb_residue_id,
            phases=("holo",),
        )


def prepare_apo_hybrid(config: FEPConfig, replicate: int = 1) -> None:
    """Build a Perses hybrid for the same mutation in the apo protein."""
    _configure_perses_runtime()
    install_openeye_shim()
    patch_perses_proline_support()

    config.validate(require_inputs=True)
    mutation = Mutation.parse(config.mutation)
    apo_source = config.leg.input_apo_pdb(replicate=replicate)
    if not apo_source.is_file():
        raise FileNotFoundError(apo_source)
    site = resolve_mutation_site(
        apo_source,
        config.leg.endpoint_apo_pdb(replicate=replicate),
        mutation,
        chain_id=config.chain_id,
    )
    input_dir = config.run_dir / "inputs"
    protein_pdb = extract_protein_only(
        apo_source,
        input_dir,
        ligand_resname=config.ligand_resname,
        output_name="apo_protein_no_ligand.pdb",
    )
    bound_ligand_sdf = _ensure_ligand_sdf(config, input_dir)
    executor = _build_point_mutation_executor(
        config,
        protein_pdb=protein_pdb,
        ligand_sdf=bound_ligand_sdf,
        pdb_residue_id=site.pdb_residue_id,
    )
    _write_phase(executor.get_apo_htf(), config, "apo", apo_source)


def prepare_hybrid_leg(
    config: FEPConfig,
    replicate: int = 1,
    phases: tuple[str, ...] = ("holo", "apo"),
) -> None:
    """Prepare one or both thermodynamic phases for a mutation leg."""
    selected = tuple(dict.fromkeys(phases))
    config.validate(require_inputs=True)
    mutation = Mutation.parse(config.mutation)
    holo_site = resolve_mutation_site(
        config.leg.input_complex_pdb(replicate=replicate),
        config.leg.endpoint_complex_pdb(replicate=replicate),
        mutation,
        chain_id=config.chain_id,
    )
    input_dir = config.run_dir / "inputs"

    if "holo" in selected:
        prepare_holo_hybrid(config, replicate=replicate, write_metadata=False)
    if "apo" in selected:
        prepare_apo_hybrid(config, replicate=replicate)

    protein_pdb = input_dir / "wt_receptor_no_ligand.pdb"
    if not protein_pdb.is_file():
        protein_pdb = input_dir / "apo_protein_no_ligand.pdb"
    ligand_sdf = input_dir / "dor_bound_pose.sdf"
    config.write(config.run_dir / "config.json")
    config.approx_protocol.write(config.run_dir / "approx_protocol.json")
    _write_prepare_metadata(
        config,
        replicate,
        protein_pdb=protein_pdb,
        ligand_sdf=ligand_sdf if ligand_sdf.is_file() else None,
        source_pdb=config.preparation_complex_pdb,
        pdb_residue_id=holo_site.pdb_residue_id,
        phases=selected,
    )


def perses_available() -> bool:
    _configure_perses_runtime()
    install_openeye_shim()
    patch_perses_proline_support()
    try:
        import openmmtools  # noqa: F401
        import perses  # noqa: F401
        from openff.toolkit import Molecule  # noqa: F401
        from rdkit import Chem  # noqa: F401
    except ImportError:
        return False
    return True
