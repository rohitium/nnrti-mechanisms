"""Perses hybrid-topology preparation for holo mutation FEP legs."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from src.analysis.cli.run_perses_point_mutation_cycle import extract_protein_and_ligand

from .openeye_shim import install_openeye_shim

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


def _write_holo_phase(htf, config: FEPConfig) -> None:
    from openmm import XmlSerializer, app

    holo_dir = config.run_dir / "holo"
    holo_dir.mkdir(parents=True, exist_ok=True)
    (holo_dir / "hybrid_system.xml").write_text(XmlSerializer.serialize(htf.hybrid_system))
    schedule = {
        "phase": "holo",
        "mutation": config.mutation,
        "start_label": config.start_label,
        "end_label": config.end_label,
        "leg_id": config.leg.leg_id,
        "lambda_values": list(config.lambda_schedule.values),
        "lambda_parameter_functions": "perses-default",
        "thermodynamic_cycle": config.approx_protocol.thermodynamic_cycle,
        "prepare_backend": "perses",
        "input_complex_pdb": str(config.preparation_complex_pdb),
    }
    (holo_dir / "schedule.json").write_text(json.dumps(schedule, indent=2) + "\n")
    pdb_path = holo_dir / "hybrid_topology.pdb"
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


def prepare_holo_hybrid(config: FEPConfig, replicate: int = 1) -> None:
    """Build a Perses HybridTopologyFactory for one holo mutation leg."""
    _configure_perses_runtime()
    install_openeye_shim()
    from openmm import MonteCarloBarostat, app, unit
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

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
    executor = PointMutationExecutor(
        protein_filename=str(protein_pdb),
        mutation_chain_id=config.chain_id,
        mutation_residue_id=site.pdb_residue_id,
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
        conduct_endstate_validation=False,
        generate_unmodified_hybrid_topology_factory=True,
        generate_rest_capable_hybrid_topology_factory=False,
    )
    _write_holo_phase(executor.get_complex_htf(), config)
    config.write(config.run_dir / "config.json")
    config.approx_protocol.write(config.run_dir / "approx_protocol.json")
    (config.run_dir / "prepare_backend.json").write_text(
        json.dumps(
            {
                "backend": "perses",
                "replicate": replicate,
                "protein_pdb": str(protein_pdb),
                "ligand_sdf": str(bound_ligand_sdf),
                "input_complex_pdb": str(config.preparation_complex_pdb),
                "auth_residue_id": config.residue_id,
                "pdb_residue_id": site.pdb_residue_id,
            },
            indent=2,
        )
        + "\n"
    )


def perses_available() -> bool:
    _configure_perses_runtime()
    install_openeye_shim()
    try:
        import openmmtools  # noqa: F401
        import perses  # noqa: F401
        from openff.toolkit import Molecule  # noqa: F401
        from rdkit import Chem  # noqa: F401
    except ImportError:
        return False
    return True
