from __future__ import annotations

import logging
import time
from pathlib import Path

from .analysis_metrics import ContactMetrics, compute_contacts, pocket_volume_proxy
from .openmm.pipeline import (
    build_forcefield,
    compute_binding_proxy,
    load_ligand_molecule,
    minimize_with_restraints,
)
from .openmm.alchemy import (
    AlchemicalConfig,
    compute_alchemical_binding_free_energy,
    write_alchemical_result,
)
from .openmm.minimizer import minimize_system
from .openmm.require import require_module


def prepare_structure(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius: float,
    restraint_k: float,
    output_path: Path,
    jitter_seed: int | None = None,
    jitter_angstrom: float = 0.0,
):
    """Prepare structure with minimization and compute all metrics including binding proxy.

    This is the legacy function that computes binding proxy in addition to structural metrics.
    """
    if output_path.exists():
        logging.info("Reusing existing minimized structure: %s", output_path)
        app = require_module("openmm.app")
        with open(output_path, "r") as handle:
            pdb = app.PDBFile(handle)
        ligand = load_ligand_molecule(ligand_sdf)
        forcefield = build_forcefield([ligand])
        energies = compute_binding_proxy(
            pdb.topology, pdb.positions, forcefield, ligand_resname=ligand_resname
        )
        contacts = compute_contacts(output_path, ligand_resname=ligand_resname)
        pocket = pocket_volume_proxy(output_path, ligand_resname=ligand_resname)
        return energies, contacts, pocket

    start = time.perf_counter()
    topology, positions, forcefield = minimize_with_restraints(
        cif_path=cif_path,
        ligand_resname=ligand_resname,
        ligand_sdf=ligand_sdf,
        restraint_radius_angstrom=restraint_radius,
        restraint_k_kj_mol_nm2=restraint_k,
        output_path=output_path,
        jitter_seed=jitter_seed,
        jitter_angstrom=jitter_angstrom,
    )
    logging.info("Minimized structure in %.2fs", time.perf_counter() - start)

    start = time.perf_counter()
    _, positions = minimize_system(
        topology,
        positions,
        forcefield,
        restraint_indices=[],
        restraint_k_kj_mol_nm2=0.0,
    )
    logging.info(
        "Unrestrained minimization completed in %.2fs", time.perf_counter() - start
    )
    app = require_module("openmm.app")
    with open(output_path, "w") as handle:
        app.PDBFile.writeFile(topology, positions, handle)

    start = time.perf_counter()
    energies = compute_binding_proxy(
        topology, positions, forcefield, ligand_resname=ligand_resname
    )
    logging.info("Energy proxy computed in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    metrics_path = output_path
    contacts = compute_contacts(metrics_path, ligand_resname=ligand_resname)
    logging.info("Contacts computed in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    pocket = pocket_volume_proxy(metrics_path, ligand_resname=ligand_resname)
    logging.info("Pocket volume computed in %.2fs", time.perf_counter() - start)
    return energies, contacts, pocket


def prepare_structure_local(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius: float,
    restraint_k: float,
    output_path: Path,
    jitter_seed: int | None = None,
    jitter_angstrom: float = 0.0,
) -> tuple[ContactMetrics, float]:
    """Prepare structure with minimization and compute structural metrics only.

    This function is for the cluster workflow where FEP is run separately.
    It returns only structural metrics (contacts, H-bonds, pocket volume),
    not binding proxy or alchemical calculations.

    Args:
        cif_path: Path to input CIF structure.
        ligand_resname: Residue name of the ligand.
        ligand_sdf: Path to ligand SDF file.
        restraint_radius: Radius for restraints during minimization (angstrom).
        restraint_k: Force constant for restraints (kJ/mol/nm^2).
        output_path: Path to write minimized PDB.
        jitter_seed: Optional seed for coordinate jitter.
        jitter_angstrom: Magnitude of coordinate jitter (angstrom).

    Returns:
        Tuple of (ContactMetrics, pocket_volume_proxy).
    """
    if output_path.exists():
        logging.info("Reusing existing minimized structure: %s", output_path)
        contacts = compute_contacts(output_path, ligand_resname=ligand_resname)
        pocket = pocket_volume_proxy(output_path, ligand_resname=ligand_resname)
        return contacts, pocket

    start = time.perf_counter()
    topology, positions, forcefield = minimize_with_restraints(
        cif_path=cif_path,
        ligand_resname=ligand_resname,
        ligand_sdf=ligand_sdf,
        restraint_radius_angstrom=restraint_radius,
        restraint_k_kj_mol_nm2=restraint_k,
        output_path=output_path,
        jitter_seed=jitter_seed,
        jitter_angstrom=jitter_angstrom,
    )
    logging.info("Minimized structure in %.2fs", time.perf_counter() - start)

    start = time.perf_counter()
    _, positions = minimize_system(
        topology,
        positions,
        forcefield,
        restraint_indices=[],
        restraint_k_kj_mol_nm2=0.0,
    )
    logging.info(
        "Unrestrained minimization completed in %.2fs", time.perf_counter() - start
    )

    app = require_module("openmm.app")
    with open(output_path, "w") as handle:
        app.PDBFile.writeFile(topology, positions, handle)

    start = time.perf_counter()
    contacts = compute_contacts(output_path, ligand_resname=ligand_resname)
    logging.info("Contacts computed in %.2fs", time.perf_counter() - start)

    start = time.perf_counter()
    pocket = pocket_volume_proxy(output_path, ligand_resname=ligand_resname)
    logging.info("Pocket volume computed in %.2fs", time.perf_counter() - start)

    return contacts, pocket


def compute_alchemical_binding_metric(
    minimized_pdb_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    config: AlchemicalConfig,
    output_json: Path | None = None,
    metadata: dict | None = None,
) -> float:
    """Compute alchemical binding free energy from a minimized structure."""
    app = require_module("openmm.app")
    with open(minimized_pdb_path, "r") as handle:
        pdb = app.PDBFile(handle)
    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])
    result = compute_alchemical_binding_free_energy(
        topology=pdb.topology,
        positions=pdb.positions,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=config,
    )
    if output_json is not None:
        write_alchemical_result(
            output_json, result=result, metadata=metadata or {}
        )
    return float(result.binding_delta_g_kj_mol)
