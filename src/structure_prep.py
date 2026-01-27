from __future__ import annotations

import logging
import time
from pathlib import Path

from .analysis_metrics import compute_contacts, pocket_volume_proxy
from .openmm.pipeline import (
    build_forcefield,
    compute_binding_proxy,
    load_ligand_molecule,
    minimize_with_restraints,
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
):
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
