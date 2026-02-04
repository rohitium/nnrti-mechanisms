#!/usr/bin/env python
"""FEP worker for SLURM array job execution.

This module is the entry point for each SLURM array task. It:
1. Minimizes the structure from input CIF (with optional jitter)
2. Runs a single FEP leg (complex or solvent)
3. Saves results to JSON

Usage:
    python -m src.cluster.fep_worker --manifest results/fep_manifest.csv --task-id 0
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

from .manifest import get_task_by_id, FEPTask


def minimize_structure(task: FEPTask, output_dir: Path) -> Path:
    """Minimize the input CIF structure.

    Args:
        task: The FEPTask containing minimization parameters.
        output_dir: Directory to write minimized PDB.

    Returns:
        Path to the minimized PDB file.
    """
    from ..openmm.pipeline import (
        build_forcefield,
        load_ligand_molecule,
        minimize_with_restraints,
    )
    from ..openmm.minimizer import minimize_system
    from ..openmm.require import require_module

    min_pdb = output_dir / f"{task.safe_label}_minimized_rep{task.replicate:02d}.pdb"

    if min_pdb.exists():
        logging.info("Reusing existing minimized structure: %s", min_pdb)
        return min_pdb

    logging.info("Minimizing structure from %s", task.input_cif)
    start = time.perf_counter()

    topology, positions, forcefield = minimize_with_restraints(
        cif_path=Path(task.input_cif),
        ligand_resname=task.ligand_resname,
        ligand_sdf=Path(task.ligand_sdf),
        restraint_radius_angstrom=task.restraint_radius,
        restraint_k_kj_mol_nm2=task.restraint_k,
        output_path=min_pdb,
        jitter_seed=task.jitter_seed,
        jitter_angstrom=task.jitter_angstrom,
    )

    logging.info("Restrained minimization done in %.1fs", time.perf_counter() - start)

    # Second unrestrained minimization
    start = time.perf_counter()
    _, positions = minimize_system(
        topology,
        positions,
        forcefield,
        restraint_indices=[],
        restraint_k_kj_mol_nm2=0.0,
    )
    logging.info("Unrestrained minimization done in %.1fs", time.perf_counter() - start)

    # Write final structure
    app = require_module("openmm.app")
    with open(min_pdb, "w") as handle:
        app.PDBFile.writeFile(topology, positions, handle)

    return min_pdb


def run_fep_task(
    task: FEPTask,
    equil_steps: int = 10_000,
    prod_steps: int = 25_000,
    sample_interval: int = 200,
) -> dict:
    """Execute a single FEP leg task.

    Args:
        task: The FEPTask to execute.
        equil_steps: Equilibration steps per lambda window.
        prod_steps: Production steps per lambda window.
        sample_interval: Sample interval for energy evaluations.

    Returns:
        Dictionary with task result including delta_g_kj_mol.
    """
    from ..openmm.alchemy import AlchemicalConfig, run_single_leg

    output_path = Path(task.output_json)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Minimize structure (if not already done)
    min_pdb = minimize_structure(task, output_dir)

    # Step 2: Run FEP leg
    config = AlchemicalConfig(
        equilibration_steps=equil_steps,
        production_steps=prod_steps,
        sample_interval=sample_interval,
    )

    metadata = {
        "task_id": task.task_id,
        "structure": task.structure,
        "mutation": task.mutation,
        "safe_label": task.safe_label,
        "replicate": task.replicate,
        "leg": task.leg,
        "fold_reduction": task.fold_reduction,
        "minimized_pdb": str(min_pdb),
    }

    logging.info(
        "Running FEP task %d: %s rep%d %s leg",
        task.task_id,
        task.mutation,
        task.replicate,
        task.leg,
    )

    start_time = time.perf_counter()

    result = run_single_leg(
        minimized_pdb_path=min_pdb,
        ligand_resname=task.ligand_resname,
        ligand_sdf=Path(task.ligand_sdf),
        leg=task.leg,
        config=config,
        output_json=output_path,
        metadata=metadata,
    )

    elapsed = time.perf_counter() - start_time
    logging.info(
        "Task %d completed in %.1f s: dG = %.2f kJ/mol",
        task.task_id,
        elapsed,
        result.delta_g_kj_mol,
    )

    return {
        "task_id": task.task_id,
        "delta_g_kj_mol": result.delta_g_kj_mol,
        "elapsed_seconds": elapsed,
        "minimized_pdb": str(min_pdb),
    }


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the FEP worker."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Run a single FEP leg for SLURM array job"
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Path to FEP manifest CSV",
    )
    parser.add_argument(
        "--task-id",
        type=int,
        required=True,
        help="Task ID from SLURM_ARRAY_TASK_ID",
    )
    parser.add_argument(
        "--equil-steps",
        type=int,
        default=10_000,
        help="Equilibration steps per lambda window (default: 10000)",
    )
    parser.add_argument(
        "--prod-steps",
        type=int,
        default=25_000,
        help="Production steps per lambda window (default: 25000)",
    )
    parser.add_argument(
        "--sample-interval",
        type=int,
        default=200,
        help="Sample interval for energy evaluations (default: 200)",
    )

    args = parser.parse_args(argv)

    try:
        task = get_task_by_id(args.manifest, args.task_id)
    except ValueError as e:
        logging.error("Failed to load task: %s", e)
        return 1

    try:
        run_fep_task(
            task,
            equil_steps=args.equil_steps,
            prod_steps=args.prod_steps,
            sample_interval=args.sample_interval,
        )
    except Exception as e:
        logging.error("Task %d failed: %s", args.task_id, e, exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
