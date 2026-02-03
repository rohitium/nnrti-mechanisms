#!/usr/bin/env python
"""FEP worker for SLURM array job execution.

This module is the entry point for each SLURM array task. It runs a single
FEP leg (complex or solvent) for a specific mutation/replicate combination.

Usage:
    python -m src.cluster.fep_worker --manifest results/fep_manifest.csv --task-id 0
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from .manifest import get_task_by_id, FEPTask


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
    }

    output_path = Path(task.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logging.info(
        "Running FEP task %d: %s rep%d %s leg",
        task.task_id,
        task.mutation,
        task.replicate,
        task.leg,
    )

    start_time = time.perf_counter()

    result = run_single_leg(
        minimized_pdb_path=Path(task.minimized_pdb),
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
