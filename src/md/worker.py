#!/usr/bin/env python
"""Cluster worker for explicit-MD task execution (no alchemical protocol)."""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

from .manifest import MDTask, get_task_by_id


def run_md_task(
    task: MDTask,
    heating_ps: float,
    production_ns: float,
    report_interval: int,
    checkpoint_interval: int,
    resume_from_checkpoint: bool,
    force: bool = False,
) -> dict:
    from .openmm.md_protocol import MDProtocolConfig, run_prepared_md

    output_path = Path(task.output_json)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    final_pdb = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_md_final.pdb"
    state_csv = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_md_state.csv"
    analysis_dcd = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_analysis.dcd"
    analysis_topo_pdb = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_analysis_topology.pdb"
    checkpoint_path = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_md.chk"

    if not force and output_path.exists():
        try:
            existing = json.loads(output_path.read_text())
            if str(existing.get("status", "")).lower() == "ok":
                logging.info(
                    "Task %d already completed (%s rep%d); skipping (use --force to rerun).",
                    task.task_id,
                    task.mutation,
                    task.replicate,
                )
                return existing
        except Exception:
            pass

    # ~200 analysis frames (enough for ~100 after 25% discard).
    timestep_fs = 2.0
    production_steps = max(1, int(round((production_ns * 1_000_000.0) / timestep_fs)))
    analysis_interval = max(1, production_steps // 200)

    cfg = MDProtocolConfig(
        heating_ps=heating_ps,
        production_ns=production_ns,
        report_interval_steps=report_interval,
        analysis_report_interval_steps=analysis_interval,
    )

    if not task.prepared_system_xml or not task.prepared_topology_pdb:
        raise ValueError("Task is missing prepared MD assets (prepared_system_xml/prepared_topology_pdb).")

    result = run_prepared_md(
        prepared_topology_pdb=Path(task.prepared_topology_pdb),
        prepared_system_xml=Path(task.prepared_system_xml),
        final_pdb_path=final_pdb,
        state_csv_path=state_csv,
        config=cfg,
        analysis_dcd_path=analysis_dcd,
        analysis_topology_pdb_path=analysis_topo_pdb,
        checkpoint_path=checkpoint_path,
        checkpoint_interval_steps=max(1, int(checkpoint_interval)),
        resume_from_checkpoint=bool(resume_from_checkpoint),
    )

    payload = {
        "task_id": task.task_id,
        "structure": task.structure,
        "mutation": task.mutation,
        "safe_label": task.safe_label,
        "replicate": task.replicate,
        "fold_reduction": task.fold_reduction,
        "ligand_resname": task.ligand_resname,
        "ligand_sdf": task.ligand_sdf,
        "minimized_pdb": task.minimized_pdb,
        "prepared_topology_pdb": task.prepared_topology_pdb,
        "prepared_system_xml": task.prepared_system_xml,
        "analysis_dcd": str(analysis_dcd),
        "analysis_topology_pdb": str(analysis_topo_pdb),
        "state_csv": str(state_csv),
        "final_pdb": str(final_pdb),
        "checkpoint_path": str(checkpoint_path),
        "md_heating_steps": result.heating_steps,
        "md_production_steps": result.production_steps,
        "md_production_steps_completed": result.production_steps_completed,
        "md_total_steps": result.total_steps,
        "resumed_from_checkpoint": result.resumed_from_checkpoint,
        "elapsed_seconds": result.elapsed_seconds,
        "status": "ok",
    }
    output_path.write_text(json.dumps(payload, indent=2))
    return payload


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Run one explicit-MD task from manifest")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--heating-ps", type=float, default=25.0)
    parser.add_argument("--production-ns", type=float, default=2.0)
    parser.add_argument("--report-interval", type=int, default=2000)
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=int(os.environ.get("MD_CHECKPOINT_INTERVAL", "5000")),
        help="Write an OpenMM checkpoint every N steps.",
    )
    parser.add_argument(
        "--resume-from-checkpoint",
        type=int,
        default=int(os.environ.get("MD_RESUME_FROM_CHECKPOINT", "1")),
        help="If 1, resume from existing checkpoint when present.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore existing successful output JSON and rerun the task.",
    )
    args = parser.parse_args(argv)

    try:
        task = get_task_by_id(args.manifest, args.task_id)
    except Exception as exc:
        logging.error("Failed to load task %s: %s", args.task_id, exc)
        return 1

    try:
        out = run_md_task(
            task,
            heating_ps=args.heating_ps,
            production_ns=args.production_ns,
            report_interval=args.report_interval,
            checkpoint_interval=max(1, int(args.checkpoint_interval)),
            resume_from_checkpoint=bool(int(args.resume_from_checkpoint)),
            force=bool(args.force),
        )
        logging.info(
            "Completed task %d (%s rep%d) in %.1fs",
            task.task_id,
            task.mutation,
            task.replicate,
            out["elapsed_seconds"],
        )
    except Exception as exc:
        logging.error("Task %d failed: %s", task.task_id, exc, exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
