#!/usr/bin/env python
"""Cluster worker for explicit-MD task execution (no alchemical protocol)."""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .manifest import MDTask, get_task_by_id


def run_md_task(
    task: MDTask,
    heating_ps: float,
    production_ns: float,
    report_interval: int,
) -> dict:
    from ..openmm.md_protocol import MDProtocolConfig, run_prepared_md

    output_path = Path(task.output_json)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    final_pdb = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_md_final.pdb"
    state_csv = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_md_state.csv"
    analysis_dcd = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_analysis.dcd"
    analysis_topo_pdb = output_dir / f"{task.safe_label}_rep{task.replicate:02d}_analysis_topology.pdb"

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
        "md_heating_steps": result.heating_steps,
        "md_production_steps": result.production_steps,
        "md_total_steps": result.total_steps,
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
