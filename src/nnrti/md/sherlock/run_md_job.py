#!/usr/bin/env python3
"""Run a single MD job from command-line arguments.

Called by SLURM submission scripts (submit_md_batched.sh, submit_all_md.sh)
to avoid shell heredoc quoting issues with inline Python.

Usage:
    python3 -m nnrti.md.sherlock.run_md_job \
        --mutation K103N_P225H --replicate 1 --task-id 0 \
        --system-xml results/md_runs/K103N_P225H/rep_01/assets/K103N_P225H_md_rep01_system.xml \
        --topology-pdb results/md_runs/K103N_P225H/rep_01/assets/K103N_P225H_md_rep01_start.pdb \
        --minimized-pdb results/md_runs/K103N_P225H/rep_01/K103N_P225H_minimized_rep01.pdb \
        --output-json results/md_runs/K103N_P225H/rep_01/K103N_P225H_rep01.json
"""
from __future__ import annotations

import argparse
import logging
import sys

from nnrti.md.manifest import MDTask
from nnrti.md.worker import run_md_task


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    parser = argparse.ArgumentParser(description="Run a single MD job")
    parser.add_argument("--mutation", required=True, help="Mutation label (e.g. K103N_P225H)")
    parser.add_argument("--replicate", type=int, required=True, help="Replicate number (integer, e.g. 1)")
    parser.add_argument("--task-id", type=int, default=0)
    parser.add_argument("--system-xml", required=True, help="Path to prepared system XML")
    parser.add_argument("--topology-pdb", required=True, help="Path to prepared topology PDB")
    parser.add_argument("--minimized-pdb", required=True, help="Path to minimized PDB")
    parser.add_argument("--output-json", required=True, help="Path for output JSON result")
    parser.add_argument("--ligand-sdf", default="data/ligands/dor.sdf")
    parser.add_argument("--ligand-resname", default="2KW")
    parser.add_argument("--heating-ps", type=float, default=25.0)
    parser.add_argument("--production-ns", type=float, default=100.0)
    parser.add_argument("--report-interval", type=int, default=2000)
    parser.add_argument("--checkpoint-interval", type=int, default=5000)
    parser.add_argument("--resume", action="store_true", default=True,
                        help="Resume from checkpoint if available (default: True)")
    parser.add_argument("--no-resume", action="store_true",
                        help="Do not resume from checkpoint")
    parser.add_argument("--force", action="store_true",
                        help="Rerun even if output JSON exists with status=ok")
    args = parser.parse_args(argv)

    safe_label = args.mutation
    # mutation field uses '+' separator for display (e.g. K103N+P225H)
    mutation_display = args.mutation.replace("_", "+")

    task = MDTask(
        task_id=args.task_id,
        structure=safe_label,
        mutation=mutation_display,
        safe_label=safe_label,
        replicate=args.replicate,
        minimized_pdb=args.minimized_pdb,
        ligand_sdf=args.ligand_sdf,
        ligand_resname=args.ligand_resname,
        fold_reduction=None,
        output_json=args.output_json,
        prepared_system_xml=args.system_xml,
        prepared_topology_pdb=args.topology_pdb,
    )

    resume = args.resume and not args.no_resume

    try:
        result = run_md_task(
            task=task,
            heating_ps=args.heating_ps,
            production_ns=args.production_ns,
            report_interval=args.report_interval,
            checkpoint_interval=args.checkpoint_interval,
            resume_from_checkpoint=resume,
            force=args.force,
        )
        logging.info("MD completed with status: %s", result.get("status"))
        return 0
    except Exception as exc:
        logging.error("MD failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
