#!/usr/bin/env python3
"""Test MD execution for a single prepared system."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.cluster.manifest import MDTask
from src.cluster.md_worker import run_md_task

# Test system: K103N_P225H rep 02
task = MDTask(
    task_id=1,
    mutation="K103N+P225H",
    safe_label="K103N_P225H",
    replicate=2,
    system_xml=Path("results/md_runs/K103N_P225H/rep_02/assets/K103N_P225H_md_rep02_system.xml"),
    topology_pdb=Path("results/md_runs/K103N_P225H/rep_02/assets/K103N_P225H_md_rep02_start.pdb"),
    output_json=Path("results/md_runs/K103N_P225H/rep_02/K103N_P225H_rep02.json"),
    ligand_sdf=Path("data/ligands/dor.sdf"),
    ligand_resname="2KW",
    fold_reduction=None,
)

print(f"Running MD for {task.mutation} replicate {task.replicate}")
print(f"System XML: {task.system_xml}")
print(f"Output JSON: {task.output_json}")
print()

result = run_md_task(
    task=task,
    heating_ps=25.0,
    production_ns=2.0,
    report_interval=2000,
    checkpoint_interval=5000,
    resume_from_checkpoint=True,
    force=False,
)

print()
print(f"Status: {result.get('status')}")
print(f"Output saved to: {task.output_json}")
