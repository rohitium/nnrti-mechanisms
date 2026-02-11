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
    structure="K103N_P225H",
    mutation="K103N+P225H",
    safe_label="K103N_P225H",
    replicate=2,
    minimized_pdb="results/md_runs/K103N_P225H/rep_02/K103N_P225H_minimized_rep02.pdb",
    ligand_sdf="data/ligands/dor.sdf",
    ligand_resname="2KW",
    fold_reduction=None,
    output_json="results/md_runs/K103N_P225H/rep_02/K103N_P225H_rep02.json",
    prepared_system_xml="results/md_runs/K103N_P225H/rep_02/assets/K103N_P225H_md_rep02_system.xml",
    prepared_topology_pdb="results/md_runs/K103N_P225H/rep_02/assets/K103N_P225H_md_rep02_start.pdb",
)

print(f"Running MD for {task.mutation} replicate {task.replicate}")
print(f"System XML: {task.prepared_system_xml}")
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
