from __future__ import annotations

from pathlib import Path

from .manifest import load_manifest


SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=nnrti_md
#SBATCH --partition={partition}
#SBATCH --gres=gpu:1
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --array=0-{max_task_id}
#SBATCH --output={log_dir}/md_%A_%a.out
#SBATCH --error={log_dir}/md_%A_%a.err

# Runtime setup
{runtime_setup}

mkdir -p {log_dir}

{python_cmd} -m src.cluster.md_worker \\
    --manifest {manifest_path} \\
    --task-id $SLURM_ARRAY_TASK_ID \\
    --heating-ps {heating_ps} \\
    --production-ns {production_ns} \\
    --report-interval {report_interval}
"""


def generate_slurm_script(
    manifest_path: Path,
    output_script: Path,
    partition: str = "gpu",
    time_limit: str = "6:00:00",
    memory: str = "16G",
    log_dir: str = "logs",
    heating_ps: float = 25.0,
    production_ns: float = 2.0,
    report_interval: int = 2000,
    conda_env: str | None = None,
    use_openmm_module: bool = False,
    **_: object,
) -> Path:
    tasks = load_manifest(manifest_path)
    max_task_id = len(tasks) - 1

    if use_openmm_module:
        runtime_setup = "ml chemistry py-openmm/8.1.1_py312"
        python_cmd = "python3"
    elif conda_env:
        runtime_setup = f"ml miniforge/24.11.0-0\nmamba activate {conda_env}"
        python_cmd = "python"
    else:
        runtime_setup = "# TODO: load OpenMM runtime stack"
        python_cmd = "python3"

    script_content = SLURM_TEMPLATE.format(
        partition=partition,
        time_limit=time_limit,
        memory=memory,
        max_task_id=max_task_id,
        log_dir=log_dir,
        manifest_path=manifest_path,
        heating_ps=heating_ps,
        production_ns=production_ns,
        report_interval=report_interval,
        runtime_setup=runtime_setup,
        python_cmd=python_cmd,
    )

    output_script.parent.mkdir(parents=True, exist_ok=True)
    output_script.write_text(script_content)
    return output_script


def get_task_count(manifest_path: Path) -> int:
    return len(load_manifest(manifest_path))
