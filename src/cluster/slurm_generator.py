from __future__ import annotations

from pathlib import Path

from .manifest import load_manifest


SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=nnrti_fep
#SBATCH --partition={partition}
#SBATCH --gres=gpu:1
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --array=0-{max_task_id}
#SBATCH --output={log_dir}/fep_%A_%a.out
#SBATCH --error={log_dir}/fep_%A_%a.err

# Load modules
module load python/3.11
module load cuda/12.2

# Activate conda environment
source activate nnrti

# Force CUDA platform for OpenMM
export OPENMM_PLATFORM=CUDA

# Run FEP worker for this array task
python -m src.cluster.fep_worker \\
    --manifest {manifest_path} \\
    --task-id $SLURM_ARRAY_TASK_ID \\
    --equil-steps {equil_steps} \\
    --prod-steps {prod_steps} \\
    --sample-interval {sample_interval}
"""


def generate_slurm_script(
    manifest_path: Path,
    output_script: Path,
    partition: str = "gpu",
    time_limit: str = "4:00:00",
    memory: str = "16G",
    log_dir: str = "logs",
    equil_steps: int = 10_000,
    prod_steps: int = 25_000,
    sample_interval: int = 200,
) -> Path:
    """Generate a SLURM array job submission script.

    Args:
        manifest_path: Path to the FEP manifest CSV.
        output_script: Path to write the generated SLURM script.
        partition: SLURM partition to use.
        time_limit: Job time limit.
        memory: Memory allocation per task.
        log_dir: Directory for log files (relative to job submission directory).
        equil_steps: Equilibration steps per lambda window.
        prod_steps: Production steps per lambda window.
        sample_interval: Sample interval for energy evaluations.

    Returns:
        Path to the generated script.
    """
    tasks = load_manifest(manifest_path)
    max_task_id = len(tasks) - 1

    script_content = SLURM_TEMPLATE.format(
        partition=partition,
        time_limit=time_limit,
        memory=memory,
        max_task_id=max_task_id,
        log_dir=log_dir,
        manifest_path=manifest_path,
        equil_steps=equil_steps,
        prod_steps=prod_steps,
        sample_interval=sample_interval,
    )

    output_script.parent.mkdir(parents=True, exist_ok=True)
    output_script.write_text(script_content)

    return output_script


def get_task_count(manifest_path: Path) -> int:
    """Return the total number of tasks in a manifest."""
    return len(load_manifest(manifest_path))
