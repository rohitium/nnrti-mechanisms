#!/bin/bash
#SBATCH --job-name=nnrti_fep
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --array=0-83
#SBATCH --output=logs/fep_%A_%a.out
#SBATCH --error=logs/fep_%A_%a.err

# Load miniforge and activate conda environment
ml miniforge/24.11.0-0
mamba activate nnrti

# Force CUDA platform for OpenMM
export OPENMM_PLATFORM=CUDA

# Create log directory if needed
mkdir -p logs

# Run FEP worker for this array task
python -m src.cluster.fep_worker \
    --manifest results/fep_manifest.csv \
    --task-id $SLURM_ARRAY_TASK_ID \
    --equil-steps 10000 \
    --prod-steps 25000 \
    --sample-interval 200
