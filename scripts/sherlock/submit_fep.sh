#!/bin/bash
#SBATCH --job-name=nnrti_fep
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --array=0-83
#SBATCH --output=logs/fep_%A_%a.out
#SBATCH --error=logs/fep_%A_%a.err

# NNRTI FEP Array Job for Sherlock GPU Cluster
#
# This script runs alchemical FEP calculations for the DOR pipeline.
# Each array task processes one (mutation, replicate, leg) combination.
#
# Total tasks: 14 structures × 3 replicates × 2 legs = 84 tasks
#
# Usage:
#   1. Transfer project to Sherlock:
#      rsync -avz --exclude='.venv' --exclude='.git' \
#          . sherlock:/scratch/users/$USER/nnrti-mechanisms/
#
#   2. Submit job:
#      sbatch scripts/sherlock/submit_fep.sh
#
#   3. Monitor progress:
#      squeue -u $USER
#      tail -f logs/fep_*.out
#
#   4. Transfer results back:
#      rsync -avz sherlock:/scratch/users/$USER/nnrti-mechanisms/results/fep_runs/ \
#          results/fep_runs/

# Exit on error
set -e

# Print job info
echo "=========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on: $(hostname)"
echo "Started at: $(date)"
echo "=========================================="

# Load required modules
module load python/3.11
module load cuda/12.2

# Activate conda environment
# Adjust this path as needed for your Sherlock setup
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate nnrti

# Force CUDA platform for OpenMM
export OPENMM_PLATFORM=CUDA

# Create log directory if it doesn't exist
mkdir -p logs

# Run FEP worker for this array task
python -m src.cluster.fep_worker \
    --manifest results/fep_manifest.csv \
    --task-id $SLURM_ARRAY_TASK_ID \
    --equil-steps 10000 \
    --prod-steps 25000 \
    --sample-interval 200

echo "=========================================="
echo "Finished at: $(date)"
echo "=========================================="
