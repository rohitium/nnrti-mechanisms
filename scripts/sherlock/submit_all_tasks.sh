#!/bin/bash
#SBATCH --job-name=nnrti_fep
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/fep_%A_%a.out
#SBATCH --error=logs/fep_%A_%a.err
#
# Each array task runs BOTH legs (complex + solvent) for one
# (structure × replicate) combination.
#
# The manifest stores tasks as consecutive pairs:
#   array index 0 → task IDs 0 (complex) + 1 (solvent)
#   array index 1 → task IDs 2 (complex) + 3 (solvent)
#   array index N → task IDs 2N (complex) + 2N+1 (solvent)
#
# With 14 structures (WT + 13 mutations) × 3 replicates = 42 pairs.
# Array range: 0-41.  Adjust if manifest changes.
#SBATCH --array=0-41

# Use Sherlock OpenMM module (no conda needed)
ml chemistry py-openmm/8.1.1_py312

# Force CUDA platform for OpenMM
export OPENMM_PLATFORM=CUDA

# Create log directory if needed
mkdir -p logs

COMPLEX_ID=$((SLURM_ARRAY_TASK_ID * 2))
SOLVENT_ID=$((SLURM_ARRAY_TASK_ID * 2 + 1))

FEP_ARGS="--manifest results/fep_manifest.csv --equil-steps 10000 --prod-steps 25000 --sample-interval 200"

echo "=== Array index ${SLURM_ARRAY_TASK_ID}: tasks ${COMPLEX_ID} (complex) + ${SOLVENT_ID} (solvent) ==="

# Run complex leg first (~2 hrs), then solvent leg (~4 min)
python3 -m src.cluster.fep_worker $FEP_ARGS --task-id $COMPLEX_ID
COMPLEX_RC=$?

python3 -m src.cluster.fep_worker $FEP_ARGS --task-id $SOLVENT_ID
SOLVENT_RC=$?

if [ $COMPLEX_RC -ne 0 ] || [ $SOLVENT_RC -ne 0 ]; then
    echo "ERROR: complex_rc=${COMPLEX_RC}, solvent_rc=${SOLVENT_RC}"
    exit 1
fi
echo "=== Both legs completed successfully ==="
