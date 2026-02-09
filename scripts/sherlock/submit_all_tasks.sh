#!/bin/bash
#SBATCH --job-name=nnrti_md
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/md_%A_%a.out
#SBATCH --error=logs/md_%A_%a.err
#SBATCH --array=0-0

set -euo pipefail

source /etc/profile.d/modules.sh 2>/dev/null || true
ml chemistry py-openmm/8.1.1_py312

export OPENMM_PLATFORM=CUDA
mkdir -p logs
MANIFEST_PATH="${MANIFEST_PATH:-results/fep_manifest.csv}"

python3 -m src.cluster.fep_worker \
  --manifest "${MANIFEST_PATH}" \
  --task-id ${SLURM_ARRAY_TASK_ID} \
  --heating-ps 25 \
  --production-ns 2.0 \
  --report-interval 2000
