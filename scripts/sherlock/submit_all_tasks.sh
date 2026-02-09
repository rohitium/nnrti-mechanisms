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

# Optional override. If unset, runtime auto-selects CUDA/OpenCL/CPU.
export OPENMM_PLATFORM="${OPENMM_PLATFORM:-}"
mkdir -p logs
MANIFEST_PATH="${MANIFEST_PATH:-results/fep_manifest.csv}"

if [ "${OPENMM_PLATFORM}" = "CUDA" ] || [ -z "${OPENMM_PLATFORM}" ]; then
  if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo "ERROR: nvidia-smi not found; GPU node/runtime not available." >&2
    exit 1
  fi
  if ! nvidia-smi -L >/dev/null 2>&1; then
    echo "ERROR: GPU not visible to job (nvidia-smi -L failed)." >&2
    exit 1
  fi
fi

python3 -m src.cluster.fep_worker \
  --manifest "${MANIFEST_PATH}" \
  --task-id ${SLURM_ARRAY_TASK_ID} \
  --heating-ps 25 \
  --production-ns 2.0 \
  --report-interval 2000
