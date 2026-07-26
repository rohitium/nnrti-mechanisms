#!/bin/bash
#
# Request an interactive Sherlock GPU session for pmx NEQ smoke tests.
#
# Usage:
#   bash scripts/fep_pmx/salloc_neq_gpu.sh
#
# Once allocated:
#   bash scripts/fep_pmx/smoke_neq_em.sh
#
set -euo pipefail

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
# Prefer dev QOS for interactive smoke tests (not production batch queue).
SHERLOCK_QOS="${SHERLOCK_QOS:-dev}"
# EM smoke fits in 2h; full equil+extract+switch smoke needs 6–8h.
SHERLOCK_TIME="${SHERLOCK_TIME:-08:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-4}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"

CMD=(
    salloc
    --partition="${SHERLOCK_PARTITION}"
    --gres="${SHERLOCK_GRES}"
    --time="${SHERLOCK_TIME}"
    --mem="${SHERLOCK_MEM}"
    --cpus-per-task="${SHERLOCK_CPUS_PER_TASK}"
)

if [[ -n "${SHERLOCK_QOS}" ]]; then
    CMD+=(--qos="${SHERLOCK_QOS}")
fi

echo "Requesting interactive GPU allocation for pmx NEQ smoke test:"
echo "  ${CMD[*]}"
echo ""
echo "Once allocated (dev GPU — smoke test before any sbatch):"
echo "  cd /scratch/users/rsatija/nnrti-mechanisms-git"
echo "  git pull"
echo "  bash scripts/fep_pmx/smoke_neq_em.sh"
echo "  bash scripts/fep_pmx/smoke_neq_task.sh          # equil λ0 holo rep1"
echo "  STAGE=extract bash scripts/fep_pmx/smoke_neq_task.sh"
echo "  STAGE=switch SNAPSHOT=0 bash scripts/fep_pmx/smoke_neq_task.sh"
exec "${CMD[@]}"
