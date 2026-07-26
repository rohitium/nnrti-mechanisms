#!/bin/bash
#
# Request an interactive Sherlock GPU session for pmx NEQ smoke tests.
#
# Uses partition dev (separate QOS: 2 GPUs, 2 h) — does not consume gpu QOS submit quota.
#
# Usage:
#   bash scripts/fep_pmx/salloc_neq_gpu.sh
#
# Once allocated:
#   bash scripts/fep_pmx/smoke_neq_em.sh
#
set -euo pipefail

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-dev}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-02:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-8}"

CMD=(
    salloc
    --partition="${SHERLOCK_PARTITION}"
    --gres="${SHERLOCK_GRES}"
    --time="${SHERLOCK_TIME}"
    --mem="${SHERLOCK_MEM}"
    --cpus-per-task="${SHERLOCK_CPUS_PER_TASK}"
)

echo "Requesting interactive GPU allocation for pmx NEQ smoke test (partition=${SHERLOCK_PARTITION}):"
echo "  ${CMD[*]}"
echo ""
echo "Once allocated:"
echo "  cd /scratch/users/rsatija/nnrti-mechanisms-git"
echo "  git pull"
echo "  bash scripts/fep_pmx/smoke_neq_em.sh"
echo "  bash scripts/fep_pmx/smoke_neq_task.sh          # equil λ0 holo rep1"
echo "  STAGE=extract bash scripts/fep_pmx/smoke_neq_task.sh"
echo "  STAGE=switch SNAPSHOT=0 bash scripts/fep_pmx/smoke_neq_task.sh"
exec "${CMD[@]}"
