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
# EM smoke fits in 2h; equil (5 ns GPU) needs 4–6h if run in the same session.
SHERLOCK_TIME="${SHERLOCK_TIME:-04:00:00}"
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
echo "Once allocated:"
echo "  cd /scratch/users/rsatija/nnrti-mechanisms-git"
echo "  git pull   # if needed"
echo "  bash scripts/fep_pmx/smoke_neq_em.sh"
exec "${CMD[@]}"
