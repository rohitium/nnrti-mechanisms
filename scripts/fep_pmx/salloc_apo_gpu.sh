#!/bin/bash
#
# Request an interactive Sherlock GPU session for apo MD smoke tests.
#
# Usage:
#   bash scripts/fep_pmx/salloc_apo_gpu.sh
#
# Optional:
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#   SHERLOCK_TIME=01:00:00
#   SHERLOCK_MEM=32G

set -euo pipefail

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-01:00:00}"
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

echo "Requesting interactive GPU allocation for apo MD smoke test:"
echo "  ${CMD[*]}"
echo ""
echo "Once allocated, run:"
echo "  cd /scratch/users/rsatija/nnrti-mechanisms-git"
echo "  bash scripts/fep_pmx/test_y188l_apo_gpu.sh"
exec "${CMD[@]}"
