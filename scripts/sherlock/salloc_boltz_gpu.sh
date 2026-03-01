#!/bin/bash
#
# Request an interactive Sherlock GPU session suitable for Boltz-2 runs.
#
# Usage:
#   bash scripts/sherlock/salloc_boltz_gpu.sh
#
# Optional env vars:
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#   SHERLOCK_TIME=04:00:00
#   SHERLOCK_MEM=64G
#   SHERLOCK_CPUS_PER_TASK=8
#   SHERLOCK_QOS=<qos-name>
#

set -euo pipefail

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-04:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-64G}"
SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-8}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"

CMD=(
    salloc
    --partition="${SHERLOCK_PARTITION}"
    --gres="${SHERLOCK_GRES}"
    --time="${SHERLOCK_TIME}"
    --mem="${SHERLOCK_MEM}"
    --cpus-per-task="${SHERLOCK_CPUS_PER_TASK}"
)

if [ -n "${SHERLOCK_QOS}" ]; then
    CMD+=(--qos="${SHERLOCK_QOS}")
fi

echo "Requesting interactive allocation:"
echo "  ${CMD[*]}"
exec "${CMD[@]}"

