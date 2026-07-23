#!/bin/bash
#
# Request an interactive Sherlock GPU session for FEP lambda-window testing.
#
# Usage (from your laptop):
#   ssh sherlock
#   cd $SCRATCH/nnrti-mechanisms
#   bash scripts/sherlock/salloc_fep_jorgensen_gpu.sh
#
# Optional env vars:
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#   SHERLOCK_TIME=04:00:00
#   SHERLOCK_MEM=32G
#   SHERLOCK_CPUS_PER_TASK=4
#   SHERLOCK_QOS=<qos-name>
#

set -euo pipefail

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
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

if [ -n "${SHERLOCK_QOS}" ]; then
    CMD+=(--qos="${SHERLOCK_QOS}")
fi

echo "Requesting interactive GPU allocation for FEP window tests:"
echo "  ${CMD[*]}"
echo
echo "After the shell starts:"
echo "  cd \"\${SCRATCH:-\$HOME}/nnrti-mechanisms-git\"   # or your repo path"
echo "  export PROJECT_ROOT=\$PWD"
echo "  bash scripts/sherlock/run_fep_jorgensen_pilot.sh   # loads py-openmm module internally"
exec "${CMD[@]}"
