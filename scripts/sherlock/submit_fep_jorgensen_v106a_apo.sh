#!/bin/bash
#
# Submit all 11 lambda windows for WT -> V106A apo FEP on Sherlock GPUs.
#
# Holo windows (manifest task_id 0-10) should already be done before you analyze
# ΔΔG_bind. This script only runs apo tasks 11-21.
#
# Prerequisites:
#   1. Perses apo hybrid prep done locally (--phase apo or --phase all)
#   2. results/analysis/fep_jorgensen/legs/wt_to_V106A/apo/ rsync'd to Sherlock
#   3. worker_manifest_v106a.csv present on Sherlock
#
# Usage (on Sherlock):
#   export PROJECT_ROOT=$PWD
#   ./scripts/sherlock/submit_fep_jorgensen_v106a_apo.sh
#
# Optional env vars:
#   SHERLOCK_OPENMM_MODULE="chemistry py-openmm/8.1.1_py312"
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#   SHERLOCK_TIME=24:00:00
#   SHERLOCK_MEM=32G
#   SHERLOCK_MAX_CONCURRENT=11
#   FEP_APO_TASK_START=11   # default first apo manifest task_id
#   FEP_APO_TASK_END=21     # default last apo manifest task_id
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
MANIFEST="${MANIFEST:-$PROJECT_ROOT/results/analysis/fep_jorgensen/worker_manifest_v106a.csv}"
LEG_DIR="$PROJECT_ROOT/results/analysis/fep_jorgensen/legs/wt_to_V106A/apo"
FEP_APO_TASK_START="${FEP_APO_TASK_START:-11}"
FEP_APO_TASK_END="${FEP_APO_TASK_END:-21}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Missing manifest: $MANIFEST" >&2
    echo "Generate locally:" >&2
    echo "  PYTHONPATH=. python -m scripts.fep_jorgensen.panel --mutation V106A" >&2
    exit 1
fi

for required in "$LEG_DIR/hybrid_system.xml" "$LEG_DIR/hybrid_topology.pdb" "$LEG_DIR/schedule.json"; do
    if [[ ! -f "$required" ]]; then
        echo "Missing prepared apo artifact: $required" >&2
        echo "Rsync from laptop:" >&2
        echo "  SHERLOCK_USER=<user> bash scripts/rsync_fep_jorgensen.sh push V106A" >&2
        exit 1
    fi
done

export MANIFEST
export SHERLOCK_ARRAY_TASK="${FEP_APO_TASK_START}-${FEP_APO_TASK_END}%${SHERLOCK_MAX_CONCURRENT:-11}"
export SHERLOCK_TIME="${SHERLOCK_TIME:-24:00:00}"
export SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"

echo "Submitting V106A apo windows (manifest task_id ${FEP_APO_TASK_START}-${FEP_APO_TASK_END})."
exec "$SCRIPT_DIR/submit_fep_jorgensen_windows.sh"
