#!/bin/bash
#
# Submit ONE full-production lambda window for WT -> V106A on Sherlock.
# Use this after the short interactive GPU pilot, before the full 11-window array.
#
# Default task is manifest row 1 (state_index=1, λ≈0.1) to avoid reusing the
# short pilot artifacts on state 0.
#
# Usage (on Sherlock):
#   export PROJECT_ROOT=$PWD
#   ./scripts/sherlock/submit_fep_jorgensen_v106a_single.sh
#
# Optional env vars:
#   FEP_TASK_ID=1              manifest task_id / SLURM array index (default 1)
#   SHERLOCK_TIME=24:00:00
#   SHERLOCK_MEM=32G
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#
# Success check after the job finishes:
#   wc -l results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows/state_01_energies.csv
#   # expect 1001 lines (header + 1000 samples)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
MANIFEST="${MANIFEST:-$PROJECT_ROOT/results/analysis/fep_jorgensen/worker_manifest_v106a.csv}"
LEG_DIR="$PROJECT_ROOT/results/analysis/fep_jorgensen/legs/wt_to_V106A/holo"
FEP_TASK_ID="${FEP_TASK_ID:-1}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Missing manifest: $MANIFEST" >&2
    exit 1
fi

for required in "$LEG_DIR/hybrid_system.xml" "$LEG_DIR/hybrid_topology.pdb" "$LEG_DIR/schedule.json"; do
    if [[ ! -f "$required" ]]; then
        echo "Missing prepared holo artifact: $required" >&2
        exit 1
    fi
done

export MANIFEST
export SHERLOCK_ARRAY_TASK="$FEP_TASK_ID"
export SHERLOCK_MAX_CONCURRENT=1
export SHERLOCK_TIME="${SHERLOCK_TIME:-24:00:00}"
export SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"

echo "Submitting one full-production V106A window (manifest task_id=${FEP_TASK_ID})."
exec "$SCRIPT_DIR/submit_fep_jorgensen_windows.sh"
