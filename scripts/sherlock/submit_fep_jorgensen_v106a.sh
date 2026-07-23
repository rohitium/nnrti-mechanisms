#!/bin/bash
#
# Submit all lambda windows for WT -> V106A holo FEP on Sherlock GPUs.
#
# Prerequisites:
#   1. Perses hybrid prep done locally for V106A
#   2. results/analysis/fep_jorgensen/legs/wt_to_V106A/ rsync'd to Sherlock
#   3. worker_manifest_v106a.csv present (generate locally with panel.py --mutation V106A)
#
# Usage (on Sherlock):
#   ./scripts/sherlock/submit_fep_jorgensen_v106a.sh
#
# Optional env vars (passed through to submit_fep_jorgensen_windows.sh):
#   CONDA_ENV=nnrti-openmm
#   SHERLOCK_PARTITION=gpu
#   SHERLOCK_GRES=gpu:1
#   SHERLOCK_TIME=24:00:00
#   SHERLOCK_MEM=32G
#   SHERLOCK_MAX_CONCURRENT=11
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
MANIFEST="${MANIFEST:-$PROJECT_ROOT/results/analysis/fep_jorgensen/worker_manifest_v106a.csv}"
LEG_DIR="$PROJECT_ROOT/results/analysis/fep_jorgensen/legs/wt_to_V106A/holo"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Missing manifest: $MANIFEST" >&2
    echo "Generate locally:" >&2
    echo "  PYTHONPATH=. python -m scripts.fep_jorgensen.panel --mutation V106A" >&2
    exit 1
fi

for required in "$LEG_DIR/hybrid_system.xml" "$LEG_DIR/hybrid_topology.pdb" "$LEG_DIR/schedule.json"; do
    if [[ ! -f "$required" ]]; then
        echo "Missing prepared holo artifact: $required" >&2
        echo "Rsync from laptop:" >&2
        echo "  rsync -av results/analysis/fep_jorgensen/legs/wt_to_V106A/ \\" >&2
        echo "    sherlock:\$SCRATCH/nnrti-mechanisms/results/analysis/fep_jorgensen/legs/wt_to_V106A/" >&2
        exit 1
    fi
done

export MANIFEST
export SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-11}"
exec "$SCRIPT_DIR/submit_fep_jorgensen_windows.sh"
