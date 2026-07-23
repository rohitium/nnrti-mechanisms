#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
MANIFEST="${MANIFEST:-$PROJECT_ROOT/results/analysis/fep_jorgensen/worker_manifest.csv}"
CONDA_ENV="${CONDA_ENV:-nnrti-openmm}"
SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-24:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-20}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Missing worker manifest: $MANIFEST" >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is required; run this on Sherlock" >&2
    exit 1
fi

task_count=$(( $(wc -l < "$MANIFEST") - 1 ))
if (( task_count < 1 )); then
    echo "No tasks in $MANIFEST" >&2
    exit 1
fi
last_task=$((task_count - 1))
mkdir -p "$PROJECT_ROOT/logs"

sbatch \
    --job-name=fep_jorgensen \
    --partition="$SHERLOCK_PARTITION" \
    --gres="$SHERLOCK_GRES" \
    --time="$SHERLOCK_TIME" \
    --mem="$SHERLOCK_MEM" \
    --array="0-${last_task}%${SHERLOCK_MAX_CONCURRENT}" \
    --output="$PROJECT_ROOT/logs/fep_jorgensen.%A_%a.out" \
    --error="$PROJECT_ROOT/logs/fep_jorgensen.%A_%a.err" \
    --wrap="source ~/.bashrc && conda activate '$CONDA_ENV' && cd '$PROJECT_ROOT' && PYTHONPATH=. python -m scripts.fep_jorgensen.run_manifest_task --manifest '$MANIFEST' --task-id \$SLURM_ARRAY_TASK_ID"
