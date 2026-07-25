#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
MANIFEST="${MANIFEST:-$PROJECT_ROOT/results/analysis/fep_jorgensen/worker_manifest.csv}"
SHERLOCK_OPENMM_MODULE="${SHERLOCK_OPENMM_MODULE:-chemistry py-openmm/8.1.1_py312}"
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

if [[ -n "${SHERLOCK_ARRAY_TASK:-}" ]]; then
    if (( SHERLOCK_ARRAY_TASK < 0 || SHERLOCK_ARRAY_TASK > last_task )); then
        echo "SHERLOCK_ARRAY_TASK=${SHERLOCK_ARRAY_TASK} out of range 0-${last_task}" >&2
        exit 1
    fi
    ARRAY_SPEC="${SHERLOCK_ARRAY_TASK}"
else
    ARRAY_SPEC="0-${last_task}%${SHERLOCK_MAX_CONCURRENT}"
fi

echo "Manifest:  $MANIFEST"
echo "Array:     $ARRAY_SPEC"
echo "Partition: $SHERLOCK_PARTITION  GRES: $SHERLOCK_GRES  TIME: $SHERLOCK_TIME  MEM: $SHERLOCK_MEM"
echo

sbatch \
    --job-name=fep_jorgensen \
    --partition="$SHERLOCK_PARTITION" \
    --gres="$SHERLOCK_GRES" \
    --time="$SHERLOCK_TIME" \
    --mem="$SHERLOCK_MEM" \
    --array="$ARRAY_SPEC" \
    --output="$PROJECT_ROOT/logs/fep_jorgensen.%A_%a.out" \
    --error="$PROJECT_ROOT/logs/fep_jorgensen.%A_%a.err" \
    <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

module load ${SHERLOCK_OPENMM_MODULE}

cd ${PROJECT_ROOT}

python3 -m scripts.fep_jorgensen.run_manifest_task \
    --manifest ${MANIFEST} \
    --task-id \${SLURM_ARRAY_TASK_ID}
SBATCH_EOF
