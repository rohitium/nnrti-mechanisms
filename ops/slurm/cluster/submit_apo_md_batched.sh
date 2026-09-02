#!/bin/bash
#
# Submit apo (ligand-free) MD jobs in batches with automatic queue monitoring.
#
# Runs identically to submit_md_batched.sh but targets results/md_runs/apo/
# and uses the _apo_ file naming convention for topology/system assets.
#
# Usage:
#   ./ops/slurm/cluster/submit_apo_md_batched.sh [batch_size] [max_concurrent]
#
# Extension reruns:
#   MD_PRODUCTION_NS=100.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 \
#   ./ops/slurm/cluster/submit_apo_md_batched.sh 6 12
#

set -euo pipefail

SLEEP_PID=""
cleanup_and_exit() {
    echo ""
    echo "Received interrupt/termination signal. Stopping submit loop."
    if [ -n "${SLEEP_PID}" ]; then
        kill "${SLEEP_PID}" 2>/dev/null || true
    fi
    exit 130
}
trap cleanup_and_exit INT TERM

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

LOG_DIR="${PROJECT_ROOT}/logs/apo_md"
mkdir -p "${LOG_DIR}"

BATCH_SIZE="${1:-${BATCH_SIZE:-10}}"
MAX_CONCURRENT="${2:-${MAX_CONCURRENT:-15}}"
POLL_INTERVAL="${POLL_INTERVAL:-300}"

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-12:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"

MD_HEATING_PS="${MD_HEATING_PS:-25}"
MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-100.0}"
MD_REPORT_INTERVAL="${MD_REPORT_INTERVAL:-2000}"
MD_CHECKPOINT_INTERVAL="${MD_CHECKPOINT_INTERVAL:-5000}"
MD_RESUME_FROM_CHECKPOINT="${MD_RESUME_FROM_CHECKPOINT:-1}"
MD_FORCE_RERUN="${MD_FORCE_RERUN:-1}"
SKIP_IF_AT_TARGET="${SKIP_IF_AT_TARGET:-1}"
SKIP_IF_RUNNING="${SKIP_IF_RUNNING:-1}"
MUTATION_ALLOWLIST="${MUTATION_ALLOWLIST:-}"

APO_RUNS_ROOT="${APO_RUNS_ROOT:-results/md_runs/apo}"

# Login-node default python3 may be 3.6; reconcile_md_metadata needs >=3.7
if command -v module >/dev/null 2>&1; then
    module load python/3.9.0 2>/dev/null || module load python/3.12.1 2>/dev/null || true
fi
PYTHON="${PYTHON:-python3}"

if ! command -v jq >/dev/null 2>&1; then
    echo "ERROR: jq is required for JSON status checks." >&2
    exit 1
fi

TARGET_STEPS="$("${PYTHON}" - <<PY
ns = float("${MD_PRODUCTION_NS}")
print(max(1, int(round((ns * 1_000_000.0) / 2.0))))
PY
)"

echo "=========================================="
echo "Batched Apo MD Submission"
echo "=========================================="
echo "Apo runs root:   $APO_RUNS_ROOT"
echo "Batch size:      $BATCH_SIZE jobs"
echo "Max concurrent:  $MAX_CONCURRENT jobs"
echo "Poll interval:   ${POLL_INTERVAL}s"
echo "Partition:       $SHERLOCK_PARTITION"
echo "GRES:            $SHERLOCK_GRES"
echo "Time:            $SHERLOCK_TIME"
echo "Mem:             $SHERLOCK_MEM"
echo "Production ns:   $MD_PRODUCTION_NS"
echo "Target steps:    $TARGET_STEPS"
echo "Force rerun:     $MD_FORCE_RERUN"
echo "Skip at target:  $SKIP_IF_AT_TARGET"
if [ -n "${MUTATION_ALLOWLIST}" ]; then
    echo "Mutation allowlist: ${MUTATION_ALLOWLIST}"
fi
echo ""

# Read-only by design. This used to pass --write, which rewrote 25 run JSONs from
# stale state.csv files on 2026-08-17 (halving md_production_steps_completed and so
# every analysis time axis derived from it). A submission must never mutate run
# metadata as a side effect; run reconcile deliberately if you want it to write.
"${PYTHON}" ops/slurm/cluster/reconcile_md_metadata.py \
    --root . \
    --include-apo \
    --target-ns "${MD_PRODUCTION_NS}" >/dev/null

SYSTEMS_TO_RUN=()
SKIPPED_DONE=0
SKIPPED_AT_TARGET=0
SKIPPED_MISSING=0
SKIPPED_RUNNING=0
SKIPPED_FILTERED=0
TOTAL_PREPARED=0

declare -A ALLOWED_MUTATIONS=()
if [ -n "${MUTATION_ALLOWLIST}" ]; then
    IFS=',' read -r -a _allow_tokens <<< "${MUTATION_ALLOWLIST}"
    for tok in "${_allow_tokens[@]}"; do
        t="$(echo "$tok" | tr -d '[:space:]')"
        if [ -n "$t" ]; then
            ALLOWED_MUTATIONS["$t"]=1
        fi
    done
fi

declare -A ACTIVE_JOB_NAMES=()
if [ "$SKIP_IF_RUNNING" = "1" ]; then
    while IFS= read -r job_name; do
        if [[ "$job_name" == apo_* ]]; then
            ACTIVE_JOB_NAMES["$job_name"]=1
        fi
    done < <(squeue -u "$USER" -h -t PD,R -o "%j" 2>/dev/null || true)
fi

# Apo system XMLs are named *_apo_md_rep*_system.xml
while IFS= read -r -d '' SYSTEM_XML; do
    TOTAL_PREPARED=$((TOTAL_PREPARED + 1))
    DIR="$(dirname "$SYSTEM_XML")"          # .../rep_NN/assets
    PARENT="$(dirname "$DIR")"              # .../rep_NN
    MUTATION="$(basename "$(dirname "$PARENT")")"  # mutation label (safe_label)
    REP="$(basename "$PARENT" | sed 's/^rep_//')"  # e.g. "01"
    REP_INT=$((10#$REP))

    # Apo topology PDB uses _apo_md_ infix.
    TOPOLOGY_PDB="${DIR}/${MUTATION}_apo_md_rep${REP}_start.pdb"
    RESULT_JSON="${PARENT}/${MUTATION}_apo_rep${REP}.json"
    JOB_NAME="apo_${MUTATION}_${REP}"

    if [ -n "${MUTATION_ALLOWLIST}" ] && [ -z "${ALLOWED_MUTATIONS[$MUTATION]+x}" ]; then
        SKIPPED_FILTERED=$((SKIPPED_FILTERED + 1))
        continue
    fi

    if [ ! -f "$TOPOLOGY_PDB" ]; then
        echo "⚠ Skip $MUTATION rep $REP (missing apo topology PDB)"
        SKIPPED_MISSING=$((SKIPPED_MISSING + 1))
        continue
    fi

    STATUS=""
    COMPLETED_STEPS=0
    if [ -f "$RESULT_JSON" ]; then
        STATUS="$(jq -r '.status // ""' "$RESULT_JSON" 2>/dev/null || echo "")"
        COMPLETED_STEPS="$(jq -r '.md_production_steps_completed // .md_production_steps // 0' "$RESULT_JSON" 2>/dev/null || echo "0")"
        if ! [[ "$COMPLETED_STEPS" =~ ^[0-9]+$ ]]; then
            COMPLETED_STEPS=0
        fi
    fi

    if [ "$MD_FORCE_RERUN" != "1" ]; then
        if [ "$STATUS" = "ok" ]; then
            SKIPPED_DONE=$((SKIPPED_DONE + 1))
            continue
        fi
    else
        if [ "$SKIP_IF_AT_TARGET" = "1" ] && [ "$STATUS" = "ok" ] && [ "$COMPLETED_STEPS" -ge "$TARGET_STEPS" ]; then
            SKIPPED_AT_TARGET=$((SKIPPED_AT_TARGET + 1))
            continue
        fi
    fi

    if [ "$SKIP_IF_RUNNING" = "1" ] && [ -n "${ACTIVE_JOB_NAMES[$JOB_NAME]+x}" ]; then
        SKIPPED_RUNNING=$((SKIPPED_RUNNING + 1))
        continue
    fi

    SYSTEMS_TO_RUN+=("${MUTATION}:${REP}:${PARENT}:${SYSTEM_XML}:${DIR}:${REP_INT}:${RESULT_JSON}")
done < <(find "${APO_RUNS_ROOT}" -name "*_apo_md_rep*_system.xml" -print0 | sort -z)

TOTAL=${#SYSTEMS_TO_RUN[@]}
echo "Prepared systems:      $TOTAL_PREPARED"
echo "Will submit:           $TOTAL"
echo "Skipped done:          $SKIPPED_DONE"
echo "Skipped at target:     $SKIPPED_AT_TARGET"
echo "Skipped running:       $SKIPPED_RUNNING"
echo "Skipped by filter:     $SKIPPED_FILTERED"
echo "Skipped missing input: $SKIPPED_MISSING"
echo "Found $TOTAL systems to run"
echo ""

if [ "$TOTAL" -eq 0 ]; then
    echo "No apo systems need MD runs. All done!"
    exit 0
fi

count_my_jobs() {
    squeue -u "$USER" -h -t PD,R 2>/dev/null | grep -c "apo_" || true
}

submit_job() {
    local MUTATION="$1"
    local REP="$2"
    local PARENT="$3"
    local SYSTEM_XML="$4"
    local DIR="$5"
    local TASK_ID="$6"
    local REP_INT="$7"
    local RESULT_JSON="$8"

    local RESUME_FLAG="--resume"
    local FORCE_FLAG=""
    if [ "$MD_RESUME_FROM_CHECKPOINT" != "1" ]; then
        RESUME_FLAG="--no-resume"
    fi
    if [ "$MD_FORCE_RERUN" = "1" ]; then
        FORCE_FLAG="--force"
    fi
    local EXTRA_FLAGS="${RESUME_FLAG}"
    if [ -n "${FORCE_FLAG}" ]; then
        EXTRA_FLAGS="${EXTRA_FLAGS} ${FORCE_FLAG}"
    fi

    local SBATCH_ARGS=(
        --parsable
        --job-name="apo_${MUTATION}_${REP}"
        --partition="${SHERLOCK_PARTITION}"
        --time="${SHERLOCK_TIME}"
        --mem="${SHERLOCK_MEM}"
        --output="${LOG_DIR}/apo_${MUTATION}_rep${REP}_%j.log"
    )
    if [ -n "${SHERLOCK_GRES}" ]; then
        SBATCH_ARGS+=(--gres="${SHERLOCK_GRES}")
    fi
    if [ -n "${SHERLOCK_QOS}" ]; then
        SBATCH_ARGS+=(--qos="${SHERLOCK_QOS}")
    fi

    sbatch "${SBATCH_ARGS[@]}" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

module load chemistry py-openmm/8.1.1_py312

cd ${PROJECT_ROOT}
export PYTHONPATH="${PROJECT_ROOT}/src:\${PYTHONPATH:-}"

python3 -m nnrti.md.sherlock.run_md_job \
    --mutation "${MUTATION}" \
    --replicate ${REP_INT} \
    --task-id ${TASK_ID} \
    --system-xml "${SYSTEM_XML}" \
    --topology-pdb "${DIR}/${MUTATION}_apo_md_rep${REP}_start.pdb" \
    --minimized-pdb "${PARENT}/${MUTATION}_minimized_rep${REP}.pdb" \
    --output-json "${RESULT_JSON}" \
    --ligand-sdf "" \
    --ligand-resname "" \
    --heating-ps "${MD_HEATING_PS}" \
    --production-ns "${MD_PRODUCTION_NS}" \
    --report-interval "${MD_REPORT_INTERVAL}" \
    --checkpoint-interval "${MD_CHECKPOINT_INTERVAL}" \
    ${EXTRA_FLAGS}
SBATCH_EOF
}

SUBMITTED=0
BATCH_NUM=1
JOBS_IN_BATCH=0

for SYSTEM_INFO in "${SYSTEMS_TO_RUN[@]}"; do
    IFS=':' read -r MUTATION REP PARENT SYSTEM_XML DIR REP_INT RESULT_JSON <<< "$SYSTEM_INFO"

    echo "→ Submitting apo $MUTATION rep $REP"
    JOBID=$(submit_job "$MUTATION" "$REP" "$PARENT" "$SYSTEM_XML" "$DIR" "$SUBMITTED" "$REP_INT" "$RESULT_JSON")
    echo "  Job ID: $JOBID"

    SUBMITTED=$((SUBMITTED + 1))
    JOBS_IN_BATCH=$((JOBS_IN_BATCH + 1))

    sleep 0.5

    if [ "$JOBS_IN_BATCH" -ge "$BATCH_SIZE" ]; then
        echo ""
        echo "✓ Batch $BATCH_NUM complete ($SUBMITTED/$TOTAL submitted)"
        sleep 5
        CURRENT_JOBS=$(count_my_jobs)
        echo "  Current queue: $CURRENT_JOBS apo jobs"
        JOBS_IN_BATCH=0
        BATCH_NUM=$((BATCH_NUM + 1))

        if [ "$SUBMITTED" -lt "$TOTAL" ]; then
            echo "Waiting for queue to drop below $MAX_CONCURRENT before next batch..."
            while true; do
                CURRENT_JOBS=$(count_my_jobs)
                if [ "$CURRENT_JOBS" -lt "$MAX_CONCURRENT" ]; then
                    echo "  Queue at $CURRENT_JOBS jobs - proceeding with next batch"
                    echo ""
                    break
                fi
                echo "[$(date '+%H:%M:%S')] Queue at $CURRENT_JOBS/$MAX_CONCURRENT jobs. Checking again in ${POLL_INTERVAL}s..."
                sleep "$POLL_INTERVAL" &
                SLEEP_PID=$!
                wait "$SLEEP_PID" || true
                SLEEP_PID=""
            done
        fi
    fi
done

echo ""
echo "=========================================="
echo "Apo Submission Complete"
echo "=========================================="
echo "Total submitted: $SUBMITTED/$TOTAL"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  jq -r '[.safe_label,.replicate,.status,(.md_production_steps_completed // 0)] | @tsv' \\"
echo "    ${APO_RUNS_ROOT}/*/rep_*/*_apo_rep[0-9][0-9].json | awk '\$3!=\"ok\" || \$4<${TARGET_STEPS} {print}'"
echo ""
