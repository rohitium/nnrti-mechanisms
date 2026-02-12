#!/bin/bash
#
# Submit MD jobs in batches with automatic queue monitoring
# Usage: ./submit_md_batched.sh [batch_size] [max_concurrent]
#
# Example: ./submit_md_batched.sh 10 15
#   - Submits 10 jobs at a time
#   - Waits until you have <15 jobs running/pending before submitting next batch
#

set -e

PROJECT_ROOT="/scratch/users/rsatija/nnrti-mechanisms"
cd "$PROJECT_ROOT"

# Configurable parameters
BATCH_SIZE=${1:-10}           # How many to submit per batch (default: 10)
MAX_CONCURRENT=${2:-15}       # Max jobs to have in queue at once (default: 15)
POLL_INTERVAL=300             # Check queue every 5 minutes (300s)

echo "=========================================="
echo "Batched MD Submission"
echo "=========================================="
echo "Batch size:      $BATCH_SIZE jobs"
echo "Max concurrent:  $MAX_CONCURRENT jobs"
echo "Poll interval:   ${POLL_INTERVAL}s"
echo ""

# Find all systems needing MD
SYSTEMS_TO_RUN=()
for SYSTEM_XML in $(find results/md_runs -name "*_system.xml" | sort); do
    DIR=$(dirname "$SYSTEM_XML")
    PARENT=$(dirname "$DIR")
    MUTATION=$(basename "$(dirname "$PARENT")")
    REP=$(basename "$PARENT" | sed 's/rep_//')

    RESULT_JSON="$PARENT/${MUTATION}_rep${REP}.json"

    # Only queue systems without completed results
    if [ ! -f "$RESULT_JSON" ]; then
        SYSTEMS_TO_RUN+=("$MUTATION:$REP:$PARENT:$SYSTEM_XML:$DIR")
    fi
done

TOTAL=${#SYSTEMS_TO_RUN[@]}
echo "Found $TOTAL systems to run"
echo ""

if [ $TOTAL -eq 0 ]; then
    echo "No systems need MD runs. All done!"
    exit 0
fi

# Function to count current jobs in queue
count_my_jobs() {
    # Count jobs with md_* prefix in name
    local count=$(squeue -u $USER -h -t PD,R 2>/dev/null | grep -c "md_" || echo "0")
    echo "$count"
}

# Function to submit a single job
submit_job() {
    local MUTATION=$1
    local REP=$2
    local PARENT=$3
    local SYSTEM_XML=$4
    local DIR=$5
    local TASK_ID=$6
    local REP_INT=$7
    local RESULT_JSON="$PARENT/${MUTATION}_rep${REP}.json"

    sbatch --parsable <<SBATCH_EOF
#!/bin/bash
#SBATCH --job-name=md_${MUTATION}_${REP}
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --output=${PROJECT_ROOT}/logs/md_${MUTATION}_rep${REP}_%j.log

module load chemistry py-openmm/8.1.1_py312

cd ${PROJECT_ROOT}

python3 -m src.md.sherlock.run_md_job \
    --mutation "${MUTATION}" \
    --replicate ${REP_INT} \
    --task-id ${TASK_ID} \
    --system-xml "${SYSTEM_XML}" \
    --topology-pdb "${DIR}/${MUTATION}_md_rep${REP}_start.pdb" \
    --minimized-pdb "${PARENT}/${MUTATION}_minimized_rep${REP}.pdb" \
    --output-json "${RESULT_JSON}" \
    --resume
SBATCH_EOF
}

# Main submission loop
SUBMITTED=0
BATCH_NUM=1
JOBS_IN_BATCH=0

for SYSTEM_INFO in "${SYSTEMS_TO_RUN[@]}"; do
    IFS=':' read -r MUTATION REP PARENT SYSTEM_XML DIR <<< "$SYSTEM_INFO"
    REP_INT=$((10#$REP))  # strip leading zero for Python literal

    # Submit the job
    echo "→ Submitting $MUTATION rep $REP"
    JOBID=$(submit_job "$MUTATION" "$REP" "$PARENT" "$SYSTEM_XML" "$DIR" $SUBMITTED $REP_INT)
    echo "  Job ID: $JOBID"

    SUBMITTED=$((SUBMITTED + 1))
    JOBS_IN_BATCH=$((JOBS_IN_BATCH + 1))

    sleep 0.5  # Rate limit submissions

    # When we've submitted a full batch, pause and wait
    if [ $JOBS_IN_BATCH -ge $BATCH_SIZE ]; then
        echo ""
        echo "✓ Batch $BATCH_NUM complete ($SUBMITTED/$TOTAL submitted)"

        # Wait for squeue to catch up
        sleep 5

        CURRENT_JOBS=$(count_my_jobs)
        echo "  Current queue: $CURRENT_JOBS jobs"

        # Reset batch counter
        JOBS_IN_BATCH=0
        BATCH_NUM=$((BATCH_NUM + 1))

        # If we have more to submit, wait until queue drains below threshold
        if [ $SUBMITTED -lt $TOTAL ]; then
            echo ""
            echo "Waiting for queue to drop below $MAX_CONCURRENT before next batch..."

            while true; do
                CURRENT_JOBS=$(count_my_jobs)
                if [ $CURRENT_JOBS -lt $MAX_CONCURRENT ]; then
                    echo "  Queue at $CURRENT_JOBS jobs - proceeding with next batch"
                    echo ""
                    break
                fi
                echo "[$(date '+%H:%M:%S')] Queue at $CURRENT_JOBS/$MAX_CONCURRENT jobs. Checking again in ${POLL_INTERVAL}s..."
                sleep $POLL_INTERVAL
            done
        fi
    fi
done

echo ""
echo "=========================================="
echo "Submission Complete"
echo "=========================================="
echo "Total submitted: $SUBMITTED/$TOTAL"
echo "Current queue:   $(count_my_jobs) jobs"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  watch -n 30 'squeue -u \$USER'"
echo "  tail -f logs/md_*.log"
echo ""
echo "Check completion:"
echo "  ls results/md_runs/*/rep_*/*.json | wc -l"
echo ""
