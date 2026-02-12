#!/bin/bash
#
# Submit MD jobs for all prepared systems (have *_system.xml)
# Bypasses Snakemake - directly submits SLURM jobs
#

set -e

PROJECT_ROOT="/scratch/users/rsatija/nnrti-mechanisms"
cd "$PROJECT_ROOT"

echo "Finding prepared systems without completed MD..."

SUBMITTED=0
SKIPPED=0

# Find all system XML files
for SYSTEM_XML in $(find results/md_runs -name "*_system.xml" | sort); do
    # Extract mutation and replicate from path
    # e.g., results/md_runs/K103N_P225H/rep_02/assets/K103N_P225H_md_rep02_system.xml
    DIR=$(dirname "$SYSTEM_XML")
    PARENT=$(dirname "$DIR")
    MUTATION=$(basename "$(dirname "$PARENT")")
    REP=$(basename "$PARENT" | sed 's/rep_//')
    REP_INT=$((10#$REP))

    RESULT_JSON="$PARENT/${MUTATION}_rep${REP}.json"

    # Skip if result JSON already exists
    if [ -f "$RESULT_JSON" ]; then
        echo "✓ Skip $MUTATION rep $REP (already complete)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "→ Submitting $MUTATION rep $REP"

    # Submit SLURM job
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
    --task-id ${SUBMITTED} \
    --system-xml "${SYSTEM_XML}" \
    --topology-pdb "${DIR}/${MUTATION}_md_rep${REP}_start.pdb" \
    --minimized-pdb "${PARENT}/${MUTATION}_minimized_rep${REP}.pdb" \
    --output-json "${RESULT_JSON}" \
    --resume
SBATCH_EOF

    SUBMITTED=$((SUBMITTED + 1))
    sleep 0.5  # Rate limit
done

echo ""
echo "Summary:"
echo "  Submitted: $SUBMITTED"
echo "  Skipped:   $SKIPPED"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check logs:   tail -f logs/md_*"
