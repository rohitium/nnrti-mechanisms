#!/bin/bash
#
# Submit MD jobs for systems that are already prepped (have system.xml files)
# This script finds all prepared systems and submits SLURM jobs only for MD execution
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

echo "Finding prepared systems..."
SYSTEMS=$(find results/md_runs -name "*_system.xml" | sort)

if [ -z "$SYSTEMS" ]; then
    echo "No prepared systems found!"
    exit 1
fi

TOTAL=$(echo "$SYSTEMS" | wc -l | tr -d ' ')
echo "Found $TOTAL prepared systems"

SUBMITTED=0
SKIPPED=0

for SYSTEM_XML in $SYSTEMS; do
    # Extract mutation and replicate from path
    # e.g., results/md_runs/Y188L/rep_01/assets/Y188L_md_rep01_system.xml
    DIR=$(dirname "$SYSTEM_XML")
    PARENT=$(dirname "$DIR")
    MUTATION=$(basename "$(dirname "$PARENT")")
    REP=$(basename "$PARENT" | sed 's/rep_//')

    RESULT_JSON="$PARENT/${MUTATION}_rep${REP}.json"

    # Skip if JSON result already exists
    if [ -f "$RESULT_JSON" ]; then
        echo "✓ Skip $MUTATION rep $REP (already complete)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "→ Submitting $MUTATION rep $REP"

    # Submit SLURM job
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=md_${MUTATION}_${REP}
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/md_${MUTATION}_rep${REP}_%j.log

module load chemistry
module load py-openmm/8.1.1_py312

cd $PROJECT_ROOT

python3 -m src.cluster.md_worker \
  --system-xml "$SYSTEM_XML" \
  --output-dir "$PARENT" \
  --mutation "$MUTATION" \
  --replicate "$REP"
EOF

    SUBMITTED=$((SUBMITTED + 1))
    sleep 1  # Rate limit submissions
done

echo ""
echo "Summary:"
echo "  Total systems: $TOTAL"
echo "  Submitted:     $SUBMITTED"
echo "  Skipped:       $SKIPPED"
