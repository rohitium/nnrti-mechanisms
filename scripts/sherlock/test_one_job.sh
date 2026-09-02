#!/bin/bash
#
# Test a single MD job interactively on a GPU dev node.
# Runs a very short simulation (0.01 ns) to verify the pipeline works.
#
# Usage:
#   salloc -p gpu --gres=gpu:1 --time=1:00:00 --mem=32G
#   ./scripts/sherlock/test_one_job.sh [mutation] [rep]
#
# Examples:
#   ./scripts/sherlock/test_one_job.sh              # auto-picks first incomplete system
#   ./scripts/sherlock/test_one_job.sh G190E 01     # test specific system
#

set -e

PROJECT_ROOT="/scratch/users/rsatija/nnrti-mechanisms"
cd "$PROJECT_ROOT"
export PYTHONPATH="${PROJECT_ROOT:-$PWD}/src:${PYTHONPATH:-}"

module load chemistry py-openmm/8.1.1_py312

MUTATION=${1:-}
REP=${2:-}

if [ -z "$MUTATION" ]; then
    echo "No mutation specified — auto-picking first incomplete system..."
    for SYSTEM_XML in $(find results/md_runs -name "*_system.xml" | sort); do
        DIR=$(dirname "$SYSTEM_XML")
        PARENT=$(dirname "$DIR")
        M=$(basename "$(dirname "$PARENT")")
        R=$(basename "$PARENT" | sed 's/rep_//')
        RESULT_JSON="$PARENT/${M}_rep${R}.json"
        if [ ! -f "$RESULT_JSON" ]; then
            MUTATION="$M"
            REP="$R"
            break
        fi
    done
    if [ -z "$MUTATION" ]; then
        echo "All systems already have results!"
        exit 0
    fi
    echo "Selected: $MUTATION rep $REP"
fi

REP_INT=$((10#$REP))

# Derive paths
PARENT="results/md_runs/${MUTATION}/rep_${REP}"
DIR="${PARENT}/assets"
SYSTEM_XML="${DIR}/${MUTATION}_md_rep${REP}_system.xml"
TOPOLOGY_PDB="${DIR}/${MUTATION}_md_rep${REP}_start.pdb"
MINIMIZED_PDB="${PARENT}/${MUTATION}_minimized_rep${REP}.pdb"
OUTPUT_JSON="${PARENT}/${MUTATION}_rep${REP}_TEST.json"

echo "=========================================="
echo "Test MD Job"
echo "=========================================="
echo "Mutation:     $MUTATION"
echo "Replicate:    $REP (int: $REP_INT)"
echo "System XML:   $SYSTEM_XML"
echo "Topology PDB: $TOPOLOGY_PDB"
echo "Minimized PDB: $MINIMIZED_PDB"
echo "Output JSON:  $OUTPUT_JSON (test output, will be deleted)"
echo ""

# Verify files exist
MISSING=0
for f in "$SYSTEM_XML" "$TOPOLOGY_PDB" "$MINIMIZED_PDB"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing file: $f"
        MISSING=1
    else
        echo "✓ Found: $f"
    fi
done

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "Cannot proceed — missing input files."
    exit 1
fi

echo ""
echo "Running short test MD (0.01 ns, ~5000 steps)..."
echo ""

python3 -m nnrti.md.sherlock.run_md_job \
    --mutation "$MUTATION" \
    --replicate $REP_INT \
    --task-id 0 \
    --system-xml "$SYSTEM_XML" \
    --topology-pdb "$TOPOLOGY_PDB" \
    --minimized-pdb "$MINIMIZED_PDB" \
    --output-json "$OUTPUT_JSON" \
    --production-ns 0.01 \
    --no-resume \
    --force

STATUS=$?

echo ""
if [ $STATUS -eq 0 ]; then
    echo "=========================================="
    echo "✓ Test PASSED"
    echo "=========================================="
    echo ""
    echo "Output JSON:"
    cat "$OUTPUT_JSON"
    echo ""
    echo ""
    echo "Cleaning up test output..."
    rm -f "$OUTPUT_JSON"
    # Clean up the short test output files
    rm -f "${PARENT}/${MUTATION}_rep${REP_INT:+$(printf '%02d' $REP_INT)}_md_final.pdb"
    rm -f "${PARENT}/${MUTATION}_rep${REP_INT:+$(printf '%02d' $REP_INT)}_md_state.csv"
    rm -f "${PARENT}/${MUTATION}_rep${REP_INT:+$(printf '%02d' $REP_INT)}_md.chk"
    rm -f "${PARENT}/${MUTATION}_rep${REP_INT:+$(printf '%02d' $REP_INT)}_analysis.dcd"
    rm -f "${PARENT}/${MUTATION}_rep${REP_INT:+$(printf '%02d' $REP_INT)}_analysis_topology.pdb"
    echo "Ready to submit batch jobs with: ./scripts/sherlock/submit_md_batched.sh 6 12"
else
    echo "=========================================="
    echo "✗ Test FAILED (exit code: $STATUS)"
    echo "=========================================="
    echo "Check the output above for errors."
    rm -f "$OUTPUT_JSON"
fi

exit $STATUS
