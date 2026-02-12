#!/bin/bash
# Robust analysis pipeline with checkpointing
# Usage: ./scripts/run_analysis.sh [--force]

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

FORCE_FLAG=""
if [[ "$1" == "--force" ]]; then
    FORCE_FLAG="--force"
    echo "Force mode: will recompute all steps"
fi

echo "============================================"
echo "Robust Analysis Pipeline with Checkpointing"
echo "============================================"
echo

# Step 1: Collect MD metadata
echo "STEP 1: Collecting MD metadata..."
python -m src.analysis.cli.analyze_incremental --step collect $FORCE_FLAG
echo

# Step 2: Compute MM/GBSA (incremental, safe)
echo "STEP 2: Computing MM/GBSA metrics (incremental)..."
python -m src.analysis.cli.compute_mmgbsa_safe $FORCE_FLAG
echo

# Step 3: Compute structural metrics
echo "STEP 3: Computing structural metrics..."
python -m src.analysis.cli.analyze_incremental --step metrics $FORCE_FLAG
echo

# Step 4: Generate plots
echo "STEP 4: Generating plots..."
python -m src.analysis.cli.analyze_incremental --step plots $FORCE_FLAG
echo

echo "============================================"
echo "Analysis complete!"
echo "Results: results/"
echo "Plots: results/plots/"
echo "============================================"
