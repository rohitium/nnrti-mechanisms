#!/bin/bash
#
# Run one short FEP lambda window on the current node (interactive GPU or CPU).
# Use this inside salloc before submitting the full V106A array.
#
# Usage (on Sherlock GPU node, after loading OpenMM module):
#   source scripts/sherlock/load_openmm_module.sh
#   bash scripts/sherlock/run_fep_jorgensen_pilot.sh
#
# Optional env vars:
#   FEP_LEG_ID=wt_to_V106A
#   FEP_STATE_INDEX=0
#   FEP_PLATFORM=CUDA
#   FEP_EQUIL_STEPS=5000
#   FEP_PROD_STEPS=25000
#   FEP_ENERGY_INTERVAL=2500
#   FEP_CHECKPOINT_INTERVAL=25000
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

# shellcheck source=load_openmm_module.sh
source "$SCRIPT_DIR/load_openmm_module.sh"

FEP_LEG_ID="${FEP_LEG_ID:-wt_to_V106A}"
FEP_STATE_INDEX="${FEP_STATE_INDEX:-0}"
FEP_PLATFORM="${FEP_PLATFORM:-CUDA}"
FEP_EQUIL_STEPS="${FEP_EQUIL_STEPS:-5000}"
FEP_PROD_STEPS="${FEP_PROD_STEPS:-25000}"
FEP_ENERGY_INTERVAL="${FEP_ENERGY_INTERVAL:-2500}"
FEP_CHECKPOINT_INTERVAL="${FEP_CHECKPOINT_INTERVAL:-25000}"

PHASE_DIR="$PROJECT_ROOT/results/analysis/fep_jorgensen/legs/${FEP_LEG_ID}/holo"
WINDOW_DIR="$PHASE_DIR/windows"

for required in "$PHASE_DIR/hybrid_system.xml" "$PHASE_DIR/hybrid_topology.pdb" "$PHASE_DIR/schedule.json"; do
    if [[ ! -f "$required" ]]; then
        echo "Missing prepared holo artifact: $required" >&2
        echo "Run Perses prepare locally and rsync legs/${FEP_LEG_ID}/ to Sherlock." >&2
        exit 1
    fi
done

echo "Project:      $PROJECT_ROOT"
echo "Leg:          $FEP_LEG_ID"
echo "State index:  $FEP_STATE_INDEX"
echo "Platform:     $FEP_PLATFORM"
echo "Equil steps:  $FEP_EQUIL_STEPS"
echo "Prod steps:   $FEP_PROD_STEPS"
echo

python3 -m scripts.fep_jorgensen.worker \
    --phase-dir "$PHASE_DIR" \
    --output-dir "$WINDOW_DIR" \
    --state-index "$FEP_STATE_INDEX" \
    --platform "$FEP_PLATFORM" \
    --equilibration-steps "$FEP_EQUIL_STEPS" \
    --production-steps "$FEP_PROD_STEPS" \
    --energy-interval "$FEP_ENERGY_INTERVAL" \
    --checkpoint-interval "$FEP_CHECKPOINT_INTERVAL"

CSV="$WINDOW_DIR/state_$(printf '%02d' "$FEP_STATE_INDEX")_energies.csv"
echo
echo "Wrote: $CSV"
echo "Rows:  $(tail -n +2 "$CSV" | wc -l | tr -d ' ')"
