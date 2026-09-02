#!/bin/bash
#
# Submit the full P0 NEQ pipeline in one command (SLURM dependency chain).
#
# Stages run automatically in order: em → equil → extract → switch
# EM/extract use normal partition (CPU); equil/switch use gpu.
#
# Usage:
#   NEQ_SNAPSHOTS=100 REPLICATES=3 FORCE=1 bash ops/slurm/fep/prepare_p0_neq.sh
#   bash ops/slurm/fep/submit_p0_neq_pipeline.sh
#
# Optional:
#   SHERLOCK_MAX_CONCURRENT=20
#   scancel the EM job and re-run pipeline if manifest changes mid-flight
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

_submit() {
    local stage="$1"
    local dep="${2:-}"
    if [[ -n "${dep}" ]]; then
        DEPENDENCY="${dep}" STAGE="${stage}" bash ops/slurm/fep/submit_p0_neq.sh | tail -1
    else
        STAGE="${stage}" bash ops/slurm/fep/submit_p0_neq.sh | tail -1
    fi
}

echo "=========================================="
echo "P0 NEQ pipeline submit (em → equil → extract → switch)"
echo "=========================================="

EM_JOB="$(_submit em)"
echo "  em:      ${EM_JOB}  (normal/CPU)"

EQUIL_JOB="$(_submit equil "afterok:${EM_JOB}")"
echo "  equil:   ${EQUIL_JOB}  (gpu, after em)"

EXTRACT_JOB="$(_submit extract "afterok:${EQUIL_JOB}")"
echo "  extract: ${EXTRACT_JOB}  (normal/CPU, after equil)"

SWITCH_JOBS="$(_submit switch "afterok:${EXTRACT_JOB}")"
echo "  switch:  ${SWITCH_JOBS}  (gpu, after extract; may be multiple array jobs)"

echo ""
echo "Pipeline queued. Monitor:"
echo "  squeue -u \$USER"
echo "  sacct -j ${EM_JOB},${EQUIL_JOB},${EXTRACT_JOB} --format=JobID,JobName,State,Elapsed"
echo "  # switch job ids: ${SWITCH_JOBS}"
echo ""
echo "Cancel entire pipeline:"
echo "  scancel ${EM_JOB} ${EQUIL_JOB} ${EXTRACT_JOB} ${SWITCH_JOBS}"
