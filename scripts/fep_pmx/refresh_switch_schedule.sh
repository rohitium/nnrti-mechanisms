#!/bin/bash
#
# Shrink switch schedule mid-pipeline without canceling em/equil/extract jobs.
#
# Safe while 35807105/106/112 are running:
#   1) Updates per-leg neq_prepare.json → 5 snapshot times (extract pulls 5 not 100)
#   2) Does NOT overwrite neq_panel_manifest.csv (pending array jobs keep valid task ids)
#
# After extract finishes:
#   NEQ_SNAPSHOTS=5 REPLICATES=3 bash scripts/fep_pmx/prepare_p0_neq.sh --rebuild-panel-only
#   DEPENDENCY=afterok:<extract_job> STAGE=switch bash scripts/fep_pmx/submit_p0_neq.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

NEQ_SNAPSHOTS="${NEQ_SNAPSHOTS:-5}"
REPLICATES="${REPLICATES:-3}"

if command -v module >/dev/null 2>&1; then
    module load python/3.9.0 2>/dev/null || true
fi
PYTHON="${PYTHON:-python3}"

"${PYTHON}" scripts/fep_pmx/prepare_neq.py \
    --refresh-switch-only \
    --replicates "${REPLICATES}" \
    --n-snapshots "${NEQ_SNAPSHOTS}"

echo ""
echo "Per-leg switch schedule → ${NEQ_SNAPSHOTS} snapshots (120 switch tasks for P0)."
echo "Panel manifest NOT updated — pending em/equil/extract array ids still valid."
echo ""
echo "After extract stage completes:"
echo "  NEQ_SNAPSHOTS=${NEQ_SNAPSHOTS} REPLICATES=${REPLICATES} bash scripts/fep_pmx/prepare_p0_neq.sh --rebuild-panel-only"
echo "  DEPENDENCY=afterok:<extract_job_id> STAGE=switch bash scripts/fep_pmx/submit_p0_neq.sh"
