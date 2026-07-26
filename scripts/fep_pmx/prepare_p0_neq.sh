#!/bin/bash
#
# Prepare NEQ inputs + panel manifest for P0 legs.
#
# Usage:
#   bash scripts/fep_pmx/prepare_p0_neq.sh
#   NEQ_SNAPSHOTS=5 REPLICATES=3 FORCE=1 bash scripts/fep_pmx/prepare_p0_neq.sh  # 120 switch tasks
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

NEQ_SNAPSHOTS="${NEQ_SNAPSHOTS:-5}"
REPLICATES="${REPLICATES:-1}"
FORCE="${FORCE:-0}"

if command -v module >/dev/null 2>&1; then
    module load python/3.9.0 2>/dev/null || module load python/3.12.1 2>/dev/null || true
fi
PYTHON="${PYTHON:-python3}"

args=(--replicates "${REPLICATES}" --n-snapshots "${NEQ_SNAPSHOTS}")
if [[ "${FORCE}" == "1" ]]; then
    args+=(--force)
fi

"${PYTHON}" scripts/fep_pmx/prepare_neq.py "${args[@]}"
