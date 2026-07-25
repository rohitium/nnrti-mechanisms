#!/bin/bash
#
# Prepare NEQ inputs + panel manifest for P0 legs.
#
# Usage:
#   bash scripts/fep_pmx/prepare_p0_neq.sh
#   NEQ_SNAPSHOTS=10 REPLICATES=1 FORCE=1 bash scripts/fep_pmx/prepare_p0_neq.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

NEQ_SNAPSHOTS="${NEQ_SNAPSHOTS:-100}"
REPLICATES="${REPLICATES:-1}"
FORCE="${FORCE:-0}"

args=(--replicates "${REPLICATES}" --n-snapshots "${NEQ_SNAPSHOTS}")
if [[ "${FORCE}" == "1" ]]; then
    args+=(--force)
fi

python scripts/fep_pmx/prepare_neq.py "${args[@]}"
