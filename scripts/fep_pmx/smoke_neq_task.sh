#!/bin/bash
#
# Run one NEQ manifest task inside an interactive GPU allocation (any stage).
#
# Usage (after salloc_neq_gpu.sh):
#   bash scripts/fep_pmx/smoke_neq_task.sh 1          # equil lambda0, V106A holo rep1
#   STAGE=em LEG=wt_to_V106A PHASE=apo bash scripts/fep_pmx/smoke_neq_task.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_panel_manifest.csv}"
LEG="${LEG:-wt_to_V106A}"
PHASE="${PHASE:-holo}"
REP="${REP:-1}"
STAGE="${STAGE:-}"
TASK_ID="${1:-${TASK_ID:-}}"

source scripts/sherlock/load_gromacs_module.sh
source scripts/sherlock/activate_pmx_env.sh

PYTHON="${PYTHON:-python3}"

if [[ -z "${TASK_ID}" ]]; then
    STAGE="${STAGE:-equil}"
    TASK_ID="$("${PYTHON}" - <<PY
import csv
from pathlib import Path
manifest = Path("${MANIFEST}")
leg, phase, rep, stage = "${LEG}", "${PHASE}", int("${REP}"), "${STAGE}"
with manifest.open(newline="") as handle:
    reader = csv.DictReader(handle)
    key = "panel_task_id" if "panel_task_id" in reader.fieldnames else "task_id"
    for row in reader:
        if (
            row["leg_id"] == leg
            and row["phase"] == phase
            and int(row["replicate"]) == rep
            and row["stage"] == stage
        ):
            print(row[key])
            break
    else:
        raise SystemExit(f"No {stage} task for {leg} {phase} rep{rep}")
PY
)"
fi

echo "Running NEQ task_id=${TASK_ID} (manifest=${MANIFEST})"
"${PYTHON}" src/nnrti/fep/run_neq_task.py --manifest "${MANIFEST}" --task-id "${TASK_ID}"
