#!/bin/bash
#
# Interactive smoke test: one NEQ EM task on a GPU node (validate before sbatch).
#
# Usage (inside salloc from salloc_neq_gpu.sh):
#   bash scripts/fep_pmx/smoke_neq_em.sh
#   LEG=wt_to_V106A PHASE=holo REP=1 TASK_ID=0 bash scripts/fep_pmx/smoke_neq_em.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

LEG="${LEG:-wt_to_V106A}"
PHASE="${PHASE:-holo}"
REP="${REP:-1}"
TASK_ID="${TASK_ID:-}"
MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_panel_manifest.csv}"

# GROMACS first (module purge); pmx python module must load after that.
source scripts/sherlock/load_gromacs_module.sh
source scripts/sherlock/activate_pmx_env.sh

PYTHON="${PYTHON:-python3}"

# Sync topology + mdps (no --force unless FORCE=1; never wipes completed outputs).
PREP_ARGS=(--leg "${LEG}" --phase "${PHASE}" --replicate "${REP}" --n-snapshots "${NEQ_SNAPSHOTS:-5}")
if [[ "${FORCE:-0}" == "1" ]]; then
    PREP_ARGS+=(--force)
fi
"${PYTHON}" scripts/fep_pmx/prepare_neq.py "${PREP_ARGS[@]}"

if [[ -z "${TASK_ID}" ]]; then
    TASK_ID="$("${PYTHON}" - <<PY
import csv
from pathlib import Path
manifest = Path("${MANIFEST}")
leg, phase, rep, stage = "${LEG}", "${PHASE}", int("${REP}"), "em"
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
        raise SystemExit("EM task not found in manifest")
PY
)"
fi

echo "Running NEQ smoke: leg=${LEG} phase=${PHASE} rep=${REP} task_id=${TASK_ID}"
"${PYTHON}" scripts/fep_pmx/run_neq_task.py --manifest "${MANIFEST}" --task-id "${TASK_ID}"

OUT="results/analysis/fep_pmx/legs/${LEG}/${PHASE}/rep_$(printf '%02d' "${REP}")/neq/em/em.gro"
if [[ -f "${OUT}" ]]; then
    echo "OK: ${OUT}"
else
    echo "FAIL: missing ${OUT}" >&2
    exit 1
fi
