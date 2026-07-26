#!/bin/bash
#
# Submit pmx NEQ jobs for P0 legs on Sherlock.
#
# Single stage:
#   STAGE=em      bash scripts/fep_pmx/submit_p0_neq.sh
#   STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh
#
# Full pipeline (one command, SLURM dependency chain):
#   bash scripts/fep_pmx/submit_p0_neq_pipeline.sh
#
# Smoke test:
#   NEQ_SNAPSHOTS=10 REPLICATES=1 bash scripts/fep_pmx/prepare_p0_neq.sh
#   STAGE=switch bash scripts/fep_pmx/submit_p0_neq.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

STAGE="${STAGE:-${1:-switch}}"
MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_panel_manifest.csv}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-20}"
SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-4}"
DEPENDENCY="${DEPENDENCY:-}"

case "${STAGE}" in
    em)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-normal}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-02:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
        ;;
    equil)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-08:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
        ;;
    extract)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-normal}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-02:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
        ;;
    switch)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-02:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
        ;;
    *)
        echo "Unknown STAGE=${STAGE}. Use em|equil|extract|switch." >&2
        exit 1
        ;;
esac

if [[ ! -f "${MANIFEST}" ]]; then
    echo "Missing manifest: ${MANIFEST}. Run prepare_p0_neq.sh first." >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch required; run on Sherlock." >&2
    exit 1
fi

if command -v module >/dev/null 2>&1; then
    module load python/3.9.0 2>/dev/null || module load python/3.12.1 2>/dev/null || true
fi
PYTHON="${PYTHON:-python3}"

TASK_ID_FILE="${PROJECT_ROOT}/results/analysis/fep_pmx/neq_${STAGE}_task_ids.txt"
TASK_COUNT="$(
"${PYTHON}" - <<PY
import csv
from pathlib import Path

manifest = Path("${MANIFEST}")
stage = "${STAGE}"
out = Path("${TASK_ID_FILE}")
ids = []
with manifest.open(newline="") as handle:
    reader = csv.DictReader(handle)
    key = "panel_task_id" if "panel_task_id" in reader.fieldnames else "task_id"
    for row in reader:
        if row["stage"] == stage:
            ids.append(int(row[key]))
if not ids:
    raise SystemExit(f"No tasks for stage={stage} in {manifest}")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(str(i) for i in ids) + "\n")
print(len(ids))
PY
)"
ARRAY_MAX=$((TASK_COUNT - 1))
ARRAY_SPEC="0-${ARRAY_MAX}%${SHERLOCK_MAX_CONCURRENT}"

PMX_VENV="${PMX_VENV:-$HOME/.venvs/pmx}"
if [[ -z "${GMXLIB:-}" ]] && [[ -f "${PMX_VENV}/bin/activate" ]]; then
    # shellcheck disable=SC1090
    source "${PMX_VENV}/bin/activate"
fi
if [[ -z "${GMXLIB:-}" ]]; then
    GMXLIB="$("${PYTHON}" - <<'PY'
import os
import pmx
print(os.path.join(os.path.dirname(pmx.__file__), "data", "mutff"))
PY
)"
fi
if [[ -z "${GMXLIB}" || ! -d "${GMXLIB}/amber14sbmut.ff" ]]; then
    echo "ERROR: GMXLIB not set to pmx mutff (got: ${GMXLIB:-<empty>})" >&2
    echo "  module load python/3.9.0 && source ~/.venvs/pmx/bin/activate" >&2
    exit 1
fi

mkdir -p logs
echo "Manifest:  ${MANIFEST}"
echo "Stage:     ${STAGE}  (${TASK_COUNT} tasks)"
echo "Task ids:  ${TASK_ID_FILE}"
echo "Array:     ${ARRAY_SPEC}"
echo "GMXLIB:    ${GMXLIB}"
echo "Partition: ${SHERLOCK_PARTITION}  GRES: ${SHERLOCK_GRES:-<none>}  TIME: ${SHERLOCK_TIME}"
if [[ -n "${DEPENDENCY}" ]]; then
    echo "Depends:   ${DEPENDENCY}"
fi

SBATCH_ARGS=(
    --parsable
    --job-name="pmx_neq_${STAGE}"
    --partition="${SHERLOCK_PARTITION}"
    --time="${SHERLOCK_TIME}"
    --mem="${SHERLOCK_MEM}"
    --cpus-per-task="${SHERLOCK_CPUS_PER_TASK}"
    --array="${ARRAY_SPEC}"
    --output="${PROJECT_ROOT}/logs/pmx_neq_${STAGE}.%A_%a.out"
    --error="${PROJECT_ROOT}/logs/pmx_neq_${STAGE}.%A_%a.err"
)
if [[ -n "${SHERLOCK_GRES}" ]]; then
    SBATCH_ARGS+=(--gres="${SHERLOCK_GRES}")
fi
if [[ -n "${SHERLOCK_QOS}" ]]; then
    SBATCH_ARGS+=(--qos="${SHERLOCK_QOS}")
fi
if [[ -n "${DEPENDENCY}" ]]; then
    SBATCH_ARGS+=(--dependency="${DEPENDENCY}")
fi

JOB_ID="$(
sbatch "${SBATCH_ARGS[@]}" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

source ${PROJECT_ROOT}/scripts/sherlock/load_gromacs_module.sh
export GMXLIB=${GMXLIB}

cd ${PROJECT_ROOT}

module load python/3.9.0 2>/dev/null || true
if [[ -f "${HOME}/.venvs/pmx/bin/activate" ]]; then
  source "${HOME}/.venvs/pmx/bin/activate"
fi

TASK_ID=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" ${TASK_ID_FILE})
if [[ -z "\${TASK_ID}" ]]; then
  echo "ERROR: no task id for array index \${SLURM_ARRAY_TASK_ID} in ${TASK_ID_FILE}" >&2
  exit 1
fi

python3 scripts/fep_pmx/run_neq_task.py \
    --manifest ${MANIFEST} \
    --task-id \${TASK_ID}
SBATCH_EOF
)"

echo "Submitted batch job ${JOB_ID}"
echo "${JOB_ID}"
