#!/bin/bash
#
# Submit pmx NEQ jobs for P0 legs on Sherlock.
#
# GPU QOS limits (partition gpu → qos gpu):
#   MaxSubmitPU=100, MaxTRESPU cpu=128 + gres/gpu=16, MaxWall=48h
# Array elements count individually against MaxSubmitPU — bundle switches per task.
#
# Single stage:
#   STAGE=em      bash scripts/fep_pmx/submit_p0_neq.sh
#   STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh
#
# Full pipeline (one command, SLURM dependency chain):
#   bash scripts/fep_pmx/submit_p0_neq_pipeline.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

STAGE="${STAGE:-${1:-switch}}"
MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_panel_manifest.csv}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
DEPENDENCY="${DEPENDENCY:-}"

case "${STAGE}" in
    em)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-normal}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-04:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
        SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-4}"
        SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-64}"
        SHERLOCK_MAX_ARRAY_SIZE="${SHERLOCK_MAX_ARRAY_SIZE:-512}"
        SHERLOCK_CHAIN_CHUNKS="${SHERLOCK_CHAIN_CHUNKS:-0}"
        ;;
    equil)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-12:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
        SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-8}"
        SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-16}"
        SHERLOCK_MAX_ARRAY_SIZE="${SHERLOCK_MAX_ARRAY_SIZE:-90}"
        SHERLOCK_CHAIN_CHUNKS="${SHERLOCK_CHAIN_CHUNKS:-0}"
        ;;
    extract)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-normal}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-04:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
        SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-4}"
        SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-64}"
        SHERLOCK_MAX_ARRAY_SIZE="${SHERLOCK_MAX_ARRAY_SIZE:-512}"
        SHERLOCK_CHAIN_CHUNKS="${SHERLOCK_CHAIN_CHUNKS:-0}"
        ;;
    switch)
        SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
        SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
        SHERLOCK_TIME="${SHERLOCK_TIME:-48:00:00}"
        SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
        SHERLOCK_CPUS_PER_TASK="${SHERLOCK_CPUS_PER_TASK:-8}"
        SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-16}"
        SHERLOCK_MAX_ARRAY_SIZE="${SHERLOCK_MAX_ARRAY_SIZE:-90}"
        SHERLOCK_CHAIN_CHUNKS="${SHERLOCK_CHAIN_CHUNKS:-1}"
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
CHUNK_MANIFEST="$(
"${PYTHON}" - <<PY
import csv
import json
from pathlib import Path

manifest = Path("${MANIFEST}")
stage = "${STAGE}"
out = Path("${TASK_ID_FILE}")
chunk_size = int("${SHERLOCK_MAX_ARRAY_SIZE}")
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

chunks = []
for start in range(0, len(ids), chunk_size):
    chunk = ids[start : start + chunk_size]
    chunk_path = out.parent / f"{out.stem}_chunk{len(chunks):03d}{out.suffix}"
    chunk_path.write_text("\n".join(str(i) for i in chunk) + "\n")
    chunks.append(
        {
            "file": str(chunk_path),
            "count": len(chunk),
            "array": f"0-{len(chunk) - 1}",
        }
    )
print(json.dumps({"total": len(ids), "chunks": chunks}))
PY
)"

TASK_COUNT="$(echo "${CHUNK_MANIFEST}" | "${PYTHON}" -c "import json,sys; print(json.load(sys.stdin)['total'])")"
CHUNK_COUNT="$(echo "${CHUNK_MANIFEST}" | "${PYTHON}" -c "import json,sys; print(len(json.load(sys.stdin)['chunks']))")"

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

LOG_DIR="${PROJECT_ROOT}/logs/pmx_neq/${STAGE}"
mkdir -p "${LOG_DIR}"
echo "Logs:      ${LOG_DIR}/"
echo "Manifest:  ${MANIFEST}"
echo "Stage:     ${STAGE}  (${TASK_COUNT} tasks, ${CHUNK_COUNT} array job(s))"
echo "Task ids:  ${TASK_ID_FILE}"
echo "GMXLIB:    ${GMXLIB}"
echo "Partition: ${SHERLOCK_PARTITION}  GRES: ${SHERLOCK_GRES:-<none>}  TIME: ${SHERLOCK_TIME}"
echo "Resources: cpus/task=${SHERLOCK_CPUS_PER_TASK}  mem=${SHERLOCK_MEM}  array%=${SHERLOCK_MAX_CONCURRENT}"
if [[ -n "${DEPENDENCY}" ]]; then
    echo "Depends:   ${DEPENDENCY}"
fi
if [[ "${SHERLOCK_CHAIN_CHUNKS}" == "1" && "${CHUNK_COUNT}" -gt 1 ]]; then
    echo "Chaining:  ${CHUNK_COUNT} switch arrays sequentially (afterany)"
fi

_submit_chunk() {
    local chunk_json="$1"
    local chunk_dep="${2:-}"

    local chunk_file chunk_array chunk_n chunk_idx
    chunk_file="$(echo "${chunk_json}" | "${PYTHON}" -c "import json,sys; print(json.load(sys.stdin)['file'])")"
    chunk_array="$(echo "${chunk_json}" | "${PYTHON}" -c "import json,sys; print(json.load(sys.stdin)['array'])")"
    chunk_n="$(echo "${chunk_json}" | "${PYTHON}" -c "import json,sys; print(json.load(sys.stdin)['count'])")"
    chunk_idx="${3:-0}"
    local array_spec="${chunk_array}%${SHERLOCK_MAX_CONCURRENT}"

    echo "Chunk ${chunk_idx}: ${chunk_n} tasks  array=${array_spec}  file=$(basename "${chunk_file}")"
    if [[ -n "${chunk_dep}" ]]; then
        echo "  depends: ${chunk_dep}"
    fi

    local SBATCH_ARGS=(
        --parsable
        --job-name="pmx_neq_${STAGE}"
        --partition="${SHERLOCK_PARTITION}"
        --time="${SHERLOCK_TIME}"
        --mem="${SHERLOCK_MEM}"
        --cpus-per-task="${SHERLOCK_CPUS_PER_TASK}"
        --array="${array_spec}"
        --requeue
        --open-mode=append
        --output="${LOG_DIR}/pmx_neq_${STAGE}_c${chunk_idx}.%A_%a.out"
        --error="${LOG_DIR}/pmx_neq_${STAGE}_c${chunk_idx}.%A_%a.err"
    )
    if [[ -n "${SHERLOCK_GRES}" ]]; then
        SBATCH_ARGS+=(--gres="${SHERLOCK_GRES}")
    fi
    if [[ -n "${SHERLOCK_QOS}" ]]; then
        SBATCH_ARGS+=(--qos="${SHERLOCK_QOS}")
    fi
    if [[ -n "${chunk_dep}" ]]; then
        SBATCH_ARGS+=(--dependency="${chunk_dep}")
    fi

    sbatch "${SBATCH_ARGS[@]}" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

source ${PROJECT_ROOT}/scripts/sherlock/load_gromacs_module.sh
export GMXLIB=${GMXLIB}
export GMX_NTOMP=${SHERLOCK_CPUS_PER_TASK}

cd ${PROJECT_ROOT}

module load python/3.9.0 2>/dev/null || true
if [[ -f "${HOME}/.venvs/pmx/bin/activate" ]]; then
  source "${HOME}/.venvs/pmx/bin/activate"
fi

TASK_ID=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" ${chunk_file})
if [[ -z "\${TASK_ID}" ]]; then
  echo "ERROR: no task id for array index \${SLURM_ARRAY_TASK_ID} in ${chunk_file}" >&2
  exit 1
fi

python3 scripts/fep_pmx/run_neq_task.py \
    --manifest ${MANIFEST} \
    --task-id \${TASK_ID}
SBATCH_EOF
}

JOB_IDS=()
PREV_DEP="${DEPENDENCY}"
chunk_idx=0
while IFS= read -r chunk_json; do
    job_id="$(_submit_chunk "${chunk_json}" "${PREV_DEP}" "${chunk_idx}")"
    echo "  submitted ${job_id}"
    JOB_IDS+=("${job_id}")
    if [[ "${SHERLOCK_CHAIN_CHUNKS}" == "1" ]]; then
        PREV_DEP="afterany:${job_id}"
    fi
    chunk_idx=$((chunk_idx + 1))
done < <(echo "${CHUNK_MANIFEST}" | "${PYTHON}" -c "import json,sys; print('\n'.join(json.dumps(c) for c in json.load(sys.stdin)['chunks']))")

echo "Submitted ${#JOB_IDS[@]} batch job(s): ${JOB_IDS[*]}"
echo "${JOB_IDS[*]}"
