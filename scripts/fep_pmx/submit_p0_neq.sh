#!/bin/bash
#
# Submit pmx NEQ jobs for P0 legs on Sherlock (GPU).
#
# Run in order:
#   bash scripts/fep_pmx/prepare_p0_neq.sh
#   STAGE=em      bash scripts/fep_pmx/submit_p0_neq.sh
#   STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh
#   STAGE=extract bash scripts/fep_pmx/submit_p0_neq.sh
#   STAGE=switch  bash scripts/fep_pmx/submit_p0_neq.sh
#
# Smoke test (few snapshots):
#   NEQ_SNAPSHOTS=10 REPLICATES=1 bash scripts/fep_pmx/prepare_p0_neq.sh
#   STAGE=switch NEQ_SNAPSHOTS=10 bash scripts/fep_pmx/submit_p0_neq.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

STAGE="${STAGE:-${1:-switch}}"
MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_panel_manifest.csv}"
SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
SHERLOCK_MAX_CONCURRENT="${SHERLOCK_MAX_CONCURRENT:-20}"

case "${STAGE}" in
    em)      SHERLOCK_TIME="${SHERLOCK_TIME:-01:00:00}"; SHERLOCK_MEM="${SHERLOCK_MEM:-16G}" ;;
    equil)   SHERLOCK_TIME="${SHERLOCK_TIME:-06:00:00}"; SHERLOCK_MEM="${SHERLOCK_MEM:-32G}" ;;
    extract) SHERLOCK_TIME="${SHERLOCK_TIME:-01:00:00}"; SHERLOCK_MEM="${SHERLOCK_MEM:-16G}" ;;
    switch)  SHERLOCK_TIME="${SHERLOCK_TIME:-01:00:00}"; SHERLOCK_MEM="${SHERLOCK_MEM:-16G}" ;;
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

TASK_IDS="$(
"${PYTHON}" - <<PY
import csv
from pathlib import Path

manifest = Path("${MANIFEST}")
stage = "${STAGE}"
ids = []
with manifest.open(newline="") as handle:
    reader = csv.DictReader(handle)
    key = "panel_task_id" if "panel_task_id" in reader.fieldnames else "task_id"
    for row in reader:
        if row["stage"] == stage:
            ids.append(int(row[key]))
if not ids:
    raise SystemExit(f"No tasks for stage={stage} in {manifest}")
print(",".join(str(i) for i in ids))
PY
)"

mkdir -p logs
echo "Manifest:  ${MANIFEST}"
echo "Stage:     ${STAGE}"
echo "Tasks:     ${TASK_IDS}"
echo "Partition: ${SHERLOCK_PARTITION}  GRES: ${SHERLOCK_GRES}  TIME: ${SHERLOCK_TIME}"

sbatch \
    --job-name="pmx_neq_${STAGE}" \
    --partition="${SHERLOCK_PARTITION}" \
    --gres="${SHERLOCK_GRES}" \
    ${SHERLOCK_QOS:+--qos="${SHERLOCK_QOS}"} \
    --time="${SHERLOCK_TIME}" \
    --mem="${SHERLOCK_MEM}" \
    --array="${TASK_IDS}%${SHERLOCK_MAX_CONCURRENT}" \
    --output="${PROJECT_ROOT}/logs/pmx_neq_${STAGE}.%A_%a.out" \
    --error="${PROJECT_ROOT}/logs/pmx_neq_${STAGE}.%A_%a.err" \
    <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

source ${PROJECT_ROOT}/scripts/sherlock/load_gromacs_module.sh

cd ${PROJECT_ROOT}

python3 scripts/fep_pmx/run_neq_task.py \
    --manifest ${MANIFEST} \
    --task-id \${SLURM_ARRAY_TASK_ID}
SBATCH_EOF
