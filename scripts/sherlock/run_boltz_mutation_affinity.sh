#!/bin/bash
#
# Run a single Boltz-2 affinity prediction for a mutation (or WT) on Sherlock.
#
# Usage:
#   bash scripts/sherlock/run_boltz_mutation_affinity.sh K103N
#   bash scripts/sherlock/run_boltz_mutation_affinity.sh WT
#
# Optional env vars:
#   BOLTZ_MUTATION=K103N
#   BOLTZ_PYTHON=python3
#   BOLTZ_CHAINS=A,B
#   BOLTZ_INPUT_YAML=inputs/boltz/K103N_rt_dor_affinity.yaml
#   BOLTZ_OUT_DIR=/scratch/users/$USER/nnrti-mechanisms/results/boltz/K103N_affinity
#   BOLTZ_CACHE_DIR=/scratch/users/$USER/.boltz
#   BOLTZ_EXTRA_ARGS="--seed 1001"
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "${PROJECT_ROOT}"

normalize_label() {
    local raw="$1"
    raw="${raw//+/_}"
    raw="${raw//\//_}"
    raw="${raw// /}"
    raw="${raw//-/_}"
    printf '%s' "${raw}" | tr -cd '[:alnum:]_'
}

MUTATION="${1:-${BOLTZ_MUTATION:-}}"
if [ -z "${MUTATION}" ]; then
    echo "Usage: bash scripts/sherlock/run_boltz_mutation_affinity.sh <MUTATION|WT>" >&2
    exit 2
fi

BOLTZ_PYTHON="${BOLTZ_PYTHON:-python3}"
BOLTZ_CHAINS="${BOLTZ_CHAINS:-A,B}"

MUT_UPPER="$(printf '%s' "${MUTATION}" | tr '[:lower:]' '[:upper:]')"
if [ "${MUT_UPPER}" = "WT" ]; then
    SAFE_LABEL="wt"
    DEFAULT_INPUT_CIF="data/prepared/dor_4ncg/wt_4ncg.cif"
else
    SAFE_LABEL="$(normalize_label "${MUTATION}")"
    SAFE_LABEL="$(printf '%s' "${SAFE_LABEL}" | tr '[:lower:]' '[:upper:]')"
    if [ -z "${SAFE_LABEL}" ]; then
        echo "ERROR: mutation label '${MUTATION}' is invalid after normalization." >&2
        exit 2
    fi
    DEFAULT_INPUT_CIF="data/prepared/dor_4ncg/mut_${SAFE_LABEL}.cif"
fi

if [ ! -f "${DEFAULT_INPUT_CIF}" ]; then
    echo "ERROR: input CIF not found: ${DEFAULT_INPUT_CIF}" >&2
    echo "Available prepared CIF files:" >&2
    find data/prepared/dor_4ncg -maxdepth 1 -type f -name "*.cif" | sort >&2
    exit 1
fi

TARGET_ID="${SAFE_LABEL}_rt_dor_affinity"
if [ "${SAFE_LABEL}" = "wt" ]; then
    TARGET_ID="wt_rt_dor_affinity"
fi

if [ -d "/scratch/users/${USER}" ]; then
    DEFAULT_BASE="/scratch/users/${USER}/nnrti-mechanisms"
    DEFAULT_CACHE="/scratch/users/${USER}/.boltz"
else
    DEFAULT_BASE="${PROJECT_ROOT}"
    DEFAULT_CACHE="${PROJECT_ROOT}/.cache/boltz"
fi

BOLTZ_INPUT_CIF="${BOLTZ_INPUT_CIF:-${DEFAULT_INPUT_CIF}}"
BOLTZ_INPUT_YAML="${BOLTZ_INPUT_YAML:-inputs/boltz/${TARGET_ID}.yaml}"
BOLTZ_OUT_DIR="${BOLTZ_OUT_DIR:-${DEFAULT_BASE}/results/boltz/${SAFE_LABEL}_affinity}"
BOLTZ_CACHE_DIR="${BOLTZ_CACHE_DIR:-${DEFAULT_CACHE}}"

mkdir -p "$(dirname "${BOLTZ_INPUT_YAML}")"

echo "Preparing Boltz input for mutation: ${MUTATION}"
echo "Input CIF:  ${BOLTZ_INPUT_CIF}"
echo "Input YAML: ${BOLTZ_INPUT_YAML}"

"${BOLTZ_PYTHON}" scripts/sherlock/make_boltz_wt_input.py \
    --input-cif "${BOLTZ_INPUT_CIF}" \
    --output-yaml "${BOLTZ_INPUT_YAML}" \
    --chains "${BOLTZ_CHAINS}"

SHERLOCK_MODULES="${SHERLOCK_MODULES:-}" \
BOLTZ_PYTHON="${BOLTZ_PYTHON}" \
BOLTZ_BIN="${BOLTZ_BIN:-boltz}" \
BOLTZ_INPUT_YAML="${BOLTZ_INPUT_YAML}" \
BOLTZ_OUT_DIR="${BOLTZ_OUT_DIR}" \
BOLTZ_CACHE_DIR="${BOLTZ_CACHE_DIR}" \
BOLTZ_ACCELERATOR="${BOLTZ_ACCELERATOR:-gpu}" \
BOLTZ_LOW_MEM="${BOLTZ_LOW_MEM:-1}" \
BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:-}" \
bash scripts/sherlock/run_boltz_wt_affinity.sh

PRED_DIR="${BOLTZ_OUT_DIR}/boltz_results_${TARGET_ID}/predictions/${TARGET_ID}"
echo ""
echo "Expected prediction directory:"
echo "  ${PRED_DIR}"
find "${PRED_DIR}" -maxdepth 1 -type f \( -name "*affinity*.json" -o -name "*confidence*.json" \) 2>/dev/null | sort || true
