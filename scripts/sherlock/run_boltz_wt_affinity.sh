#!/bin/bash
#
# Run a single Boltz-2 affinity prediction for WT RT + doravirine on Sherlock.
#
# Intended usage:
#   1) Get an interactive GPU allocation
#        bash scripts/sherlock/salloc_boltz_gpu.sh
#   2) On the allocated node
#        bash scripts/sherlock/run_boltz_wt_affinity.sh
#
# Optional env vars:
#   BOLTZ_ENV_NAME=boltz2
#   BOLTZ_INPUT_YAML=inputs/boltz/wt_rt_dor_affinity.yaml
#   BOLTZ_OUT_DIR=/scratch/users/$USER/nnrti-mechanisms/results/boltz/wt_affinity
#   BOLTZ_CACHE_DIR=/scratch/users/$USER/.boltz
#   BOLTZ_ACCELERATOR=gpu
#   BOLTZ_EXTRA_ARGS="--sampling_steps 50 --diffusion_samples 1"
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

BOLTZ_ENV_NAME="${BOLTZ_ENV_NAME:-boltz2}"
CONDA_HOME="${CONDA_HOME:-$HOME/miniconda3}"
BOLTZ_INPUT_YAML="${BOLTZ_INPUT_YAML:-inputs/boltz/wt_rt_dor_affinity.yaml}"
BOLTZ_ACCELERATOR="${BOLTZ_ACCELERATOR:-gpu}"
BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:-}"

if [ -d "/scratch/users/${USER}" ]; then
    DEFAULT_BASE="/scratch/users/${USER}/nnrti-mechanisms"
    DEFAULT_CACHE="/scratch/users/${USER}/.boltz"
else
    DEFAULT_BASE="${PROJECT_ROOT}"
    DEFAULT_CACHE="${PROJECT_ROOT}/.cache/boltz"
fi

BOLTZ_OUT_DIR="${BOLTZ_OUT_DIR:-${DEFAULT_BASE}/results/boltz/wt_affinity}"
BOLTZ_CACHE_DIR="${BOLTZ_CACHE_DIR:-${DEFAULT_CACHE}}"

mkdir -p "${BOLTZ_OUT_DIR}" "${BOLTZ_CACHE_DIR}" logs

if [ -f "${CONDA_HOME}/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "${CONDA_HOME}/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    echo "ERROR: conda was not found. Set CONDA_HOME or load your conda module first." >&2
    exit 1
fi

conda activate "${BOLTZ_ENV_NAME}"

if [ ! -f "${BOLTZ_INPUT_YAML}" ]; then
    echo "Boltz input YAML not found; generating it from WT CIF..."
    python scripts/sherlock/make_boltz_wt_input.py --output-yaml "${BOLTZ_INPUT_YAML}"
fi

export BOLTZ_CACHE="${BOLTZ_CACHE_DIR}"
export HF_HOME="${HF_HOME:-${BOLTZ_CACHE_DIR}/hf}"

echo "=========================================="
echo "Boltz WT Affinity Run"
echo "=========================================="
echo "Project root:   ${PROJECT_ROOT}"
echo "Input YAML:     ${BOLTZ_INPUT_YAML}"
echo "Output dir:     ${BOLTZ_OUT_DIR}"
echo "Cache dir:      ${BOLTZ_CACHE_DIR}"
echo "Conda env:      ${BOLTZ_ENV_NAME}"
echo "Accelerator:    ${BOLTZ_ACCELERATOR}"
if [ -n "${BOLTZ_EXTRA_ARGS}" ]; then
    echo "Extra args:     ${BOLTZ_EXTRA_ARGS}"
fi
echo ""

CMD=(
    boltz
    predict
    "${BOLTZ_INPUT_YAML}"
    --out_dir
    "${BOLTZ_OUT_DIR}"
    --cache
    "${BOLTZ_CACHE_DIR}"
    --override
    --use_msa_server
    --affinity_predictor
)

if [ -n "${BOLTZ_ACCELERATOR}" ]; then
    CMD+=(--accelerator "${BOLTZ_ACCELERATOR}")
fi

if [ -n "${BOLTZ_EXTRA_ARGS}" ]; then
    # Intentionally split user-provided flags into argv tokens.
    # shellcheck disable=SC2206
    EXTRA_ARR=(${BOLTZ_EXTRA_ARGS})
    CMD+=("${EXTRA_ARR[@]}")
fi

echo "Running:"
echo "  ${CMD[*]}"
echo ""

time "${CMD[@]}"

echo ""
echo "Run complete. Affinity outputs (if present):"
find "${BOLTZ_OUT_DIR}" -type f -name "*affinity*.json" | sort || true
