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
#   SHERLOCK_MODULES="python/3.12.1"
#   BOLTZ_PYTHON=python3
#   BOLTZ_BIN=boltz
#   BOLTZ_INPUT_YAML=inputs/boltz/wt_rt_dor_affinity.yaml
#   BOLTZ_OUT_DIR=/scratch/users/$USER/nnrti-mechanisms/results/boltz/wt_affinity
#   BOLTZ_CACHE_DIR=/scratch/users/$USER/.boltz
#   BOLTZ_ACCELERATOR=gpu
#   BOLTZ_LOW_MEM=0
#   BOLTZ_EXTRA_ARGS="--sampling_steps 50 --diffusion_samples 1"
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

SHERLOCK_MODULES="${SHERLOCK_MODULES:-}"
BOLTZ_PYTHON="${BOLTZ_PYTHON:-python3}"
BOLTZ_BIN="${BOLTZ_BIN:-boltz}"
BOLTZ_INPUT_YAML="${BOLTZ_INPUT_YAML:-inputs/boltz/wt_rt_dor_affinity.yaml}"
BOLTZ_ACCELERATOR="${BOLTZ_ACCELERATOR:-gpu}"
BOLTZ_LOW_MEM="${BOLTZ_LOW_MEM:-0}"
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

# Some Sherlock setups require explicit module loads for Python/CUDA.
if [ -n "${SHERLOCK_MODULES}" ] && command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    module load ${SHERLOCK_MODULES}
fi

export PATH="$HOME/.local/bin:$PATH"

if [ ! -f "${BOLTZ_INPUT_YAML}" ]; then
    echo "Boltz input YAML not found; generating it from WT CIF..."
    "${BOLTZ_PYTHON}" scripts/sherlock/make_boltz_wt_input.py --output-yaml "${BOLTZ_INPUT_YAML}"
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
echo "Python:         ${BOLTZ_PYTHON}"
echo "Boltz bin:      ${BOLTZ_BIN}"
echo "Accelerator:    ${BOLTZ_ACCELERATOR}"
if [ -n "${BOLTZ_EXTRA_ARGS}" ]; then
    echo "Extra args:     ${BOLTZ_EXTRA_ARGS}"
fi
echo ""

if command -v "${BOLTZ_BIN}" >/dev/null 2>&1; then
    BOLTZ_CMD=("${BOLTZ_BIN}")
elif [ "${BOLTZ_BIN}" = "boltz" ]; then
    BOLTZ_CMD=("${BOLTZ_PYTHON}" -m boltz)
else
    echo "ERROR: Boltz binary not found: ${BOLTZ_BIN}" >&2
    exit 1
fi

CMD=(
    "${BOLTZ_CMD[@]}"
    predict
    "${BOLTZ_INPUT_YAML}"
    --out_dir
    "${BOLTZ_OUT_DIR}"
    --cache
    "${BOLTZ_CACHE_DIR}"
    --override
    --use_msa_server
)

# Compatibility shim:
# - Older Boltz CLIs required --affinity_predictor.
# - Newer CLIs infer affinity from YAML properties and do not expose this flag.
HELP_TEXT="$("${BOLTZ_CMD[@]}" predict --help 2>&1 || true)"
if echo "${HELP_TEXT}" | grep -q -- "--affinity_predictor"; then
    CMD+=(--affinity_predictor)
fi

if [ -n "${BOLTZ_ACCELERATOR}" ]; then
    CMD+=(--accelerator "${BOLTZ_ACCELERATOR}")
fi

if [ "${BOLTZ_LOW_MEM}" = "1" ]; then
    # Lower-memory configuration for older/smaller GPUs (e.g., 12 GB cards).
    CMD+=(
        --recycling_steps 1
        --sampling_steps 50
        --diffusion_samples 1
        --max_parallel_samples 1
        --max_msa_seqs 512
        --subsample_msa
        --num_subsampled_msa 256
        --sampling_steps_affinity 50
        --diffusion_samples_affinity 1
        --num_workers 0
        --preprocessing-threads 4
    )
    export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
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

set +e
time "${CMD[@]}"
RC=$?
set -e

if [ "${RC}" -ne 0 ]; then
    echo ""
    echo "Boltz exited with code ${RC}."
    echo "If you saw 'WARNING: ran out of memory, skipping batch', the GPU likely OOM'd during structure prediction."
    echo "In that case, retry with BOLTZ_LOW_MEM=1 or request a larger-memory GPU."
    exit "${RC}"
fi

echo ""
echo "Run complete. Affinity outputs (if present):"
find "${BOLTZ_OUT_DIR}" -type f -name "*affinity*.json" | sort || true
