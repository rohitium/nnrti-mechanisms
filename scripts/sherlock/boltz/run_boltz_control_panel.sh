#!/bin/bash
#
# Run Boltz-2 affinity replicates for a mutation panel on Sherlock.
#
# Default panel:
#   K103N Y181C V106A Y318F
#
# Usage:
#   bash scripts/sherlock/boltz/run_boltz_control_panel.sh
#
# Optional env vars:
#   BOLTZ_MUTATIONS="K103N Y181C V106A Y318F"
#   BOLTZ_REPLICATES=10
#   BOLTZ_SEED_START=1001
#   BOLTZ_PANEL_ROOT=/scratch/users/$USER/nnrti-mechanisms/results/boltz/control_panel
#   BOLTZ_LOW_MEM=1
#   BOLTZ_BASE_EXTRA_ARGS="--sampling_steps 50 --diffusion_samples 1"
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
cd "${PROJECT_ROOT}"

normalize_label() {
    local raw="$1"
    raw="${raw//+/_}"
    raw="${raw//\//_}"
    raw="${raw// /}"
    raw="${raw//-/_}"
    printf '%s' "${raw}" | tr -cd '[:alnum:]_'
}

if [ -d "/scratch/users/${USER}" ]; then
    DEFAULT_BASE="/scratch/users/${USER}/nnrti-mechanisms"
else
    DEFAULT_BASE="${PROJECT_ROOT}"
fi

BOLTZ_MUTATIONS="${BOLTZ_MUTATIONS:-K103N Y181C V106A Y318F}"
BOLTZ_REPLICATES="${BOLTZ_REPLICATES:-10}"
BOLTZ_SEED_START="${BOLTZ_SEED_START:-1001}"
BOLTZ_PANEL_ROOT="${BOLTZ_PANEL_ROOT:-${DEFAULT_BASE}/results/boltz/control_panel}"
BOLTZ_BASE_EXTRA_ARGS="${BOLTZ_BASE_EXTRA_ARGS:-${BOLTZ_EXTRA_ARGS:-}}"

read -r -a MUTATIONS_ARR <<< "${BOLTZ_MUTATIONS}"
if [ "${#MUTATIONS_ARR[@]}" -eq 0 ]; then
    echo "ERROR: no mutations provided in BOLTZ_MUTATIONS." >&2
    exit 2
fi

echo "=========================================="
echo "Boltz Control Panel Run"
echo "=========================================="
echo "Mutations:      ${BOLTZ_MUTATIONS}"
echo "Replicates:     ${BOLTZ_REPLICATES}"
echo "Seed start:     ${BOLTZ_SEED_START}"
echo "Panel root:     ${BOLTZ_PANEL_ROOT}"
echo "Low mem mode:   ${BOLTZ_LOW_MEM:-1}"
echo ""

for mutation in "${MUTATIONS_ARR[@]}"; do
    safe_label="$(normalize_label "${mutation}")"
    safe_label="$(printf '%s' "${safe_label}" | tr '[:lower:]' '[:upper:]')"
    if [ -z "${safe_label}" ]; then
        echo "Skipping invalid mutation label: ${mutation}" >&2
        continue
    fi

    target_id="${safe_label}_rt_dor_affinity"
    if [ "$(printf '%s' "${mutation}" | tr '[:lower:]' '[:upper:]')" = "WT" ]; then
        safe_label="wt"
        target_id="wt_rt_dor_affinity"
    fi

    run_root="${BOLTZ_PANEL_ROOT}/${safe_label}"
    work_dir="${run_root}/work"
    save_dir="${run_root}/replicates"
    mkdir -p "${save_dir}"

    echo "------------------------------------------"
    echo "Mutation: ${mutation} (${safe_label})"
    echo "Work dir: ${work_dir}"
    echo "Save dir: ${save_dir}"

    for ((rep = 1; rep <= BOLTZ_REPLICATES; rep++)); do
        seed=$((BOLTZ_SEED_START + rep - 1))
        if [ -n "${BOLTZ_BASE_EXTRA_ARGS}" ]; then
            extra_args="${BOLTZ_BASE_EXTRA_ARGS} --seed ${seed}"
        else
            extra_args="--seed ${seed}"
        fi

        echo ""
        echo "=== ${mutation} replicate ${rep}/${BOLTZ_REPLICATES} (seed=${seed}) ==="

        SHERLOCK_MODULES="${SHERLOCK_MODULES:-}" \
        BOLTZ_PYTHON="${BOLTZ_PYTHON:-python3}" \
        BOLTZ_OUT_DIR="${work_dir}" \
        BOLTZ_LOW_MEM="${BOLTZ_LOW_MEM:-1}" \
        BOLTZ_EXTRA_ARGS="${extra_args}" \
        bash scripts/sherlock/boltz/run_boltz_mutation_affinity.sh "${mutation}"

        pred_dir="${work_dir}/boltz_results_${target_id}/predictions/${target_id}"
        aff_json="${pred_dir}/affinity_${target_id}.json"
        conf_json="${pred_dir}/confidence_${target_id}_model_0.json"

        if [ ! -f "${aff_json}" ]; then
            echo "ERROR: missing affinity output after replicate ${rep}: ${aff_json}" >&2
            exit 1
        fi

        cp "${aff_json}" "${save_dir}/affinity_seed${seed}.json"
        if [ -f "${conf_json}" ]; then
            cp "${conf_json}" "${save_dir}/confidence_seed${seed}.json"
        fi
    done
done

echo ""
echo "Panel run complete. Replicate outputs:"
find "${BOLTZ_PANEL_ROOT}" -type f \( -name "affinity_seed*.json" -o -name "confidence_seed*.json" \) | sort
