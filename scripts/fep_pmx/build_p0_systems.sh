#!/bin/bash
#
# Build solvated GROMACS hybrid systems for P0 legs.
# Requires GROMACS + pmx on PATH (Sherlock login recommended).
#
# Usage:
#   bash scripts/fep_pmx/build_p0_systems.sh
#   REPLICATES=1 PHASES=holo bash scripts/fep_pmx/build_p0_systems.sh
#   FORCE=1 bash scripts/fep_pmx/build_p0_systems.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

REPLICATES="${REPLICATES:-1}"
# REPLICATES is a LIST of rep INDICES, not a count. Accept "1-3" range syntax, and
# loudly warn on a bare integer >1 (the footgun: REPLICATES=3 means rep 3 ONLY, not
# reps 1-3 — this once rebuilt only rep 3 of the charge legs). See git history.
_orig_reps="${REPLICATES}"
if [[ "${REPLICATES}" =~ ^[0-9]+-[0-9]+$ ]]; then
    REPLICATES="$(seq "${REPLICATES%-*}" "${REPLICATES#*-}")"
fi
if [[ "${_orig_reps}" =~ ^[0-9]+$ ]] && (( _orig_reps > 1 )); then
    echo "NOTE: REPLICATES='${_orig_reps}' = rep index ${_orig_reps} ONLY (not reps 1-${_orig_reps})." >&2
    echo "      For a range use REPLICATES='1-${_orig_reps}' or REPLICATES='1 2 ... ${_orig_reps}'." >&2
fi
echo "Reps to process: $(echo ${REPLICATES})"
PHASES="${PHASES:-holo apo}"
FORCE="${FORCE:-0}"
# Override with e.g. LEGS="wt_to_F227C wt_to_G190A" for other legs (P1).
read -r -a LEGS <<< "${LEGS:-wt_to_V106A wt_to_Y188L}"

# GROMACS (Sherlock)
if [[ -f scripts/sherlock/load_gromacs_module.sh ]]; then
    # shellcheck disable=SC1091
    source scripts/sherlock/load_gromacs_module.sh
fi

# pmx mutff (Sherlock venv or conda)
PMX_VENV="${PMX_VENV:-$HOME/.venvs/pmx}"
if command -v module >/dev/null 2>&1; then
    module load python/3.9.0 2>/dev/null || module load python/3.12.1 2>/dev/null || true
fi
if ! command -v pmx >/dev/null 2>&1 && [[ -f "${PMX_VENV}/bin/activate" ]]; then
    # shellcheck disable=SC1090
    source "${PMX_VENV}/bin/activate"
fi
PYTHON="${PYTHON:-python3}"

if [[ -z "${GMXLIB:-}" ]] && command -v "${PYTHON}" >/dev/null 2>&1; then
    export GMXLIB="$("${PYTHON}" - <<'PY'
import os
try:
    import pmx
    print(os.path.join(os.path.dirname(pmx.__file__), "data", "mutff"))
except ImportError:
    print("")
PY
)"
fi

if ! command -v gmx >/dev/null 2>&1; then
    echo "ERROR: gmx not found. On Sherlock: source scripts/sherlock/load_gromacs_module.sh" >&2
    exit 1
fi
if ! command -v pmx >/dev/null 2>&1; then
    echo "ERROR: pmx not found. On Sherlock: source ~/.venvs/pmx/bin/activate" >&2
    echo "       Or install: bash scripts/sherlock/setup_pmx_env.sh" >&2
    exit 1
fi
if [[ -z "${GMXLIB:-}" ]]; then
    echo "ERROR: GMXLIB not set (pmx mutff path required for pdb2gmx)" >&2
    exit 1
fi

echo "=========================================="
echo "P0 GROMACS solvated system build"
echo "=========================================="
echo "gmx:        $(command -v gmx)"
echo "pmx:        $(command -v pmx)"
echo "GMXLIB:     ${GMXLIB}"
echo "Replicates: ${REPLICATES}"
echo "Phases:     ${PHASES}"
echo "Force:      ${FORCE}"
echo ""

OK=0
SKIP=0
FAIL=0

for leg in "${LEGS[@]}"; do
    for phase in ${PHASES}; do
        for rep in ${REPLICATES}; do
            out="results/analysis/fep_pmx/legs/${leg}/${phase}/rep_$(printf '%02d' "${rep}")/gromacs_build/system.gro"
            if [[ -f "${out}" && "${FORCE}" != "1" ]]; then
                echo "SKIP ${leg} ${phase} rep${rep}"
                SKIP=$((SKIP + 1))
                continue
            fi
            echo "→ ${leg} ${phase} rep${rep}"
            if [[ "${FORCE}" == "1" ]]; then
                if "${PYTHON}" src/nnrti/fep/build_solvated_system.py \
                    --leg "${leg}" \
                    --phase "${phase}" \
                    --replicate "${rep}" \
                    --gmxlib "${GMXLIB}" \
                    --force; then
                    OK=$((OK + 1))
                else
                    echo "FAILED ${leg} ${phase} rep${rep}" >&2
                    FAIL=$((FAIL + 1))
                fi
            elif "${PYTHON}" src/nnrti/fep/build_solvated_system.py \
                --leg "${leg}" \
                --phase "${phase}" \
                --replicate "${rep}" \
                --gmxlib "${GMXLIB}"; then
                OK=$((OK + 1))
            else
                echo "FAILED ${leg} ${phase} rep${rep}" >&2
                FAIL=$((FAIL + 1))
            fi
        done
    done
done

echo ""
echo "Done: ${OK} built, ${SKIP} skipped, ${FAIL} failed"
[[ "${FAIL}" -eq 0 ]]
