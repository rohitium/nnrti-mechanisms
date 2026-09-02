#!/bin/bash
#
# Build pmx hybrid PDBs for P0 legs (wt_to_V106A, wt_to_Y188L) × holo/apo × reps.
#
# Usage:
#   bash ops/slurm/fep/prepare_p0_hybrids.sh
#   REPLICATES=1 bash ops/slurm/fep/prepare_p0_hybrids.sh   # rep 1 only
#   FORCE=1 bash ops/slurm/fep/prepare_p0_hybrids.sh        # overwrite existing
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"
export PYTHONPATH="${PROJECT_ROOT:-$PWD}/src:${PYTHONPATH:-}"

REPLICATES="${REPLICATES:-1 2 3}"
# REPLICATES is a LIST of rep INDICES, not a count. Accept "1-3" range syntax, and
# loudly warn on a bare integer >1 (REPLICATES=3 means rep 3 ONLY, not reps 1-3).
_orig_reps="${REPLICATES}"
if [[ "${REPLICATES}" =~ ^[0-9]+-[0-9]+$ ]]; then
    REPLICATES="$(seq "${REPLICATES%-*}" "${REPLICATES#*-}")"
fi
if [[ "${_orig_reps}" =~ ^[0-9]+$ ]] && (( _orig_reps > 1 )); then
    echo "NOTE: REPLICATES='${_orig_reps}' = rep index ${_orig_reps} ONLY (not reps 1-${_orig_reps})." >&2
    echo "      For a range use REPLICATES='1-${_orig_reps}' or REPLICATES='1 2 ... ${_orig_reps}'." >&2
fi
echo "Reps to process: $(echo ${REPLICATES})"
FORCE="${FORCE:-0}"
PMX_CONDA_ENV="${PMX_CONDA_ENV:-pmx}"

if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${PMX_CONDA_ENV}"
fi

if ! command -v pmx >/dev/null 2>&1; then
    echo "ERROR: pmx not on PATH. Run: bash ops/slurm/fep/setup_pmx_env.sh" >&2
    exit 1
fi

if [[ -z "${GMXLIB:-}" ]]; then
    export GMXLIB="$(python - <<'PY'
import os
import pmx
print(os.path.join(os.path.dirname(pmx.__file__), "data", "mutff"))
PY
)"
fi

# Override with e.g. LEGS="wt_to_F227C wt_to_G190A" for other legs (P1).
read -r -a LEGS <<< "${LEGS:-wt_to_V106A wt_to_Y188L}"
PHASES=(holo apo)

echo "=========================================="
echo "P0 pmx hybrid prep"
echo "=========================================="
echo "GMXLIB:     ${GMXLIB}"
echo "Replicates: ${REPLICATES}"
echo "Force:      ${FORCE}"
echo ""

OK=0
SKIP=0
FAIL=0

for leg in "${LEGS[@]}"; do
    for phase in "${PHASES[@]}"; do
        for rep in ${REPLICATES}; do
            out="results/analysis/fep_pmx/legs/${leg}/${phase}/rep_$(printf '%02d' "${rep}")/hybrid.pdb"
            if [[ -f "${out}" && "${FORCE}" != "1" ]]; then
                echo "SKIP ${leg} ${phase} rep${rep} (exists)"
                SKIP=$((SKIP + 1))
                continue
            fi
            echo "→ ${leg} ${phase} rep${rep}"
            if python -m nnrti.fep.prepare_hybrid \
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
