#!/bin/bash
#
# Build pmx hybrid PDBs for P0 legs (wt_to_V106A, wt_to_Y188L) × holo/apo × reps.
#
# Usage:
#   bash scripts/fep_pmx/prepare_p0_hybrids.sh
#   REPLICATES=1 bash scripts/fep_pmx/prepare_p0_hybrids.sh   # rep 1 only
#   FORCE=1 bash scripts/fep_pmx/prepare_p0_hybrids.sh        # overwrite existing
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

REPLICATES="${REPLICATES:-1 2 3}"
FORCE="${FORCE:-0}"
PMX_CONDA_ENV="${PMX_CONDA_ENV:-pmx}"

if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${PMX_CONDA_ENV}"
fi

if ! command -v pmx >/dev/null 2>&1; then
    echo "ERROR: pmx not on PATH. Run: bash scripts/fep_pmx/setup_pmx_env.sh" >&2
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

LEGS=(wt_to_V106A wt_to_Y188L)
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
            if python scripts/fep_pmx/prepare_hybrid.py \
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
