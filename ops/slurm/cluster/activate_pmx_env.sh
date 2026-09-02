#!/bin/bash
#
# Activate pmx venv + GMXLIB on Sherlock.
#
# Call AFTER load_gromacs_module.sh (gromacs load purges python modules).
#
# Usage:
#   source ops/slurm/cluster/load_gromacs_module.sh
#   source ops/slurm/cluster/activate_pmx_env.sh
#
# Optional:
#   PMX_VENV=$HOME/.venvs/pmx
#   SHERLOCK_PYTHON_MODULE=python/3.9.0

PMX_VENV="${PMX_VENV:-$HOME/.venvs/pmx}"
SHERLOCK_PYTHON_MODULE="${SHERLOCK_PYTHON_MODULE:-python/3.9.0}"

_sherlock_init_modules() {
    if command -v module >/dev/null 2>&1; then
        return 0
    fi
    for init in \
        /etc/profile.d/modules.sh \
        /usr/share/lmod/lmod/init/bash \
        /share/software/user/open/lmod/lmod/init/bash \
        "${HOME}/.bashrc"
    do
        if [[ -f "${init}" ]]; then
            set +u
            # shellcheck disable=SC1090
            source "${init}"
            set -u
            if command -v module >/dev/null 2>&1; then
                return 0
            fi
        fi
    done
    return 1
}

_sherlock_init_modules || {
    echo "Could not initialize Lmod 'module' command." >&2
    return 1 2>/dev/null || exit 1
}

if command -v module >/dev/null 2>&1; then
    if ! module load "${SHERLOCK_PYTHON_MODULE}" 2>/dev/null; then
        module load python/3.12.1 2>/dev/null || {
            echo "ERROR: could not load ${SHERLOCK_PYTHON_MODULE} (required for pmx venv)" >&2
            return 1 2>/dev/null || exit 1
        }
    fi
fi

if [[ ! -f "${PMX_VENV}/bin/activate" ]]; then
    echo "ERROR: pmx venv not found at ${PMX_VENV}" >&2
    echo "  bash ops/slurm/cluster/setup_pmx_env.sh" >&2
    return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC1090
source "${PMX_VENV}/bin/activate"

PYTHON="${PYTHON:-python3}"
if [[ -z "${GMXLIB:-}" ]]; then
    export GMXLIB="$("${PYTHON}" - <<'PY'
import os
import pmx
print(os.path.join(os.path.dirname(pmx.__file__), "data", "mutff"))
PY
)"
fi

if [[ -z "${GMXLIB}" || ! -d "${GMXLIB}/amber14sbmut.ff" ]]; then
    echo "ERROR: GMXLIB not set to pmx mutff (got: ${GMXLIB:-<empty>})" >&2
    return 1 2>/dev/null || exit 1
fi

echo "[pmx] python: $(command -v "${PYTHON}") ($("${PYTHON}" --version 2>&1))"
echo "[pmx] GMXLIB=${GMXLIB}"
