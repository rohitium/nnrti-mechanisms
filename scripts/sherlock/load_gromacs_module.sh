#!/bin/bash
#
# Load Sherlock GROMACS stack for pmx NEQ production (verified 2026-03).
#
# Usage:
#   source scripts/sherlock/load_gromacs_module.sh
#
# GPU MD jobs use gmx_cuda; prep (grompp, solvate) uses gmx.
# Do NOT load py-openmm in the same shell — CUDA stacks conflict.
#
# Override:
#   SHERLOCK_GROMACS_MODULE="chemistry gromacs/2023.1" source ...
#   GMX_MDRUN=gmx_cuda   # default for production
#   GMX=gmx              # default for grompp/prep

SHERLOCK_GROMACS_MODULE="${SHERLOCK_GROMACS_MODULE:-chemistry gromacs/2023.1}"
export GMX="${GMX:-gmx}"
export GMX_MDRUN="${GMX_MDRUN:-gmx_cuda}"

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

# shellcheck disable=SC2086
module purge 2>/dev/null || true
# shellcheck disable=SC2086
module load ${SHERLOCK_GROMACS_MODULE}

if ! command -v "${GMX}" >/dev/null 2>&1; then
    echo "gmx not found after: module load ${SHERLOCK_GROMACS_MODULE}" >&2
    return 1 2>/dev/null || exit 1
fi

if ! command -v "${GMX_MDRUN}" >/dev/null 2>&1; then
    echo "${GMX_MDRUN} not found — falling back to ${GMX} (CPU only)" >&2
    export GMX_MDRUN="${GMX}"
fi

# Quick sanity banner (non-fatal if grep misses on login node)
_gpu_line="$("${GMX_MDRUN}" --version 2>/dev/null | grep -i 'GPU support' || true)"
echo "[gromacs] module: ${SHERLOCK_GROMACS_MODULE}"
echo "[gromacs] prep: ${GMX}  mdrun: ${GMX_MDRUN}"
if [[ -n "${_gpu_line}" ]]; then
    echo "[gromacs] ${_gpu_line}"
fi
