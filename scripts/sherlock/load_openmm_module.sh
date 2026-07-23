#!/bin/bash
#
# Load Sherlock's OpenMM module stack (same as submit_md_batched.sh).
#
# Usage:
#   source scripts/sherlock/load_openmm_module.sh
#
# Override module string if needed:
#   SHERLOCK_OPENMM_MODULE="chemistry py-openmm/8.1.1_py312" \
#     source scripts/sherlock/load_openmm_module.sh
#

SHERLOCK_OPENMM_MODULE="${SHERLOCK_OPENMM_MODULE:-chemistry py-openmm/8.1.1_py312}"

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
            # shellcheck disable=SC1090
            source "${init}"
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
module load ${SHERLOCK_OPENMM_MODULE}

if ! python3 -c "import openmm" >/dev/null 2>&1; then
    echo "openmm is not importable after: module load ${SHERLOCK_OPENMM_MODULE}" >&2
    echo "python3: $(command -v python3 || echo missing)" >&2
    return 1 2>/dev/null || exit 1
fi
