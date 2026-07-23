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

if ! command -v module >/dev/null 2>&1; then
    echo "Sherlock 'module' command not found; are you on a compute node?" >&2
    return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC2086
module load ${SHERLOCK_OPENMM_MODULE}
