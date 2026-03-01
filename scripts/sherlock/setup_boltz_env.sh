#!/bin/bash
#
# Install/update Boltz-2 on Sherlock without creating a conda env.
#
# Usage:
#   bash scripts/sherlock/setup_boltz_env.sh
#
# Optional env vars:
#   SHERLOCK_MODULES="python/3.12.1"
#   BOLTZ_PYTHON=python3
#   BOLTZ_PIP_SPEC="boltz[cuda]"
#   BOLTZ_USER_INSTALL=1
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

SHERLOCK_MODULES="${SHERLOCK_MODULES:-}"
BOLTZ_PYTHON="${BOLTZ_PYTHON:-python3}"
BOLTZ_PIP_SPEC="${BOLTZ_PIP_SPEC:-boltz[cuda]}"
BOLTZ_USER_INSTALL="${BOLTZ_USER_INSTALL:-1}"

# Some Sherlock setups require explicit module loads for Python/CUDA.
if [ -n "${SHERLOCK_MODULES}" ] && command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    module load ${SHERLOCK_MODULES}
fi

# Ensure user-site scripts (e.g. ~/.local/bin/boltz) are reachable.
export PATH="$HOME/.local/bin:$PATH"

# Boltz/PyTorch require a modern but supported Python.
PY_VER="$("${BOLTZ_PYTHON}" - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
PY_MAJOR="${PY_VER%%.*}"
PY_MINOR="${PY_VER##*.}"
if [ "${PY_MAJOR}" -ne 3 ] || [ "${PY_MINOR}" -lt 10 ] || [ "${PY_MINOR}" -gt 13 ]; then
    echo "ERROR: ${BOLTZ_PYTHON} is Python ${PY_VER}; expected Python 3.10-3.13." >&2
    echo "Load a supported Python module first, for example:" >&2
    echo "  module avail python" >&2
    echo "  module load python/3.12.1" >&2
    echo "Then rerun this script (optionally set SHERLOCK_MODULES)." >&2
    exit 1
fi

PIP_ARGS=()
if [ "${BOLTZ_USER_INSTALL}" = "1" ]; then
    PIP_ARGS+=(--user)
fi

"${BOLTZ_PYTHON}" -m pip install -U "${PIP_ARGS[@]}" pip setuptools wheel
if ! "${BOLTZ_PYTHON}" -m pip install -U "${PIP_ARGS[@]}" "${BOLTZ_PIP_SPEC}"; then
    echo "Primary install failed for '${BOLTZ_PIP_SPEC}'. Retrying without extras..." >&2
    "${BOLTZ_PYTHON}" -m pip install -U "${PIP_ARGS[@]}" boltz
fi

echo ""
echo "Installed versions:"
"${BOLTZ_PYTHON}" - <<'PY'
import platform
import torch

print("python:", platform.python_version())
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu_count:", torch.cuda.device_count())
    print("gpu_name_0:", torch.cuda.get_device_name(0))
PY

# Smoke-check CLI import/arg parsing.
if command -v boltz >/dev/null 2>&1; then
    boltz predict --help >/dev/null
else
    "${BOLTZ_PYTHON}" -m boltz predict --help >/dev/null
fi

echo ""
echo "Boltz install is ready (no conda env created)."
echo "Using python: ${BOLTZ_PYTHON}"
