#!/bin/bash
#
# Install/update Boltz-2 on Sherlock without creating a conda env.
#
# Usage:
#   bash scripts/sherlock/setup_boltz_env.sh
#
# Optional env vars:
#   SHERLOCK_MODULES="python/3.11 cuda/12.4"
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

PIP_ARGS=()
if [ "${BOLTZ_USER_INSTALL}" = "1" ]; then
    PIP_ARGS+=(--user)
fi

"${BOLTZ_PYTHON}" -m pip install -U "${PIP_ARGS[@]}" pip setuptools wheel
"${BOLTZ_PYTHON}" -m pip install -U "${PIP_ARGS[@]}" "${BOLTZ_PIP_SPEC}"

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
