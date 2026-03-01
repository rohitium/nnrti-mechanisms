#!/bin/bash
#
# Create/update a Sherlock conda env for Boltz-2 inference.
#
# Usage:
#   bash scripts/sherlock/setup_boltz_env.sh
#
# Optional env vars:
#   BOLTZ_ENV_NAME=boltz2
#   BOLTZ_PYTHON_VERSION=3.11
#   CONDA_HOME=$HOME/miniconda3
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

BOLTZ_ENV_NAME="${BOLTZ_ENV_NAME:-boltz2}"
BOLTZ_PYTHON_VERSION="${BOLTZ_PYTHON_VERSION:-3.11}"
CONDA_HOME="${CONDA_HOME:-$HOME/miniconda3}"

if [ -f "${CONDA_HOME}/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "${CONDA_HOME}/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    echo "ERROR: conda was not found. Set CONDA_HOME or load your conda module first." >&2
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "${BOLTZ_ENV_NAME}"; then
    echo "Creating conda env '${BOLTZ_ENV_NAME}' (python=${BOLTZ_PYTHON_VERSION})..."
    conda create -y -n "${BOLTZ_ENV_NAME}" "python=${BOLTZ_PYTHON_VERSION}"
else
    echo "Conda env '${BOLTZ_ENV_NAME}' already exists; updating packages in-place."
fi

conda activate "${BOLTZ_ENV_NAME}"

python -m pip install -U pip setuptools wheel
python -m pip install -U "boltz[cuda]"

echo ""
echo "Installed versions:"
python - <<'PY'
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
boltz predict --help >/dev/null

echo ""
echo "Boltz environment is ready."
echo "Activate with: conda activate ${BOLTZ_ENV_NAME}"

