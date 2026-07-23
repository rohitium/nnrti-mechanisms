#!/bin/bash
set -euo pipefail

# Install Perses hybrid-FEP dependencies into nnrti-prep without breaking MD prep.
# Pins numpy<2.4 because openmmtools/numba currently fail on numpy 2.4+.

ENV_NAME="${ENV_NAME:-nnrti-prep}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

echo "Updating conda env: ${ENV_NAME}"
conda install -n "${ENV_NAME}" -y -c conda-forge \
    "numpy<2.4" \
    perses=0.10.3 \
    openmmtools \
    openmmforcefields \
    ambertools

echo "Verifying imports (PYTHONNOUSERSITE=1 avoids ~/.local numpy shadowing conda)..."
PYTHONNOUSERSITE=1 conda run -n "${ENV_NAME}" python - <<'PY'
import numpy
import openmmtools
import perses
print("numpy", numpy.__version__)
print("openmmtools", openmmtools.__version__)
print("perses", perses.__version__)
PY

echo "Optional: create isolated nnrti-fep env from ${REPO_ROOT}/envs/nnrti-fep.yml"
echo "  conda env create -f ${REPO_ROOT}/envs/nnrti-fep.yml"
