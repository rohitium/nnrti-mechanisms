#!/bin/bash
#
# Build the MM/GBSA scoring environment on Sherlock.
#
# Sherlock has no conda, so this layers a venv on top of the py-openmm module
# with --system-site-packages: OpenMM, numpy and scipy come from the module and
# only what the module lacks is installed. The venv lives in $GROUP_HOME
# because /etc/agents/AGENTS.md forbids environments on the purged filesystems.
#
# AmberTools is deliberately not installed. DOR's AM1-BCC charges are derived
# once off-cluster and shipped in data/ligands/dor_am1bcc.sdf, so the scoring
# job only has to read charges, never derive them.
#
set -euo pipefail

ENV_DIR="${ENV_DIR:-${GROUP_HOME}/envs/nnrti-mmgbsa}"

module load chemistry py-openmm/8.1.1_py312

if [ ! -x "${ENV_DIR}/bin/python" ]; then
    echo "creating ${ENV_DIR}"
    python3 -m venv --system-site-packages "${ENV_DIR}"
fi
# shellcheck disable=SC1091
source "${ENV_DIR}/bin/activate"

python -m pip install --upgrade pip wheel >/dev/null
python -m pip install \
    "openff-toolkit>=0.14" \
    "openmmforcefields>=0.12" \
    "MDAnalysis>=2.7" \
    "mdtraj>=1.9" \
    "pandas>=2.0"

echo "=== verification ==="
python - <<'PY'
import importlib
for m in ("openmm", "numpy", "scipy", "pandas", "MDAnalysis", "mdtraj",
          "openff.toolkit", "openmmforcefields"):
    mod = importlib.import_module(m)
    print(f"  {m:22s} {getattr(mod, '__version__', 'n/a')}")
from openff.toolkit.topology import Molecule
print("  openff Molecule import OK")
PY
echo "env ready: ${ENV_DIR}"
