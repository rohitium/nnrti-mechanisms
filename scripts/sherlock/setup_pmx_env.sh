#!/bin/bash
#
# Install de Groot lab pmx (NOT the unrelated PyPI "pmx" package) for NEQ FEP.
#
# Usage (Sherlock login node — do NOT load gromacs in this shell):
#   module purge
#   module load python/3.12.1
#   bash scripts/sherlock/setup_pmx_env.sh
#
# Optional:
#   PMX_VENV=$HOME/.venvs/pmx
#   PMX_BRANCH=develop          # required for Python 3 + PRO mutations (P225H)
#   PMX_REPO=$HOME/src/pmx
#   SHERLOCK_MODULES="python/3.12.1"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

PMX_VENV="${PMX_VENV:-$HOME/.venvs/pmx}"
PMX_BRANCH="${PMX_BRANCH:-develop}"
PMX_REPO="${PMX_REPO:-$HOME/src/pmx}"
PMX_GIT_URL="${PMX_GIT_URL:-https://github.com/deGrootLab/pmx.git}"
SHERLOCK_MODULES="${SHERLOCK_MODULES:-python/3.12.1}"

if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null || true
    # shellcheck disable=SC2086
    module load ${SHERLOCK_MODULES}
fi

PYTHON="${PYTHON:-python3}"
if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "python3 not found; try: module load python/3.12.1" >&2
    exit 1
fi

echo "[pmx] Python: $(${PYTHON} --version 2>&1)"
echo "[pmx] venv: ${PMX_VENV}"
echo "[pmx] repo: ${PMX_REPO} (branch ${PMX_BRANCH})"

if [[ ! -d "${PMX_VENV}" ]]; then
    "${PYTHON}" -m venv "${PMX_VENV}"
fi
# shellcheck disable=SC1090
source "${PMX_VENV}/bin/activate"

pip install -U pip wheel
# pmx setup.py requires setuptools~=46.0.0; pip build isolation pulls latest setuptools → conflict
pip install 'setuptools~=46.0.0' numpy scipy matplotlib future

mkdir -p "$(dirname "${PMX_REPO}")"
if [[ ! -d "${PMX_REPO}/.git" ]]; then
    git clone "${PMX_GIT_URL}" "${PMX_REPO}"
fi

# git -C requires git >= 1.8.5; Sherlock default git is older
(
    cd "${PMX_REPO}"
    git fetch origin
    git checkout "${PMX_BRANCH}"
    git pull --ff-only origin "${PMX_BRANCH}" || true
)

pip install --no-build-isolation -U "${PMX_REPO}"

if ! command -v pmx >/dev/null 2>&1; then
    echo "pmx CLI not on PATH after install" >&2
    exit 1
fi

echo ""
echo "[pmx] Installed: $(pmx -h 2>&1 | head -1 || pmx --help 2>&1 | head -1)"
echo "[pmx] Set GMXLIB (once per session, or add to ~/.bashrc):"
pmx gmxlib 2>/dev/null || true
echo ""
echo "Activate with: source ${PMX_VENV}/bin/activate"
