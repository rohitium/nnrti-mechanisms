#!/bin/bash
#
# Install de Groot lab pmx (NOT the unrelated PyPI "pmx" package) for NEQ FEP.
#
# Usage (Sherlock login node — do NOT load gromacs in this shell):
#   module purge
#   module load python/3.9.0
#   bash scripts/sherlock/setup_pmx_env.sh
#
# pmx pins setuptools~=46.0.0, which does not work on Python 3.12.
# Use Python 3.9 here; keep GROMACS jobs in a separate shell (load_gromacs_module.sh).
#
# Optional:
#   PMX_VENV=$HOME/.venvs/pmx
#   PMX_BRANCH=develop          # required for Python 3 + PRO mutations (P225H)
#   PMX_REPO=$HOME/src/pmx
#   PMX_RECREATE_VENV=1           # force fresh venv (needed after failed 3.12 attempt)
#   SHERLOCK_MODULES="python/3.9.0"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

PMX_VENV="${PMX_VENV:-$HOME/.venvs/pmx}"
PMX_BRANCH="${PMX_BRANCH:-develop}"
PMX_REPO="${PMX_REPO:-$HOME/src/pmx}"
PMX_GIT_URL="${PMX_GIT_URL:-https://github.com/deGrootLab/pmx.git}"
PMX_RECREATE_VENV="${PMX_RECREATE_VENV:-0}"
SHERLOCK_MODULES="${SHERLOCK_MODULES:-python/3.9.0}"

if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null || true
    # shellcheck disable=SC2086
    module load ${SHERLOCK_MODULES}
fi

PYTHON="${PYTHON:-python3}"
if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "python3 not found; try: module load python/3.9.0" >&2
    exit 1
fi

PY_VER="$("${PYTHON}" - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
if [[ "${PY_VER}" == "3.12" ]] || [[ "${PY_VER}" == "3.13" ]]; then
    echo "[pmx] ERROR: Python ${PY_VER} is incompatible with pmx (setuptools~=46)." >&2
    echo "[pmx] Use: module load python/3.9.0 && PMX_RECREATE_VENV=1 bash $0" >&2
    exit 1
fi

echo "[pmx] Python: $(${PYTHON} --version 2>&1)"
echo "[pmx] venv: ${PMX_VENV}"
echo "[pmx] repo: ${PMX_REPO} (branch ${PMX_BRANCH})"

if [[ "${PMX_RECREATE_VENV}" == "1" ]] && [[ -d "${PMX_VENV}" ]]; then
    echo "[pmx] Removing existing venv (PMX_RECREATE_VENV=1)"
    rm -rf "${PMX_VENV}"
fi

if [[ -d "${PMX_VENV}" ]] && [[ ! -x "${PMX_VENV}/bin/python" ]]; then
    rm -rf "${PMX_VENV}"
fi

if [[ -d "${PMX_VENV}" ]]; then
    VENV_PY="$("${PMX_VENV}/bin/python" - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
    if [[ "${VENV_PY}" != "${PY_VER}" ]]; then
        echo "[pmx] Recreating venv (was Python ${VENV_PY}, need ${PY_VER})"
        rm -rf "${PMX_VENV}"
    fi
fi

if [[ ! -d "${PMX_VENV}" ]]; then
    "${PYTHON}" -m venv "${PMX_VENV}"
fi
# shellcheck disable=SC1090
source "${PMX_VENV}/bin/activate"

pip install -U pip wheel
# pmx setup.py pins setuptools~=46.0.0
pip install 'setuptools~=46.0.0'
# Avoid compiling numpy/scipy on login nodes (old GCC); wheels only
pip install --only-binary=:all: numpy scipy matplotlib future

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
    # setup.py path avoids pip/pyproject + build-isolation setuptools conflicts
    python setup.py install
)

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
