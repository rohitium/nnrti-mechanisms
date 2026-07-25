#!/bin/bash
#
# Local pmx env for hybrid topology prep + NEQ analysis (Mac or Linux).
# Heavy GROMACS mdrun stays on Sherlock GPU nodes — not here.
#
# Usage:
#   bash scripts/fep_pmx/setup_pmx_env.sh
#
# Creates conda env "pmx" (Python 3.11) unless PMX_VENV is set for a venv path.
# Clones deGrootLab/pmx develop branch — NOT `pip install pmx` (wrong PyPI package).
#
# Optional:
#   PMX_REPO=$HOME/src/pmx
#   PMX_BRANCH=develop
#   PMX_CONDA_ENV=pmx
#   PMX_RECREATE=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

PMX_REPO="${PMX_REPO:-$HOME/src/pmx}"
PMX_BRANCH="${PMX_BRANCH:-develop}"
PMX_GIT_URL="${PMX_GIT_URL:-https://github.com/deGrootLab/pmx.git}"
PMX_CONDA_ENV="${PMX_CONDA_ENV:-pmx}"
PMX_VENV="${PMX_VENV:-}"
PMX_RECREATE="${PMX_RECREATE:-0}"
PMX_PYTHON="${PMX_PYTHON:-3.11}"

pick_python() {
    if [[ -n "${PMX_VENV}" ]]; then
        if [[ "${PMX_RECREATE}" == "1" ]] && [[ -d "${PMX_VENV}" ]]; then
            rm -rf "${PMX_VENV}"
        fi
        if [[ ! -d "${PMX_VENV}" ]]; then
            if command -v "python${PMX_PYTHON}" >/dev/null 2>&1; then
                "python${PMX_PYTHON}" -m venv "${PMX_VENV}"
            else
                python3 -m venv "${PMX_VENV}"
            fi
        fi
        # shellcheck disable=SC1090
        source "${PMX_VENV}/bin/activate"
        echo "[pmx] Using venv: ${PMX_VENV}"
        return
    fi

    if command -v conda >/dev/null 2>&1; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        if conda env list | awk '{print $1}' | grep -qx "${PMX_CONDA_ENV}"; then
            if [[ "${PMX_RECREATE}" == "1" ]]; then
                echo "[pmx] Removing conda env ${PMX_CONDA_ENV}"
                conda env remove -n "${PMX_CONDA_ENV}" -y
            fi
        fi
        if ! conda env list | awk '{print $1}' | grep -qx "${PMX_CONDA_ENV}"; then
            echo "[pmx] Creating conda env ${PMX_CONDA_ENV} (python=${PMX_PYTHON})"
            conda create -n "${PMX_CONDA_ENV}" "python=${PMX_PYTHON}" -y
        fi
        conda activate "${PMX_CONDA_ENV}"
        echo "[pmx] Using conda env: ${PMX_CONDA_ENV}"
        return
    fi

    echo "[pmx] ERROR: need conda or set PMX_VENV=/path/to/venv" >&2
    exit 1
}

pick_python

PY_VER="$(python - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
if [[ "${PY_VER}" == "3.12" ]] || [[ "${PY_VER}" == "3.13" ]]; then
    echo "[pmx] ERROR: Python ${PY_VER} incompatible with pmx (setuptools~=46)." >&2
    echo "[pmx] Use conda: PMX_PYTHON=3.11 bash $0" >&2
    exit 1
fi

echo "[pmx] Python: $(python --version)"
echo "[pmx] repo: ${PMX_REPO} (branch ${PMX_BRANCH})"

python -m pip install -U 'pip<24' wheel
python -m pip install 'setuptools==46.0.0'
# conda may leave a stale pth from a newer setuptools
find "$(python - <<'PY'
import site
print(site.getsitepackages()[0])
PY
)" -maxdepth 1 -name 'distutils-precedence.pth' -delete 2>/dev/null || true

python -m pip install numpy scipy matplotlib future

mkdir -p "$(dirname "${PMX_REPO}")"
if [[ ! -d "${PMX_REPO}/.git" ]]; then
    git clone "${PMX_GIT_URL}" "${PMX_REPO}"
fi

(
    cd "${PMX_REPO}"
    git fetch origin
    git checkout "${PMX_BRANCH}"
    git pull --ff-only origin "${PMX_BRANCH}" || true
    python setup.py build
    # setup.py install tries to resolve pmx on PyPI at the end; --skip-build avoids re-build
    python setup.py install --skip-build || true
)

if ! command -v pmx >/dev/null 2>&1; then
    echo "[pmx] ERROR: pmx CLI not on PATH after install" >&2
    exit 1
fi

GMXLIB_PATH="$(python - <<'PY'
import os
import pmx
print(os.path.join(os.path.dirname(pmx.__file__), "data", "mutff"))
PY
)"

echo ""
echo "[pmx] OK: $(pmx -h 2>&1 | head -1 || true)"
echo "[pmx] Add to ~/.zshrc (or export each session):"
echo "export GMXLIB=${GMXLIB_PATH}"
echo ""
if [[ -n "${PMX_VENV}" ]]; then
    echo "Activate: source ${PMX_VENV}/bin/activate"
else
    echo "Activate: conda activate ${PMX_CONDA_ENV}"
fi
