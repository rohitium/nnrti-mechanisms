#!/usr/bin/env bash
set -euo pipefail

# Sync FEP pmx NEQ *light* results between Sherlock and the local Mac:
# per-system analysis (analysis.json, integ_{fwd,rev}.dat, results.txt,
# work_dist.png), targets/*/summary.json, panel CSVs + plots, and manifests.
# Heavy MD artifacts (trajectories, dgdl.xvg, tpr/cpt/gro/edr/log) are never
# transferred — mirrors the .gitignore split.
#
# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull   # Sherlock -> local (default)
#   SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh push   # local -> Sherlock
#   DRY=1 SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull   # preview only

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija" >&2
  exit 1
fi

DIRECTION="${1:-pull}"
if [[ "${DIRECTION}" != "pull" && "${DIRECTION}" != "push" ]]; then
  echo "Usage: SHERLOCK_USER=<user> bash scripts/rsync_fep_pmx.sh pull|push" >&2
  exit 1
fi

REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms-git"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"
SUBPATH="results/analysis/fep_pmx/"
LOCAL_DIR="${SUBPATH}"
REMOTE_DIR="${REMOTE_HOST}:${REMOTE_BASE}/${SUBPATH}"

# ControlMaster socket so Duo auth is entered once (shared with the jorgensen helper).
SSH_CTL="${TMPDIR:-/tmp}/nnrti_sherlock_ctl_${SHERLOCK_USER}.sock"
if ! ssh -S "${SSH_CTL}" -O check "${REMOTE_HOST}" 2>/dev/null; then
  echo "[ssh] Opening ControlMaster connection to ${REMOTE_HOST} (Duo auth required)…"
  if ! ssh -M -S "${SSH_CTL}" -fN \
    -o ControlPersist=4h \
    -o ServerAliveInterval=60 \
    "${REMOTE_HOST}"; then
    echo "SSH to ${REMOTE_HOST} failed; authenticate and retry." >&2
    exit 1
  fi
fi

# Light provenance only; first-match-wins filters, --exclude='*' blocks the rest.
FILTERS=(
  --include='*/'
  --include='**/analysis/**'
  --include='targets/***'
  --include='lambda_profiles/***'
  --include='crooks_overlap/***'
  --include='panel_*.csv'
  --include='panel_*.png'
  --include='neq_panel_manifest.csv'
  --include='**/residue_map.json'
  --include='**/neq_prepare.json'
  --include='**/neq_manifest.csv'
  --exclude='*'
)
RSYNC=(rsync -avzm --partial --update -e "ssh -S ${SSH_CTL}"
  --exclude='.DS_Store' --exclude='__pycache__'
  "${FILTERS[@]}")
if [[ -n "${DRY:-}" ]]; then
  RSYNC+=(--dry-run)
  echo "[dry-run] no files will be transferred"
fi

if [[ "${DIRECTION}" == "pull" ]]; then
  mkdir -p "${LOCAL_DIR}"
  echo "[pull] ${REMOTE_DIR} -> ${LOCAL_DIR}"
  "${RSYNC[@]}" "${REMOTE_DIR}" "${LOCAL_DIR}"
else
  echo "[push] ${LOCAL_DIR} -> ${REMOTE_DIR}"
  "${RSYNC[@]}" "${LOCAL_DIR}" "${REMOTE_DIR}"
fi
echo "[done] ${DIRECTION} fep_pmx light results"
