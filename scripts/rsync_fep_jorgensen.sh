#!/usr/bin/env bash
set -euo pipefail

# Sync one FEP Jorgensen leg between local Mac and Sherlock.
#
# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_fep_jorgensen.sh push V106A
#   SHERLOCK_USER=rsatija bash scripts/rsync_fep_jorgensen.sh pull V106A
#
# Push sends the full leg directory (inputs + holo + apo prep artifacts).
# Pull brings back window energy CSVs after Sherlock production runs.

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija"
  exit 1
fi

DIRECTION="${1:-push}"
MUTATION="${2:-V106A}"
if [[ "${DIRECTION}" != "push" && "${DIRECTION}" != "pull" ]]; then
  echo "Usage: SHERLOCK_USER=<user> bash scripts/rsync_fep_jorgensen.sh push|pull [MUTATION]"
  exit 1
fi

case "${MUTATION}" in
  V106A) LEG_ID="wt_to_V106A" ;;
  *)
    echo "Unsupported mutation label: ${MUTATION} (add a case mapping in this script)"
    exit 1
    ;;
esac

REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms-git"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"
LOCAL_LEG="results/analysis/fep_jorgensen/legs/${LEG_ID}/"
LOCAL_MANIFEST="results/analysis/fep_jorgensen/worker_manifest_v106a.csv"
REMOTE_LEG="${REMOTE_HOST}:${REMOTE_BASE}/results/analysis/fep_jorgensen/legs/${LEG_ID}/"
REMOTE_MANIFEST="${REMOTE_HOST}:${REMOTE_BASE}/results/analysis/fep_jorgensen/worker_manifest_v106a.csv"

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

RSYNC_SSH="ssh -S ${SSH_CTL}"
RSYNC_FLAGS=(-av --progress --partial --inplace -e "${RSYNC_SSH}"
  --exclude='.DS_Store' --exclude='__pycache__' --exclude='holo/multistate.nc')

if [[ "${DIRECTION}" == "push" ]]; then
  if [[ ! -f "${LOCAL_LEG}apo/hybrid_system.xml" ]]; then
    echo "Missing local apo hybrid: ${LOCAL_LEG}apo/hybrid_system.xml" >&2
    exit 1
  fi
  echo "[push] ${LOCAL_LEG} -> ${REMOTE_LEG}"
  rsync "${RSYNC_FLAGS[@]}" "${LOCAL_LEG}" "${REMOTE_LEG}"
  if [[ -f "${LOCAL_MANIFEST}" ]]; then
    echo "[push] ${LOCAL_MANIFEST} -> ${REMOTE_MANIFEST}"
    rsync "${RSYNC_FLAGS[@]}" "${LOCAL_MANIFEST}" "${REMOTE_MANIFEST}"
  fi
else
  echo "[pull] ${REMOTE_LEG} -> ${LOCAL_LEG}"
  rsync "${RSYNC_FLAGS[@]}" "${REMOTE_LEG}" "${LOCAL_LEG}"
fi

echo "Done (${DIRECTION} ${MUTATION})."
