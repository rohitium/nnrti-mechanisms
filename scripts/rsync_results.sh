#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh push
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh pull
# Optional:
#   DIRECTION=push|pull
#   RETRIES=8 SLEEP_SECONDS=15
#   SYNC_LOGS=1      # also sync logs/
#
# Default direction is "push" so local artifacts can be restored to Sherlock
# in one command (results/md_runs + results/md_manifest.csv).

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija"
  exit 1
fi

DIRECTION="${1:-${DIRECTION:-push}}"
if [[ "${DIRECTION}" != "push" && "${DIRECTION}" != "pull" ]]; then
  echo "Invalid direction: ${DIRECTION} (expected: push or pull)"
  exit 1
fi

RETRIES="${RETRIES:-5}"
SLEEP_SECONDS="${SLEEP_SECONDS:-10}"
SYNC_LOGS="${SYNC_LOGS:-0}"
REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"
LOCAL_MD_RUNS="results/md_runs/"
LOCAL_MANIFEST="results/md_manifest.csv"
LOCAL_LOGS="logs/"
REMOTE_MD_RUNS="${REMOTE_HOST}:${REMOTE_BASE}/results/md_runs/"
REMOTE_MANIFEST="${REMOTE_HOST}:${REMOTE_BASE}/results/md_manifest.csv"
REMOTE_LOGS="${REMOTE_HOST}:${REMOTE_BASE}/logs/"

mkdir -p "${LOCAL_MD_RUNS}" "${LOCAL_LOGS}"

RSYNC_FLAGS=(-avzP --partial --inplace)
if rsync --help 2>/dev/null | grep -q -- "--append-verify"; then
  RSYNC_FLAGS+=(--append-verify)
elif rsync --help 2>/dev/null | grep -q -- "--append"; then
  RSYNC_FLAGS+=(--append)
fi

sync_with_retries() {
  local src="$1"
  local dst="$2"
  local label="$3"

  for ((i=1; i<=RETRIES; i++)); do
    echo "[${label}] attempt ${i}/${RETRIES}"
    if rsync "${RSYNC_FLAGS[@]}" "${src}" "${dst}"; then
      echo "[${label}] sync complete"
      return 0
    fi
    if [[ "${i}" -lt "${RETRIES}" ]]; then
      echo "[${label}] retrying in ${SLEEP_SECONDS}s..."
      sleep "${SLEEP_SECONDS}"
    fi
  done

  echo "[${label}] failed after ${RETRIES} attempts"
  return 1
}

sync_optional_with_retries() {
  local src="$1"
  local dst="$2"
  local label="$3"

  if sync_with_retries "${src}" "${dst}" "${label}"; then
    return 0
  fi
  echo "[${label}] optional sync skipped"
  return 0
}

if [[ "${DIRECTION}" = "push" ]]; then
  if [[ ! -d "${LOCAL_MD_RUNS}" ]]; then
    echo "Missing local directory: ${LOCAL_MD_RUNS}"
    exit 1
  fi

  sync_with_retries "${LOCAL_MD_RUNS}" "${REMOTE_MD_RUNS}" "push md_runs"

  if [[ -f "${LOCAL_MANIFEST}" ]]; then
    sync_with_retries "${LOCAL_MANIFEST}" "${REMOTE_MANIFEST}" "push md_manifest.csv"
  else
    echo "[push md_manifest.csv] not found locally, skipping"
  fi

  if [[ "${SYNC_LOGS}" = "1" ]]; then
    sync_optional_with_retries "${LOCAL_LOGS}" "${REMOTE_LOGS}" "push logs"
  fi
else
  sync_with_retries "${REMOTE_MD_RUNS}" "${LOCAL_MD_RUNS}" "pull md_runs"
  sync_optional_with_retries "${REMOTE_MANIFEST}" "${LOCAL_MANIFEST}" "pull md_manifest.csv"

  if [[ "${SYNC_LOGS}" = "1" ]]; then
    sync_optional_with_retries "${REMOTE_LOGS}" "${LOCAL_LOGS}" "pull logs"
  fi
fi

echo "Done (${DIRECTION})."
