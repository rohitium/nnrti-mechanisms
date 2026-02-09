#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash rsync_results.sh
# Optional:
#   RETRIES=8 SLEEP_SECONDS=15 SHERLOCK_USER=rsatija bash rsync_results.sh

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija"
  exit 1
fi

RETRIES="${RETRIES:-5}"
SLEEP_SECONDS="${SLEEP_SECONDS:-10}"
REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"

mkdir -p results/fep_runs logs

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

sync_with_retries \
  "${REMOTE_HOST}:${REMOTE_BASE}/results/fep_runs/" \
  "results/fep_runs/" \
  "results"

sync_with_retries \
  "${REMOTE_HOST}:${REMOTE_BASE}/logs/" \
  "logs/" \
  "logs"

echo "Done. You can now run: python -m src.main --collect-results"
