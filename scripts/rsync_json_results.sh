#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_json_results.sh
#   SHERLOCK_USER=rsatija bash scripts/rsync_json_results.sh push
#   SHERLOCK_USER=rsatija bash scripts/rsync_json_results.sh pull
#
# Default direction is "push" so local JSON metadata can be restored to Sherlock
# in one command.

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

REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
REMOTE_MD_RUNS="${SHERLOCK_USER}@login.sherlock.stanford.edu:${REMOTE_BASE}/results/md_runs/"
REMOTE_MANIFEST="${SHERLOCK_USER}@login.sherlock.stanford.edu:${REMOTE_BASE}/results/md_manifest.csv"
LOCAL_MD_RUNS="results/md_runs/"
LOCAL_MANIFEST="results/md_manifest.csv"

mkdir -p "${LOCAL_MD_RUNS}"

JSON_FLAGS=(-avz)
if rsync --help 2>/dev/null | grep -q -- "--append-verify"; then
  JSON_FLAGS+=(--append-verify)
elif rsync --help 2>/dev/null | grep -q -- "--append"; then
  JSON_FLAGS+=(--append)
fi

sync_json_tree() {
  local src="$1"
  local dst="$2"
  rsync "${JSON_FLAGS[@]}" \
    --include='*/' \
    --include='*.json' \
    --exclude='*' \
    "${src}" \
    "${dst}"
}

if [[ "${DIRECTION}" = "push" ]]; then
  sync_json_tree "${LOCAL_MD_RUNS}" "${REMOTE_MD_RUNS}"
  if [[ -f "${LOCAL_MANIFEST}" ]]; then
    rsync "${JSON_FLAGS[@]}" "${LOCAL_MANIFEST}" "${REMOTE_MANIFEST}"
  fi
else
  sync_json_tree "${REMOTE_MD_RUNS}" "${LOCAL_MD_RUNS}"
  rsync "${JSON_FLAGS[@]}" "${REMOTE_MANIFEST}" "${LOCAL_MANIFEST}" || true
fi

echo "JSON sync complete (${DIRECTION})."
