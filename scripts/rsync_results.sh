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
#   COMPLETE_ONLY=1  # when pulling, transfer only replicate dirs that reached target steps
#   MD_PRODUCTION_NS=10.0  # target used with COMPLETE_ONLY=1
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
COMPLETE_ONLY="${COMPLETE_ONLY:-0}"
MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-10.0}"
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

TARGET_STEPS="$(python3 - <<PY
ns = float("${MD_PRODUCTION_NS}")
print(max(1, int(round((ns * 1_000_000.0) / 2.0))))
PY
)"

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

build_complete_files_list() {
  local out_file="$1"
  ssh "${REMOTE_HOST}" "cd '${REMOTE_BASE}' && python3 - '${TARGET_STEPS}' <<'PY'
import glob
import json
import os
import re
import sys

target_steps = int(sys.argv[1])
json_pat = re.compile(r'.*_rep[0-9]{2}\.json$')
paths = set()

for jp in glob.glob('results/md_runs/*/rep_*/*.json'):
    if not json_pat.match(jp):
        continue
    try:
        with open(jp) as fh:
            payload = json.load(fh)
    except Exception:
        continue
    status = str(payload.get('status', '')).lower()
    steps_raw = payload.get('md_production_steps_completed', payload.get('md_production_steps', 0))
    try:
        steps = int(steps_raw or 0)
    except Exception:
        steps = 0
    if status != 'ok' or steps < target_steps:
        continue

    rep_dir = os.path.dirname(jp)
    for root, _dirs, files in os.walk(rep_dir):
        for fn in files:
            paths.add(os.path.join(root, fn))

if os.path.exists('results/md_manifest.csv'):
    paths.add('results/md_manifest.csv')

for p in sorted(paths):
    print(p)
PY" > "${out_file}"
}

sync_complete_only_pull_with_retries() {
  local files_list="$1"
  local label="$2"

  if [[ ! -s "${files_list}" ]]; then
    echo "[${label}] no complete files found at target (${TARGET_STEPS} steps)"
    return 0
  fi

  for ((i=1; i<=RETRIES; i++)); do
    echo "[${label}] attempt ${i}/${RETRIES}"
    if rsync "${RSYNC_FLAGS[@]}" --files-from="${files_list}" "${REMOTE_HOST}:${REMOTE_BASE}/" "./"; then
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
  if [[ "${COMPLETE_ONLY}" = "1" ]]; then
    TMP_LIST="$(mktemp /tmp/nnrti_complete_files.XXXXXX)"
    trap 'rm -f "${TMP_LIST}"' EXIT
    build_complete_files_list "${TMP_LIST}"
    sync_complete_only_pull_with_retries "${TMP_LIST}" "pull complete-only md_runs"
    rm -f "${TMP_LIST}"
    trap - EXIT
  else
    sync_with_retries "${REMOTE_MD_RUNS}" "${LOCAL_MD_RUNS}" "pull md_runs"
    sync_optional_with_retries "${REMOTE_MANIFEST}" "${LOCAL_MANIFEST}" "pull md_manifest.csv"
  fi

  if [[ "${SYNC_LOGS}" = "1" ]]; then
    sync_optional_with_retries "${REMOTE_LOGS}" "${LOCAL_LOGS}" "pull logs"
  fi
fi

echo "Done (${DIRECTION})."
