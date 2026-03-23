#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh push
#   SHERLOCK_USER=rsatija bash scripts/rsync_results.sh pull
# Optional:
#   DIRECTION=push|pull
#   RETRIES=8 SLEEP_SECONDS=15
#   SYNC_LOGS=1        # also sync logs/
#   COMPLETE_ONLY=1    # when pulling, transfer only replicate dirs that reached target steps
#   MD_PRODUCTION_NS=100.0  # target used with COMPLETE_ONLY=1
#   PARALLEL_JOBS=6    # concurrent rsync workers for push (default 6)
#
# Push uses parallel per-mutation-dir workers sharing one SSH ControlMaster
# connection (single Duo auth).  Pull uses a single rsync stream.
#
# Default direction is "push" so local artifacts can be restored to Sherlock
# in one command (results/md_runs + manifests/md_manifest.csv).

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
PARALLEL_JOBS="${PARALLEL_JOBS:-6}"
MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-100.0}"
REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"
LOCAL_MD_RUNS="results/md_runs/"
LOCAL_MANIFEST="manifests/md_manifest.csv"
LOCAL_LOGS="logs/"
REMOTE_MD_RUNS="${REMOTE_HOST}:${REMOTE_BASE}/results/md_runs/"
REMOTE_MANIFEST="${REMOTE_HOST}:${REMOTE_BASE}/manifests/md_manifest.csv"
REMOTE_LOGS="${REMOTE_HOST}:${REMOTE_BASE}/logs/"

# SSH ControlMaster socket — one Duo auth shared across all parallel rsync workers.
SSH_CTL="${TMPDIR:-/tmp}/nnrti_sherlock_ctl_${SHERLOCK_USER}.sock"

mkdir -p "${LOCAL_MD_RUNS}" "${LOCAL_LOGS}"
mkdir -p "$(dirname "${LOCAL_MANIFEST}")"

# Open a persistent master connection if one isn't already running.
if ! ssh -S "${SSH_CTL}" -O check "${REMOTE_HOST}" 2>/dev/null; then
  echo "[ssh] Opening ControlMaster connection to ${REMOTE_HOST} (Duo auth required)…"
  ssh -M -S "${SSH_CTL}" -fN \
    -o ControlPersist=4h \
    -o ServerAliveInterval=60 \
    "${REMOTE_HOST}"
fi

RSYNC_SSH="ssh -S ${SSH_CTL}"
RSYNC_FLAGS=(-avzP --partial --inplace -e "${RSYNC_SSH}"
  --exclude='*.bak' --exclude='.DS_Store' --exclude='__pycache__')
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

# Parallel per-subdirectory sync: spawns up to PARALLEL_JOBS rsync workers at once.
# Faster than a single rsync stream when the tree has many large files (e.g. DCDs).
sync_parallel_subdirs() {
  local src_root="$1"   # local dir (trailing slash optional)
  local dst_root="$2"   # remote host:path base (no trailing slash)
  local label="$3"
  src_root="${src_root%/}"

  local subdirs=()
  while IFS= read -r d; do subdirs+=("$d"); done \
    < <(find "${src_root}" -mindepth 1 -maxdepth 1 -type d | sort)

  if [[ ${#subdirs[@]} -eq 0 ]]; then
    sync_with_retries "${src_root}/" "${dst_root}/" "${label}"
    return
  fi

  echo "[${label}] parallel sync: ${#subdirs[@]} subdirs, ${PARALLEL_JOBS} workers"

  local pids=() failed=0 active=0

  for subdir in "${subdirs[@]}"; do
    local name
    name="$(basename "${subdir}")"
    rsync "${RSYNC_FLAGS[@]}" "${subdir}/" "${dst_root}/${name}/" \
      >"${TMPDIR:-/tmp}/nnrti_rsync_${name}.log" 2>&1 &
    pids+=("$!")
    (( active++ ))

    # Throttle: wait for a slot when we hit the concurrency limit.
    if [[ "${active}" -ge "${PARALLEL_JOBS}" ]]; then
      wait "${pids[0]}" || { echo "[${label}] worker failed: ${pids[0]}"; (( failed++ )); }
      pids=("${pids[@]:1}")
      (( active-- ))
    fi
  done

  # Wait for remaining workers.
  for pid in "${pids[@]}"; do
    wait "${pid}" || { echo "[${label}] worker failed: ${pid}"; (( failed++ )); }
  done

  if [[ "${failed}" -gt 0 ]]; then
    echo "[${label}] ${failed} worker(s) failed — check /tmp/nnrti_rsync_*.log"
    return 1
  fi
  echo "[${label}] parallel sync complete"
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

from pathlib import Path

from src.md.artifact_steps import infer_state_csv_path, reconcile_json_with_state_csv

target_steps = int(sys.argv[1])
json_pat = re.compile(r'.*_rep[0-9]{2}\.json$')
paths = set()

for jp in glob.glob('results/md_runs/*/rep_*/*.json'):
    if not json_pat.match(jp):
        continue
    json_path = Path(jp)
    try:
        status = str(json.loads(json_path.read_text()).get('status', '')).lower()
    except Exception:
        continue
    reconciled = reconcile_json_with_state_csv(
        json_path=json_path,
        state_csv_path=infer_state_csv_path(json_path),
        write=True,
        target_steps=target_steps,
    )
    status = reconciled.status
    steps = reconciled.json_steps
    if status != 'ok' or steps < target_steps:
        continue

    rep_dir = os.path.dirname(jp)
    for root, _dirs, files in os.walk(rep_dir):
        for fn in files:
            paths.add(os.path.join(root, fn))

if os.path.exists('manifests/md_manifest.csv'):
    paths.add('manifests/md_manifest.csv')

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

  sync_parallel_subdirs "${LOCAL_MD_RUNS}" "${REMOTE_HOST}:${REMOTE_BASE}/results/md_runs" "push md_runs"

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
