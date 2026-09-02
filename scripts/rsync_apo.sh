#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_apo.sh push
#   SHERLOCK_USER=rsatija bash scripts/rsync_apo.sh pull
# Optional:
#   COMPLETE_ONLY=1    # pull only replicates that reached target steps
#   MD_PRODUCTION_NS=100.0
#   PARALLEL_JOBS=6    # concurrent rsync workers for push (default 6)

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
COMPLETE_ONLY="${COMPLETE_ONLY:-0}"
PARALLEL_JOBS="${PARALLEL_JOBS:-6}"
MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-100.0}"
REMOTE_BASE="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
REMOTE_HOST="${SHERLOCK_USER}@login.sherlock.stanford.edu"
LOCAL_APO_RUNS="results/md_runs/apo/"
LOCAL_MANIFEST="manifests/apo_md_manifest.csv"
REMOTE_APO_RUNS="${REMOTE_HOST}:${REMOTE_BASE}/results/md_runs/apo/"
REMOTE_MANIFEST="${REMOTE_HOST}:${REMOTE_BASE}/manifests/apo_md_manifest.csv"

# SSH ControlMaster socket — one Duo auth shared across all parallel rsync workers.
SSH_CTL="${TMPDIR:-/tmp}/nnrti_sherlock_ctl_${SHERLOCK_USER}.sock"

mkdir -p "${LOCAL_APO_RUNS}"
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
  local src="$1" dst="$2" label="$3"
  for ((i=1; i<=RETRIES; i++)); do
    echo "[${label}] attempt ${i}/${RETRIES}"
    if rsync "${RSYNC_FLAGS[@]}" "${src}" "${dst}"; then
      echo "[${label}] sync complete"; return 0
    fi
    [[ "${i}" -lt "${RETRIES}" ]] && { echo "[${label}] retrying in ${SLEEP_SECONDS}s..."; sleep "${SLEEP_SECONDS}"; }
  done
  echo "[${label}] failed after ${RETRIES} attempts"; return 1
}

sync_optional_with_retries() {
  sync_with_retries "$@" && return 0
  echo "[$3] optional sync skipped"; return 0
}

sync_parallel_subdirs() {
  local src_root="$1" dst_root="$2" label="$3"
  src_root="${src_root%/}"

  local subdirs=()
  while IFS= read -r d; do subdirs+=("$d"); done \
    < <(find "${src_root}" -mindepth 1 -maxdepth 1 -type d | sort)

  if [[ ${#subdirs[@]} -eq 0 ]]; then
    sync_with_retries "${src_root}/" "${dst_root}/" "${label}"; return
  fi

  echo "[${label}] parallel sync: ${#subdirs[@]} subdirs, ${PARALLEL_JOBS} workers"
  local pids=() failed=0 active=0

  for subdir in "${subdirs[@]}"; do
    local name; name="$(basename "${subdir}")"
    rsync "${RSYNC_FLAGS[@]}" "${subdir}/" "${dst_root}/${name}/" \
      >"${TMPDIR:-/tmp}/nnrti_rsync_${name}.log" 2>&1 &
    pids+=("$!"); (( active++ ))
    if [[ "${active}" -ge "${PARALLEL_JOBS}" ]]; then
      wait "${pids[0]}" || { echo "[${label}] worker failed: ${pids[0]}"; (( failed++ )); }
      pids=("${pids[@]:1}"); (( active-- ))
    fi
  done
  for pid in "${pids[@]}"; do
    wait "${pid}" || { echo "[${label}] worker failed: ${pid}"; (( failed++ )); }
  done

  [[ "${failed}" -gt 0 ]] && { echo "[${label}] ${failed} worker(s) failed"; return 1; }
  echo "[${label}] parallel sync complete"
}

build_complete_files_list() {
  local out_file="$1"
  ssh -S "${SSH_CTL}" "${REMOTE_HOST}" \
    "cd '${REMOTE_BASE}' && python3 - '${TARGET_STEPS}' <<'PY'
import glob, json, os, re, sys
from pathlib import Path

from nnrti.md.artifact_steps import infer_state_csv_path, reconcile_json_with_state_csv
target_steps = int(sys.argv[1])
json_pat = re.compile(r'.*_apo_rep[0-9]{2}\.json$')
paths = set()
for jp in glob.glob('results/md_runs/apo/*/rep_*/*.json'):
    if not json_pat.match(jp): continue
    json_path = Path(jp)
    try:
        with open(jp) as fh: payload = json.load(fh)
    except Exception: continue
    status = str(payload.get('status', '')).lower()
    reconciled = reconcile_json_with_state_csv(
        json_path=json_path,
        state_csv_path=infer_state_csv_path(json_path),
        write=True,
        target_steps=target_steps,
    )
    status = reconciled.status
    steps = reconciled.json_steps
    if status != 'ok' or steps < target_steps: continue
    rep_dir = os.path.dirname(jp)
    for root, _dirs, files in os.walk(rep_dir):
        for fn in files: paths.add(os.path.join(root, fn))
if os.path.exists('manifests/apo_md_manifest.csv'):
    paths.add('manifests/apo_md_manifest.csv')
for p in sorted(paths): print(p)
PY" > "${out_file}"
}

if [[ "${DIRECTION}" = "push" ]]; then
  sync_parallel_subdirs "${LOCAL_APO_RUNS}" \
    "${REMOTE_HOST}:${REMOTE_BASE}/results/md_runs/apo" "push apo_md_runs"

  if [[ -f "${LOCAL_MANIFEST}" ]]; then
    sync_with_retries "${LOCAL_MANIFEST}" "${REMOTE_MANIFEST}" "push apo_md_manifest.csv"
  else
    echo "[push apo_md_manifest.csv] not found locally, skipping"
  fi
else
  if [[ "${COMPLETE_ONLY}" = "1" ]]; then
    TMP_LIST="$(mktemp /tmp/nnrti_apo_complete_files.XXXXXX)"
    trap 'rm -f "${TMP_LIST}"' EXIT
    build_complete_files_list "${TMP_LIST}"
    if [[ -s "${TMP_LIST}" ]]; then
      echo "[pull complete-only apo_md_runs] syncing…"
      rsync "${RSYNC_FLAGS[@]}" --files-from="${TMP_LIST}" \
        "${REMOTE_HOST}:${REMOTE_BASE}/" "./" \
        || echo "[pull complete-only apo_md_runs] rsync returned non-zero"
    else
      echo "[pull apo] no complete replicates found at ${TARGET_STEPS} steps"
    fi
    rm -f "${TMP_LIST}"; trap - EXIT
  else
    sync_with_retries "${REMOTE_APO_RUNS}" "${LOCAL_APO_RUNS}" "pull apo_md_runs"
    sync_optional_with_retries "${REMOTE_MANIFEST}" "${LOCAL_MANIFEST}" "pull apo_md_manifest.csv"
  fi
fi

echo "Done (${DIRECTION})."
