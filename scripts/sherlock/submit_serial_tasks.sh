#!/bin/bash
# Submit MD tasks one-by-one from a manifest, waiting for each task to finish.
#
# Usage:
#   bash scripts/sherlock/submit_serial_tasks.sh [manifest]
#
# Optional env vars:
#   SHERLOCK_PARTITION   (default: gpu)
#   SHERLOCK_GRES        (default: gpu:1)
#   SHERLOCK_TIME        (default: 6:00:00)
#   SHERLOCK_MEM         (default: 16G)
#   SHERLOCK_QOS         (optional)
#   POLL_INTERVAL        (default: 60)
#   TASK_IDS             (optional, comma/range format, e.g. "0,1,4-7")
#   MD_HEATING_PS        (default: 25)
#   MD_PRODUCTION_NS     (default: 2.0)
#   MD_REPORT_INTERVAL   (default: 2000)
#   OPENMM_PLATFORM      (optional; CUDA/CPU/OpenCL)
#   MD_CHECKPOINT_INTERVAL (default: 5000)
#   MD_RESUME_FROM_CHECKPOINT (default: 1)

set -euo pipefail

MANIFEST_PATH="${1:-results/md_manifest.csv}"
if [ ! -f "$MANIFEST_PATH" ]; then
  echo "ERROR: Manifest not found: $MANIFEST_PATH" >&2
  exit 1
fi

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-6:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
POLL_INTERVAL="${POLL_INTERVAL:-60}"
TASK_IDS="${TASK_IDS:-}"

export MD_HEATING_PS="${MD_HEATING_PS:-25}"
export MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-2.0}"
export MD_REPORT_INTERVAL="${MD_REPORT_INTERVAL:-2000}"
export OPENMM_PLATFORM="${OPENMM_PLATFORM:-}"
export MD_CHECKPOINT_INTERVAL="${MD_CHECKPOINT_INTERVAL:-5000}"
export MD_RESUME_FROM_CHECKPOINT="${MD_RESUME_FROM_CHECKPOINT:-1}"

source /etc/profile.d/modules.sh 2>/dev/null || true
ml chemistry py-openmm/8.1.1_py312

log() { echo "$(date '+%H:%M:%S') [SERIAL] $*"; }

readarray -t TASK_LIST < <(
python3 - "$MANIFEST_PATH" "$TASK_IDS" <<'PY'
import csv
import sys

manifest = sys.argv[1]
selector = sys.argv[2].strip()

all_ids = []
with open(manifest, newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        all_ids.append(int(row["task_id"]))

if not all_ids:
    sys.exit(0)

if not selector:
    for tid in sorted(all_ids):
        print(tid)
    sys.exit(0)

selected = set()
for token in selector.split(","):
    token = token.strip()
    if not token:
        continue
    if "-" in token:
        lo_s, hi_s = token.split("-", 1)
        lo = int(lo_s)
        hi = int(hi_s)
        if hi < lo:
            lo, hi = hi, lo
        for tid in range(lo, hi + 1):
            selected.add(tid)
    else:
        selected.add(int(token))

all_set = set(all_ids)
for tid in sorted(selected):
    if tid in all_set:
        print(tid)
PY
)

if [ "${#TASK_LIST[@]}" -eq 0 ]; then
  log "No tasks selected."
  exit 0
fi

log "Selected ${#TASK_LIST[@]} task(s) from ${MANIFEST_PATH}"

for TASK_ID in "${TASK_LIST[@]}"; do
  STATUS_LINE=$(
    python3 - "$MANIFEST_PATH" "$TASK_ID" <<'PY'
import csv
import json
import os
import sys

manifest = sys.argv[1]
task_id = int(sys.argv[2])
row = None
with open(manifest, newline="") as handle:
    for r in csv.DictReader(handle):
        if int(r["task_id"]) == task_id:
            row = r
            break

if row is None:
    print("MISSING")
    sys.exit(0)

out = row["output_json"]
if os.path.exists(out):
    try:
        with open(out) as handle:
            payload = json.load(handle)
        if str(payload.get("status", "")).lower() == "ok":
            print("DONE")
            sys.exit(0)
    except Exception:
        pass

print("RUN")
PY
  )

  if [ "$STATUS_LINE" = "MISSING" ]; then
    log "Task ${TASK_ID}: missing from manifest, skipping."
    continue
  fi
  if [ "$STATUS_LINE" = "DONE" ]; then
    log "Task ${TASK_ID}: already completed, skipping."
    continue
  fi

  SBATCH_ARGS=(
    --parsable
    --array="${TASK_ID}-${TASK_ID}"
    --partition="${SHERLOCK_PARTITION}"
    --time="${SHERLOCK_TIME}"
    --mem="${SHERLOCK_MEM}"
  )
  if [ -n "${SHERLOCK_GRES}" ]; then
    SBATCH_ARGS+=(--gres="${SHERLOCK_GRES}")
  fi
  if [ -n "${SHERLOCK_QOS}" ]; then
    SBATCH_ARGS+=(--qos="${SHERLOCK_QOS}")
  fi

  log "Submitting task ${TASK_ID}..."
  JOB_ID=$(
    sbatch "${SBATCH_ARGS[@]}" \
      --export=ALL,MANIFEST_PATH="${MANIFEST_PATH}" \
      scripts/sherlock/submit_all_tasks.sh
  )
  log "Submitted job ${JOB_ID} for task ${TASK_ID}"

  while true; do
    ACTIVE=$(squeue -j "${JOB_ID}" -h | wc -l | tr -d '[:space:]')
    if [ "$ACTIVE" = "0" ]; then
      break
    fi
    sleep "${POLL_INTERVAL}"
  done

  FINAL=$(
    python3 - "$MANIFEST_PATH" "$TASK_ID" <<'PY'
import csv
import json
import os
import sys

manifest = sys.argv[1]
task_id = int(sys.argv[2])
row = None
with open(manifest, newline="") as handle:
    for r in csv.DictReader(handle):
        if int(r["task_id"]) == task_id:
            row = r
            break

if row is None:
    print("missing")
    sys.exit(0)

out = row["output_json"]
if not os.path.exists(out):
    print("missing_output")
    sys.exit(0)

try:
    with open(out) as handle:
        payload = json.load(handle)
    st = str(payload.get("status", "")).lower()
except Exception:
    st = "invalid_json"

print(st if st else "unknown")
PY
  )

  if [ "$FINAL" = "ok" ]; then
    log "Task ${TASK_ID} completed successfully."
  else
    log "Task ${TASK_ID} did not finish cleanly (status=${FINAL}). Check logs/md_${JOB_ID}_${TASK_ID}.err"
  fi
done

log "Serial submission loop complete."
