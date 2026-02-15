#!/bin/bash
# ============================================================================
# DEPRECATED: This script is superseded by the batched Sherlock submit flow.
#
# New usage (on Sherlock):
#   bash scripts/sherlock/submit_md_batched.sh 6 12
#
# See README.md for full migration instructions.
# ============================================================================
echo "WARNING: orchestrate.sh is deprecated. Use submit_md_batched.sh instead." >&2
echo "  Run:  bash scripts/sherlock/submit_md_batched.sh 6 12" >&2
echo "" >&2

# Resumable mutation-by-mutation Sherlock workflow for explicit MD.
#
# Behavior:
#  1) For each mutation, verify local prep for all replicates; prep if missing.
#  2) Sync that mutation's prepared assets + manifest slice to Sherlock.
#  3) Submit only missing replicate tasks for that mutation.
#  4) Continue to next mutation without waiting for completion.
#  5) Poll running jobs, stream progress from md_state.csv, fetch partial results.
#  6) Run local collection/analysis incrementally as completed JSONs appear.
#
# Checkpoint/Resume:
#  - MD tasks write OpenMM checkpoint files (*.chk) and resume automatically.
#  - Orchestration state persists under results/orchestrate_state/.
#  - Re-running this script resumes from where it left off.
#
# Usage:
#   ./scripts/orchestrate.sh
#   ./scripts/orchestrate.sh --test
#   ./scripts/orchestrate.sh --collect-only

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SHERLOCK_USER="${SHERLOCK_USER:?Set SHERLOCK_USER (e.g., export SHERLOCK_USER=rsatija)}"
SHERLOCK_HOST="login.sherlock.stanford.edu"
SHERLOCK_DEST="${SHERLOCK_USER}@${SHERLOCK_HOST}"
SHERLOCK_DIR="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python}"
REPLICATES="${REPLICATES:-3}"
SEED="${SEED:-42}"
JITTER="${JITTER:-0.1}"
POLL_INTERVAL="${POLL_INTERVAL:-180}"

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-08:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-}"

MD_HEATING_PS="${MD_HEATING_PS:-25}"
MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-2.0}"
MD_REPORT_INTERVAL="${MD_REPORT_INTERVAL:-2000}"
OPENMM_PLATFORM="${OPENMM_PLATFORM:-}"

# New: explicit checkpoint controls for resume/extension
MD_CHECKPOINT_INTERVAL="${MD_CHECKPOINT_INTERVAL:-5000}"
MD_RESUME_FROM_CHECKPOINT="${MD_RESUME_FROM_CHECKPOINT:-1}"

# Incremental collection controls
COLLECT_METRIC_FRAME_STRIDE="${COLLECT_METRIC_FRAME_STRIDE:-5}"
COLLECT_METRIC_MAX_FRAMES="${COLLECT_METRIC_MAX_FRAMES:-200}"
COLLECT_MMGBSA_SNAPSHOTS="${COLLECT_MMGBSA_SNAPSHOTS:-100}"
COLLECT_MMGBSA_DISCARD_FRACTION="${COLLECT_MMGBSA_DISCARD_FRACTION:-0.25}"

STATE_DIR="${PROJECT_DIR}/results/orchestrate_state"
STATE_MANIFEST_DIR="${STATE_DIR}/manifests"
STATE_JOB_DIR="${STATE_DIR}/jobs"
STATE_MUT_DIR="${STATE_DIR}/mutations"
STATE_LAST_COLLECT_COUNT="${STATE_DIR}/last_collect_count.txt"
MASTER_MANIFEST="${STATE_DIR}/master_manifest.csv"

LOCAL_SOURCE_MANIFEST="${PROJECT_DIR}/results/md_manifest.csv"

SSH_SOCKET="/tmp/ssh-sherlock-${SHERLOCK_USER}"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=${SSH_SOCKET} -o ControlPersist=10800 -o ServerAliveInterval=60"

TEST_MODE=false
COLLECT_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --test) TEST_MODE=true ;;
        --collect-only) COLLECT_ONLY=true ;;
        --help|-h)
            sed -n '1,40p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown flag: $arg" >&2
            exit 1
            ;;
    esac
done

cd "$PROJECT_DIR"

log() { echo "$(date '+%H:%M:%S') [$1] $2"; }
fail() { echo "$(date '+%H:%M:%S') [FAIL] $1" >&2; }
warn() { echo "$(date '+%H:%M:%S') [WARN] $1"; }

mkdir -p "$STATE_DIR" "$STATE_MANIFEST_DIR" "$STATE_JOB_DIR" "$STATE_MUT_DIR"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
list_mutations() {
    "$PYTHON_BIN" - "$LOCAL_SOURCE_MANIFEST" <<'PY'
import csv
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
ordered = []
seen = set()

if manifest.exists():
    with manifest.open(newline="") as h:
        for row in csv.DictReader(h):
            m = str(row.get("mutation", "")).strip()
            if m and m not in seen:
                seen.add(m)
                ordered.append(m)
else:
    from src.analysis.susceptibility import load_dor_susceptibilities
    from pathlib import Path as _Path
    df = load_dor_susceptibilities(_Path("data/DRM-susceptibilities.csv.xlsx"), default_chain="A")
    for m in df["mutation"].astype(str).tolist():
        m = m.strip()
        if m and m not in seen:
            seen.add(m)
            ordered.append(m)

if "WT" not in seen:
    ordered = ["WT"] + ordered
else:
    ordered = ["WT"] + [m for m in ordered if m != "WT"]

for m in ordered:
    print(m)
PY
}

safe_label() {
    local mutation="$1"
    "$PYTHON_BIN" - "$mutation" <<'PY'
import sys
from src.utils import sanitize_label
print(sanitize_label(sys.argv[1]))
PY
}

mutation_manifest_path() {
    local safe="$1"
    echo "${STATE_MANIFEST_DIR}/${safe}.csv"
}

extract_mutation_manifest() {
    local source_manifest="$1"
    local mutation="$2"
    local out_manifest="$3"
    local replicates="$4"

    "$PYTHON_BIN" - "$source_manifest" "$mutation" "$out_manifest" "$replicates" <<'PY'
import csv
import sys
from pathlib import Path

src = Path(sys.argv[1])
mutation = sys.argv[2]
out = Path(sys.argv[3])
replicates = int(sys.argv[4])

if not src.exists():
    raise SystemExit(f"Source manifest missing: {src}")

rows = []
with src.open(newline="") as h:
    reader = csv.DictReader(h)
    fieldnames = reader.fieldnames
    for row in reader:
        if str(row.get("mutation", "")).strip() != mutation:
            continue
        rep = int(row.get("replicate", "0") or 0)
        if rep < 1 or rep > replicates:
            continue
        rows.append(row)

if not rows:
    raise SystemExit(f"No rows found for mutation '{mutation}' in {src}")

rows.sort(key=lambda r: int(r["replicate"]))
for i, row in enumerate(rows):
    row["task_id"] = str(i)

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as h:
    writer = csv.DictWriter(h, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {out} ({len(rows)} rows)")
PY
}

local_prep_complete() {
    local safe="$1"
    local rep
    for rep in $(seq 1 "$REPLICATES"); do
        local rep_tag
        rep_tag=$(printf "%02d" "$rep")
        local run_dir="${PROJECT_DIR}/results/md_runs/${safe}/rep_${rep_tag}"
        local min_pdb="${run_dir}/${safe}_minimized_rep${rep_tag}.pdb"
        local start_pdb="${run_dir}/assets/${safe}_md_rep${rep_tag}_start.pdb"
        local system_xml="${run_dir}/assets/${safe}_md_rep${rep_tag}_system.xml"
        if [ ! -f "$min_pdb" ] || [ ! -f "$start_pdb" ] || [ ! -f "$system_xml" ]; then
            return 1
        fi
    done
    return 0
}

prepare_mutation_if_needed() {
    local mutation="$1"
    local safe="$2"
    local mut_manifest
    mut_manifest="$(mutation_manifest_path "$safe")"

    if local_prep_complete "$safe"; then
        if [ -f "$mut_manifest" ]; then
            log PREP "${mutation}: local prep already complete."
            touch "${STATE_MUT_DIR}/${safe}.prep.ok"
            return 0
        fi
        if [ -f "$LOCAL_SOURCE_MANIFEST" ]; then
            log PREP "${mutation}: prep files exist; generating mutation manifest from results/md_manifest.csv."
            extract_mutation_manifest "$LOCAL_SOURCE_MANIFEST" "$mutation" "$mut_manifest" "$REPLICATES"
            touch "${STATE_MUT_DIR}/${safe}.prep.ok"
            return 0
        fi
    fi

    local tmp_manifest="${STATE_DIR}/tmp_${safe}_prep_manifest.csv"

    if [ "$mutation" = "WT" ]; then
        log PREP "WT prep missing. Running full local prep to materialize WT assets."
        (cd "$PROJECT_DIR" && "$PYTHON_BIN" -m src.main \
            --prepare-local-openmm-only \
            --replicates "$REPLICATES" \
            --seed "$SEED" \
            --jitter-angstrom "$JITTER" \
            --manifest "$tmp_manifest")
    else
        log PREP "${mutation}: local prep incomplete. Preparing this mutation now..."
        (cd "$PROJECT_DIR" && "$PYTHON_BIN" -m src.main \
            --prepare-local-openmm-only \
            --mutation "$mutation" \
            --replicates "$REPLICATES" \
            --seed "$SEED" \
            --jitter-angstrom "$JITTER" \
            --manifest "$tmp_manifest")
    fi

    if ! local_prep_complete "$safe"; then
        fail "${mutation}: prep command finished but required assets are still missing."
        return 1
    fi

    extract_mutation_manifest "$tmp_manifest" "$mutation" "$mut_manifest" "$REPLICATES"
    touch "${STATE_MUT_DIR}/${safe}.prep.ok"
    log PREP "${mutation}: prep complete."
}

refresh_master_manifest() {
    "$PYTHON_BIN" - "$STATE_MANIFEST_DIR" "$MASTER_MANIFEST" <<'PY'
import csv
import sys
from pathlib import Path

manifest_dir = Path(sys.argv[1])
out = Path(sys.argv[2])
files = sorted(manifest_dir.glob("*.csv"))

if not files:
    # Keep an empty marker file so callers can test existence.
    out.write_text("")
    raise SystemExit(0)

rows = []
fieldnames = None
for mf in files:
    with mf.open(newline="") as h:
        reader = csv.DictReader(h)
        if fieldnames is None:
            fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

rows.sort(key=lambda r: (str(r.get("mutation", "")), int(r.get("replicate", "0") or 0)))
for i, row in enumerate(rows):
    row["task_id"] = str(i)

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as h:
    writer = csv.DictWriter(h, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY
}

pending_task_ids() {
    local manifest="$1"
    "$PYTHON_BIN" - "$manifest" <<'PY'
import csv
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
if not manifest.exists() or manifest.stat().st_size == 0:
    print("")
    raise SystemExit(0)

pending = []
with manifest.open(newline="") as h:
    for row in csv.DictReader(h):
        task_id = int(row["task_id"])
        out = Path(str(row.get("output_json", "")).strip())
        done = False
        if out.exists():
            try:
                payload = json.loads(out.read_text())
                done = str(payload.get("status", "")).lower() == "ok"
            except Exception:
                done = False
        if not done:
            pending.append(str(task_id))

print(",".join(pending))
PY
}

completed_task_count() {
    local manifest="$1"
    "$PYTHON_BIN" - "$manifest" <<'PY'
import csv
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
if not manifest.exists() or manifest.stat().st_size == 0:
    print(0)
    raise SystemExit(0)

count = 0
with manifest.open(newline="") as h:
    for row in csv.DictReader(h):
        out = Path(str(row.get("output_json", "")).strip())
        if not out.exists():
            continue
        try:
            payload = json.loads(out.read_text())
            if str(payload.get("status", "")).lower() == "ok":
                count += 1
        except Exception:
            continue

print(count)
PY
}

mutation_progress_line() {
    local mut_manifest="$1"
    local target_steps
    target_steps=$("$PYTHON_BIN" - "$MD_PRODUCTION_NS" <<'PY'
import sys
ns=float(sys.argv[1])
steps=max(1,int(round((ns*1_000_000.0)/2.0)))
print(steps)
PY
)

    "$PYTHON_BIN" - "$mut_manifest" "$target_steps" <<'PY'
import csv
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
target_steps = int(sys.argv[2])

if not manifest.exists() or manifest.stat().st_size == 0:
    raise SystemExit(0)

rows = []
mutation_name = None
with manifest.open(newline="") as h:
    for row in csv.DictReader(h):
        rows.append(row)
        if mutation_name is None:
            mutation_name = str(row.get("mutation", "")).strip()

if not rows:
    raise SystemExit(0)

parts = []
for row in sorted(rows, key=lambda r: int(r.get("replicate", "0") or 0)):
    safe = str(row.get("safe_label", "")).strip()
    rep = int(row.get("replicate", "0") or 0)
    out = Path(str(row.get("output_json", "")).strip())

    if out.exists():
        try:
            payload = json.loads(out.read_text())
            if str(payload.get("status", "")).lower() == "ok":
                parts.append(f"r{rep}:done")
                continue
        except Exception:
            pass

    state_csv = out.parent / f"{safe}_rep{rep:02d}_md_state.csv"
    if state_csv.exists() and state_csv.stat().st_size > 0:
        last_step = None
        try:
            with state_csv.open() as h:
                for line in h:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith('"Step"'):
                        continue
                    try:
                        last_step = int(float(line.split(",", 1)[0]))
                    except Exception:
                        continue
        except Exception:
            last_step = None

        if last_step is not None and target_steps > 0:
            pct = min(100.0, 100.0 * float(last_step) / float(target_steps))
            parts.append(f"r{rep}:{pct:5.1f}%")
        else:
            parts.append(f"r{rep}:started")
    else:
        parts.append(f"r{rep}:queued")

print(f"{mutation_name}: " + ", ".join(parts))
PY
}

job_file_for_safe() {
    local safe="$1"
    echo "${STATE_JOB_DIR}/${safe}.jobid"
}

job_is_active_remote() {
    local job_id="$1"
    local out
    out=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" "squeue -j ${job_id} -h -o '%T' 2>/dev/null" || true)
    [ -n "${out}" ]
}

active_job_id_for_safe() {
    local safe="$1"
    local jf
    jf="$(job_file_for_safe "$safe")"
    if [ ! -f "$jf" ]; then
        echo ""
        return 0
    fi
    local job_id
    job_id=$(tr -d '[:space:]' < "$jf")
    if [ -z "$job_id" ]; then
        rm -f "$jf"
        echo ""
        return 0
    fi
    if job_is_active_remote "$job_id"; then
        echo "$job_id"
        return 0
    fi
    rm -f "$jf"
    echo ""
}

sync_codebase_once() {
    log SYNC "Syncing codebase to Sherlock (excluding bulky trajectories)..."
    rsync -avz -e "ssh $SSH_OPTS" \
        --exclude='.git' --exclude='.venv' --exclude='__pycache__' \
        --exclude='results/md_runs/**' \
        --exclude='results/orchestrate_state/**' \
        --exclude='*_md.dcd' --exclude='*_md_final.pdb' \
        "${PROJECT_DIR}/" "${SHERLOCK_DEST}:${SHERLOCK_DIR}/"
    ssh $SSH_OPTS "${SHERLOCK_DEST}" "mkdir -p '${SHERLOCK_DIR}/results/md_runs' '${SHERLOCK_DIR}/results/orchestrate_state/manifests' '${SHERLOCK_DIR}/logs'"
    log SYNC "Code sync complete."
}

sync_mutation_to_remote() {
    local safe="$1"
    local local_mut_dir="${PROJECT_DIR}/results/md_runs/${safe}"
    local local_manifest
    local_manifest="$(mutation_manifest_path "$safe")"
    local remote_manifest_rel="results/orchestrate_state/manifests/${safe}.csv"

    if [ -d "$local_mut_dir" ]; then
        rsync -avz -e "ssh $SSH_OPTS" \
            "${local_mut_dir}/" \
            "${SHERLOCK_DEST}:${SHERLOCK_DIR}/results/md_runs/${safe}/"
    fi

    rsync -avz -e "ssh $SSH_OPTS" \
        "$local_manifest" \
        "${SHERLOCK_DEST}:${SHERLOCK_DIR}/${remote_manifest_rel}"
}

submit_pending_tasks_for_mutation() {
    local mutation="$1"
    local safe="$2"
    local pending_ids="$3"

    if [ -z "$pending_ids" ]; then
        return 0
    fi

    local manifest_rel="results/orchestrate_state/manifests/${safe}.csv"

    sync_mutation_to_remote "$safe"

    local array_spec="$pending_ids"
    if [ -n "$ARRAY_CONCURRENCY" ] && [ "$ARRAY_CONCURRENCY" != "0" ]; then
        # Concurrency cap is only valid for range specs; keep comma-list unchanged.
        if [[ "$array_spec" == *"-"* ]] && [[ "$array_spec" != *","* ]]; then
            array_spec="${array_spec}%${ARRAY_CONCURRENCY}"
        fi
    fi

    log SUBMIT "${mutation}: submitting pending task IDs [${pending_ids}]"

    local job_id
    job_id=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" bash -s -- \
        "$SHERLOCK_DIR" "$PROJECT_DIR" "$manifest_rel" "$array_spec" \
        "$SHERLOCK_PARTITION" "$SHERLOCK_GRES" "$SHERLOCK_TIME" "$SHERLOCK_MEM" "$SHERLOCK_QOS" \
        "$MD_HEATING_PS" "$MD_PRODUCTION_NS" "$MD_REPORT_INTERVAL" "$OPENMM_PLATFORM" \
        "$MD_CHECKPOINT_INTERVAL" "$MD_RESUME_FROM_CHECKPOINT" <<'REMOTE'
set -euo pipefail
WORK_DIR="$1"
LOCAL_ROOT="$2"
MANIFEST_REL="$3"
ARRAY_SPEC="$4"
PARTITION="$5"
GRES="$6"
TIME_LIMIT="$7"
MEMORY="$8"
QOS="$9"
MD_HEATING_PS="${10}"
MD_PRODUCTION_NS="${11}"
MD_REPORT_INTERVAL="${12}"
OPENMM_PLATFORM="${13}"
MD_CHECKPOINT_INTERVAL="${14}"
MD_RESUME_FROM_CHECKPOINT="${15}"

cd "$WORK_DIR"
source /etc/profile.d/modules.sh 2>/dev/null || true
ml chemistry py-openmm/8.1.1_py312

python3 -m src.md.sherlock.rewrite_manifest_paths \
  --manifest "$MANIFEST_REL" \
  --from-root "$LOCAL_ROOT" \
  --to-root "$WORK_DIR" >/dev/null

mkdir -p logs

sbatch_args=(--parsable --array="$ARRAY_SPEC" --partition="$PARTITION" --time="$TIME_LIMIT" --mem="$MEMORY")
if [ -n "$GRES" ]; then
  sbatch_args+=(--gres="$GRES")
fi
if [ -n "$QOS" ]; then
  sbatch_args+=(--qos="$QOS")
fi

export MD_HEATING_PS MD_PRODUCTION_NS MD_REPORT_INTERVAL OPENMM_PLATFORM
export MD_CHECKPOINT_INTERVAL MD_RESUME_FROM_CHECKPOINT

sbatch "${sbatch_args[@]}" --export=ALL,MANIFEST_PATH="$MANIFEST_REL" scripts/sherlock/submit_all_tasks.sh
REMOTE
)

    job_id="${job_id%%;*}"
    if [ -z "$job_id" ]; then
        fail "${mutation}: sbatch did not return a job ID."
        return 1
    fi

    echo "$job_id" > "$(job_file_for_safe "$safe")"
    log SUBMIT "${mutation}: submitted job ${job_id}"
}

sync_remote_results_for_known_mutations() {
    local mf
    for mf in "$STATE_MANIFEST_DIR"/*.csv; do
        [ -e "$mf" ] || continue
        local safe
        safe="$(basename "$mf" .csv)"
        rsync -az -e "ssh $SSH_OPTS" \
            "${SHERLOCK_DEST}:${SHERLOCK_DIR}/results/md_runs/${safe}/" \
            "${PROJECT_DIR}/results/md_runs/${safe}/" 2>/dev/null || true
    done
}

rewrite_remote_paths_to_local() {
    if [ ! -f "$MASTER_MANIFEST" ] || [ ! -s "$MASTER_MANIFEST" ]; then
        return 0
    fi

    python3 -m src.md.sherlock.rewrite_manifest_paths \
        --manifest "$MASTER_MANIFEST" \
        --from-root "$SHERLOCK_DIR" \
        --to-root "$PROJECT_DIR" \
        --rewrite-jsons "${PROJECT_DIR}/results/md_runs" >/dev/null 2>&1 || true
}

run_incremental_collect_if_needed() {
    if [ ! -f "$MASTER_MANIFEST" ] || [ ! -s "$MASTER_MANIFEST" ]; then
        return 0
    fi

    local done_count
    done_count=$(completed_task_count "$MASTER_MANIFEST")

    local prev_count=0
    if [ -f "$STATE_LAST_COLLECT_COUNT" ]; then
        prev_count=$(tr -d '[:space:]' < "$STATE_LAST_COLLECT_COUNT" || echo 0)
        prev_count=${prev_count:-0}
    fi

    if [ "$done_count" -le "$prev_count" ]; then
        return 0
    fi

    log COLLECT "Detected ${done_count} completed task(s) (previously ${prev_count}); running incremental analysis..."
    if (cd "$PROJECT_DIR" && "$PYTHON_BIN" -m src.main \
        --collect-results \
        --manifest "$MASTER_MANIFEST" \
        --metric-frame-stride "$COLLECT_METRIC_FRAME_STRIDE" \
        --metric-max-frames "$COLLECT_METRIC_MAX_FRAMES" \
        --mmgbsa-snapshots "$COLLECT_MMGBSA_SNAPSHOTS" \
        --mmgbsa-discard-fraction "$COLLECT_MMGBSA_DISCARD_FRACTION"); then
        echo "$done_count" > "$STATE_LAST_COLLECT_COUNT"
        log COLLECT "Incremental analysis complete."
    else
        warn "Incremental collect failed; will retry on next poll."
    fi
}

log_progress_snapshot() {
    local mf
    for mf in "$STATE_MANIFEST_DIR"/*.csv; do
        [ -e "$mf" ] || continue
        local line
        line=$(mutation_progress_line "$mf" || true)
        if [ -n "$line" ]; then
            log PROGRESS "$line"
        fi
    done
}

poll_and_collect_once() {
    refresh_master_manifest

    local jf
    for jf in "$STATE_JOB_DIR"/*.jobid; do
        [ -e "$jf" ] || continue
        local safe job_id
        safe="$(basename "$jf" .jobid)"
        job_id=$(tr -d '[:space:]' < "$jf")
        if [ -z "$job_id" ]; then
            rm -f "$jf"
            continue
        fi

        local sq
        sq=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" "squeue -j ${job_id} -h -o '%T' 2>/dev/null" || true)
        if [ -z "$sq" ]; then
            log WAIT "${safe}: job ${job_id} no longer in queue (finished or failed)."
            rm -f "$jf"
        else
            local summary
            summary=$(echo "$sq" | sort | uniq -c | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g')
            log WAIT "${safe}: ${job_id} ${summary}"
        fi
    done

    sync_remote_results_for_known_mutations
    rewrite_remote_paths_to_local
    run_incremental_collect_if_needed
    log_progress_snapshot
}

all_mutations_complete() {
    local mutations=("$@")
    local mutation
    for mutation in "${mutations[@]}"; do
        local safe
        safe="$(safe_label "$mutation")"
        local manifest
        manifest="$(mutation_manifest_path "$safe")"

        if [ ! -f "$manifest" ] || [ ! -s "$manifest" ]; then
            return 1
        fi

        local pending
        pending=$(pending_task_ids "$manifest")
        if [ -n "$pending" ]; then
            return 1
        fi

        local active
        active=$(active_job_id_for_safe "$safe")
        if [ -n "$active" ]; then
            return 1
        fi
    done
    return 0
}

submit_if_needed() {
    local mutation="$1"
    local safe="$2"
    local manifest
    manifest="$(mutation_manifest_path "$safe")"

    local pending
    pending=$(pending_task_ids "$manifest")
    if [ -z "$pending" ]; then
        log SUBMIT "${mutation}: all replicates already complete."
        return 0
    fi

    local active
    active=$(active_job_id_for_safe "$safe")
    if [ -n "$active" ]; then
        log SUBMIT "${mutation}: existing active job ${active}; skipping new submission."
        return 0
    fi

    submit_pending_tasks_for_mutation "$mutation" "$safe" "$pending"
}

run_connectivity_test() {
    echo ""
    echo "=============================="
    echo " Sherlock connectivity test"
    echo "=============================="
    echo ""

    if ssh $SSH_OPTS "${SHERLOCK_DEST}" "echo 'SSH OK'; hostname"; then
        log PASS "SSH connection works"
    else
        fail "SSH connection failed"
        return 1
    fi

    if rsync -az -e "ssh $SSH_OPTS" "${PROJECT_DIR}/scripts/orchestrate.sh" "${SHERLOCK_DEST}:${SHERLOCK_DIR}/scripts/orchestrate.sh"; then
        log PASS "rsync works"
    else
        fail "rsync failed"
        return 1
    fi

    if ssh $SSH_OPTS "${SHERLOCK_DEST}" "squeue -u ${SHERLOCK_USER} --noheader | head -5; true"; then
        log PASS "SLURM query works"
    else
        fail "squeue failed"
        return 1
    fi

    echo ""
    echo "All connectivity checks passed."
    echo ""
}

# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------
if [ "$TEST_MODE" = true ]; then
    run_connectivity_test
    exit 0
fi

# Start control-master session once to avoid repeated Duo prompts.
ssh $SSH_OPTS "${SHERLOCK_DEST}" "echo 'SSH session ready on $(hostname)'" >/dev/null

sync_codebase_once

readarray -t MUTATIONS < <(list_mutations)
if [ "${#MUTATIONS[@]}" -eq 0 ]; then
    fail "No mutations discovered from manifest/workbook."
    exit 1
fi

log INFO "Mutations to process (${#MUTATIONS[@]}): ${MUTATIONS[*]}"

if [ "$COLLECT_ONLY" = true ]; then
    refresh_master_manifest
    poll_and_collect_once
    log DONE "Collect-only pass finished."
    ssh $SSH_OPTS -O exit "${SHERLOCK_DEST}" 2>/dev/null || true
    exit 0
fi

# Pass 1: iterate mutations, prep+submit one-by-one without waiting.
for mutation in "${MUTATIONS[@]}"; do
    safe="$(safe_label "$mutation")"
    prepare_mutation_if_needed "$mutation" "$safe"
    refresh_master_manifest
    submit_if_needed "$mutation" "$safe"
    poll_and_collect_once
    echo
    sleep 1
done

# Pass 2: keep polling; resubmit unfinished mutations as needed.
while true; do
    all_done=true

    for mutation in "${MUTATIONS[@]}"; do
        safe="$(safe_label "$mutation")"
        manifest="$(mutation_manifest_path "$safe")"

        if [ ! -f "$manifest" ] || [ ! -s "$manifest" ]; then
            all_done=false
            continue
        fi

        pending="$(pending_task_ids "$manifest")"
        if [ -n "$pending" ]; then
            all_done=false
            submit_if_needed "$mutation" "$safe"
        fi
    done

    poll_and_collect_once

    if [ "$all_done" = true ] && all_mutations_complete "${MUTATIONS[@]}"; then
        log DONE "All mutation/replicate tasks completed."
        break
    fi

    sleep "$POLL_INTERVAL"
done

# Final collect to ensure all plots/tables are up-to-date.
refresh_master_manifest
run_incremental_collect_if_needed

ssh $SSH_OPTS -O exit "${SHERLOCK_DEST}" 2>/dev/null || true
log DONE "Pipeline complete. Results in ${PROJECT_DIR}/results/"
