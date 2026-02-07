#!/bin/bash
# End-to-end FEP pipeline: local prep → Sherlock → collect results
#
# Prerequisites:
#   - conda activate nnrti-prep  (for local OpenMM prep & analysis)
#   - SSH key auth to Sherlock    (ssh-copy-id $SHERLOCK_USER@login.sherlock.stanford.edu)
#
# Usage:
#   ./scripts/orchestrate.sh --test          # Test SSH/rsync connectivity
#   ./scripts/orchestrate.sh                 # Full pipeline
#   ./scripts/orchestrate.sh --skip-prep     # Skip local prep (already done)
#   ./scripts/orchestrate.sh --collect-only  # Just rsync results back + analyze
#
# Environment variables:
#   SHERLOCK_USER   (required) Your SUNet ID
#   REPLICATES      Number of replicates (default: 3)
#   SEED            Base random seed (default: 42)
#   POLL_INTERVAL   Seconds between job status checks (default: 300)

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SHERLOCK_USER="${SHERLOCK_USER:?Set SHERLOCK_USER to your SUNet ID (export SHERLOCK_USER=rsatija)}"
SHERLOCK_HOST="login.sherlock.stanford.edu"
SHERLOCK_DEST="${SHERLOCK_USER}@${SHERLOCK_HOST}"
SHERLOCK_DIR="/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

REPLICATES="${REPLICATES:-3}"
SEED="${SEED:-42}"
JITTER="0.1"
POLL_INTERVAL="${POLL_INTERVAL:-300}"

# SSH multiplexing — authenticate once (Duo), reuse for all subsequent calls.
# ControlPersist=10800 keeps the socket alive for 3 hours without activity,
# which covers the full cluster wait time.
SSH_SOCKET="/tmp/ssh-sherlock-${SHERLOCK_USER}"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=${SSH_SOCKET} -o ControlPersist=10800 -o ServerAliveInterval=60"

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
SKIP_PREP=false
COLLECT_ONLY=false
TEST_MODE=false

for arg in "$@"; do
    case $arg in
        --skip-prep)    SKIP_PREP=true ;;
        --collect-only) COLLECT_ONLY=true; SKIP_PREP=true ;;
        --test)         TEST_MODE=true; SKIP_PREP=true ;;
        --help|-h)
            head -18 "$0" | tail -16
            exit 0
            ;;
        *) echo "Unknown flag: $arg"; exit 1 ;;
    esac
done

log() { echo "$(date '+%H:%M:%S') [$1] $2"; }
pass() { echo "$(date '+%H:%M:%S') [PASS] $1"; }
fail() { echo "$(date '+%H:%M:%S') [FAIL] $1"; }

# ---------------------------------------------------------------------------
# --test: Verify SSH, rsync, and remote commands work
# ---------------------------------------------------------------------------
if [ "$TEST_MODE" = true ]; then
    echo ""
    echo "=============================="
    echo " Sherlock connectivity test"
    echo "=============================="
    echo ""
    echo "Target: ${SHERLOCK_DEST}:${SHERLOCK_DIR}"
    echo ""

    # Test 1: SSH connection
    echo "--- Test 1: SSH connection ---"
    echo "  (You may be prompted for Duo 2FA — this authenticates the"
    echo "   multiplexed socket so later commands won't re-prompt.)"
    echo ""
    if ssh $SSH_OPTS "${SHERLOCK_DEST}" "echo 'SSH OK'; hostname" 2>&1; then
        pass "SSH connection works"
    else
        fail "SSH connection failed"
        echo ""
        echo "Troubleshooting:"
        echo "  1. Check VPN is connected (if off-campus)"
        echo "  2. Test manually: ssh ${SHERLOCK_DEST}"
        echo "  3. Set up SSH keys: ssh-keygen && ssh-copy-id ${SHERLOCK_DEST}"
        exit 1
    fi
    echo ""

    # Test 2: Scratch directory exists
    echo "--- Test 2: Scratch directory ---"
    if ssh $SSH_OPTS "${SHERLOCK_DEST}" "ls -d /scratch/users/${SHERLOCK_USER}" 2>&1; then
        pass "Scratch directory exists"
    else
        fail "Scratch directory not found"
        echo "  Create it on Sherlock: mkdir -p ${SHERLOCK_DIR}"
        exit 1
    fi
    echo ""

    # Test 3: rsync a small test file
    echo "--- Test 3: rsync transfer ---"
    TESTFILE=$(mktemp "${PROJECT_DIR}/.rsync_test_XXXXXX")
    echo "orchestrate test $(date)" > "$TESTFILE"
    TESTBASE=$(basename "$TESTFILE")
    if rsync -avz -e "ssh $SSH_OPTS" "$TESTFILE" \
        "${SHERLOCK_DEST}:${SHERLOCK_DIR}/${TESTBASE}" 2>&1; then
        pass "rsync to Sherlock works"
        # Cleanup
        ssh $SSH_OPTS "${SHERLOCK_DEST}" "rm -f ${SHERLOCK_DIR}/${TESTBASE}" 2>/dev/null
    else
        fail "rsync failed"
        echo "  Make sure ${SHERLOCK_DIR} exists on Sherlock"
    fi
    rm -f "$TESTFILE"
    echo ""

    # Test 4: Remote module loading + python
    echo "--- Test 4: Remote OpenMM module ---"
    if ssh $SSH_OPTS "${SHERLOCK_DEST}" bash -c "'
        source /etc/profile.d/modules.sh 2>/dev/null || true
        ml chemistry py-openmm/8.1.1_py312 2>/dev/null
        python3 -c \"import openmm; print(f\\\"OpenMM {openmm.__version__}\\\")\"
    '" 2>&1; then
        pass "Remote OpenMM module works"
    else
        fail "Could not load OpenMM module on Sherlock"
    fi
    echo ""

    # Test 5: SLURM commands
    echo "--- Test 5: SLURM access ---"
    if ssh $SSH_OPTS "${SHERLOCK_DEST}" "squeue -u ${SHERLOCK_USER} --noheader | head -5; echo 'squeue OK'" 2>&1; then
        pass "SLURM squeue works"
    else
        fail "squeue failed"
    fi
    echo ""

    # Test 6: Check manifest
    echo "--- Test 6: Local manifest ---"
    MANIFEST="${PROJECT_DIR}/results/fep_manifest.csv"
    if [ -f "$MANIFEST" ]; then
        N_TASKS=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
        N_PAIRS=$((N_TASKS / 2))
        pass "Manifest found: ${N_TASKS} tasks (${N_PAIRS} array jobs)"

        echo "  Mutations in manifest:"
        tail -n +2 "$MANIFEST" | cut -d',' -f3 | sort -u | while read -r mut; do
            echo "    - $mut"
        done
    else
        fail "No manifest at ${MANIFEST}"
        echo "  Run local preparation first"
    fi
    echo ""

    echo "=============================="
    echo " All tests passed!"
    echo "=============================="
    echo ""
    echo "The SSH socket is now open and will persist for 3 hours."
    echo "Subsequent runs of this script will reuse it (no re-auth)."
    echo ""
    echo "Ready to run:"
    echo "  ./scripts/orchestrate.sh --skip-prep"
    echo ""

    exit 0
fi

# ---------------------------------------------------------------------------
# Phase 1: Local preparation
# ---------------------------------------------------------------------------
if [ "$SKIP_PREP" = false ]; then
    log PREP "Preparing OpenMM assets locally (${REPLICATES} replicates, seed ${SEED})..."
    cd "$PROJECT_DIR"
    python -m src.main \
        --prepare-local-openmm-only \
        --replicates "$REPLICATES" \
        --seed "$SEED" \
        --jitter-angstrom "$JITTER"
    log PREP "Local preparation complete."
else
    log PREP "Skipped (--skip-prep)."
fi

if [ "$COLLECT_ONLY" = true ]; then
    # Jump straight to Phase 5
    log COLLECT "Transferring results from Sherlock..."
    rsync -avz -e "ssh $SSH_OPTS" \
        "${SHERLOCK_DEST}:${SHERLOCK_DIR}/results/fep_runs/" \
        "${PROJECT_DIR}/results/fep_runs/"

    log COLLECT "Running result collection..."
    cd "$PROJECT_DIR"
    python -m src.main --collect-results

    log DONE "Results collected. Check results/ for output files."
    exit 0
fi

# ---------------------------------------------------------------------------
# Phase 2: Transfer to Sherlock
# ---------------------------------------------------------------------------
log SYNC "Transferring project to Sherlock..."
rsync -avz -e "ssh $SSH_OPTS" \
    --exclude='.venv' --exclude='.git' --exclude='__pycache__' \
    "${PROJECT_DIR}/" "${SHERLOCK_DEST}:${SHERLOCK_DIR}/"
log SYNC "Transfer complete."

# ---------------------------------------------------------------------------
# Phase 3: Rewrite manifest + submit SLURM jobs
# ---------------------------------------------------------------------------
log SUBMIT "Rewriting manifest paths and submitting jobs on Sherlock..."

# Count tasks in manifest to compute array range
N_TASKS=$(tail -n +2 "${PROJECT_DIR}/results/fep_manifest.csv" | wc -l | tr -d ' ')
N_PAIRS=$((N_TASKS / 2))
ARRAY_MAX=$((N_PAIRS - 1))
log SUBMIT "Manifest has ${N_TASKS} tasks -> ${N_PAIRS} array jobs (0-${ARRAY_MAX})"

# Run remote commands via SSH
JOB_ID=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" bash -s "$SHERLOCK_DIR" "$ARRAY_MAX" <<'REMOTE_SCRIPT'
set -euo pipefail
WORK_DIR="$1"
ARRAY_MAX="$2"
cd "$WORK_DIR"

# Rewrite manifest paths (local Mac → Sherlock scratch)
source /etc/profile.d/modules.sh 2>/dev/null || true
ml chemistry py-openmm/8.1.1_py312
python3 scripts/sherlock/rewrite_manifest_paths.py >&2

# Update array range in submit script
sed -i "s/^#SBATCH --array=.*/#SBATCH --array=0-${ARRAY_MAX}/" \
    scripts/sherlock/submit_all_tasks.sh

# Create log directory
mkdir -p logs

# Submit and capture job ID
SUBMIT_OUT=$(sbatch scripts/sherlock/submit_all_tasks.sh 2>&1)
echo "$SUBMIT_OUT" | grep -oP '\d+$'
REMOTE_SCRIPT
)

if [ -z "$JOB_ID" ]; then
    log ERROR "Failed to submit SLURM job"
    exit 1
fi
log SUBMIT "Submitted SLURM array job: ${JOB_ID} (${N_PAIRS} tasks)"

# ---------------------------------------------------------------------------
# Phase 4: Poll for completion
# ---------------------------------------------------------------------------
log WAIT "Polling job status every ${POLL_INTERVAL}s..."

while true; do
    # Count running/pending tasks for this job
    ACTIVE=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" \
        "squeue -j ${JOB_ID} -h 2>/dev/null | wc -l" || echo "0")
    ACTIVE=$(echo "$ACTIVE" | tr -d '[:space:]')

    if [ "$ACTIVE" -eq 0 ]; then
        log WAIT "All array tasks finished."
        break
    fi

    log WAIT "${ACTIVE} task(s) still running/pending..."
    sleep "$POLL_INTERVAL"
done

# Check for failures via sacct
log WAIT "Checking job exit statuses..."
FAILED=$(ssh $SSH_OPTS "${SHERLOCK_DEST}" \
    "sacct -j ${JOB_ID} --format=JobID,State,ExitCode --noheader -P 2>/dev/null | grep -c 'FAILED'" \
    || echo "0")
FAILED=$(echo "$FAILED" | tr -d '[:space:]')

if [ "$FAILED" -gt 0 ]; then
    log WARN "${FAILED} task(s) FAILED. Continuing with available results."
    ssh $SSH_OPTS "${SHERLOCK_DEST}" \
        "sacct -j ${JOB_ID} --format=JobID%-20,State,ExitCode,Elapsed --noheader 2>/dev/null | grep FAILED"
fi

# ---------------------------------------------------------------------------
# Phase 5: Transfer results back + collect
# ---------------------------------------------------------------------------
log COLLECT "Transferring results from Sherlock..."
rsync -avz -e "ssh $SSH_OPTS" \
    "${SHERLOCK_DEST}:${SHERLOCK_DIR}/results/fep_runs/" \
    "${PROJECT_DIR}/results/fep_runs/"

log COLLECT "Running result collection and analysis..."
cd "$PROJECT_DIR"
python -m src.main --collect-results

# Close SSH multiplexed connection
ssh $SSH_OPTS -O exit "${SHERLOCK_DEST}" 2>/dev/null || true

log DONE "Pipeline complete! Results in ${PROJECT_DIR}/results/"
log DONE "Key outputs:"
log DONE "  results/ddg_summary.csv"
log DONE "  results/correlation_analysis.csv"
log DONE "  results/plots/ddg_vs_fold_reduction.png"
