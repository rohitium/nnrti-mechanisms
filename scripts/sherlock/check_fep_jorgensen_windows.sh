#!/bin/bash
#
# Summarize FEP lambda-window batch job status on Sherlock.
#
# Usage:
#   export PROJECT_ROOT=$PWD
#   ./scripts/sherlock/check_fep_jorgensen_windows.sh           # latest fep_jorgensen array
#   ./scripts/sherlock/check_fep_jorgensen_windows.sh 35519693  # explicit job id
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
WINDOWS="${WINDOWS:-$PROJECT_ROOT/results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows}"
TARGET_SAMPLES=1000

job_id="${1:-}"
if [[ -z "$job_id" ]]; then
    job_id="$(ls -1 "$PROJECT_ROOT/logs"/fep_jorgensen.*_0.out 2>/dev/null | sed -E 's|.*/fep_jorgensen\.([0-9]+)_0\.out|\1|' | sort -n | tail -1 || true)"
fi
if [[ -z "$job_id" ]]; then
    echo "No fep_jorgensen batch logs found under $PROJECT_ROOT/logs" >&2
    exit 1
fi

echo "=== Job $job_id ==="
if command -v sacct >/dev/null 2>&1; then
    sacct -j "$job_id" --format=JobID,State,ExitCode,Elapsed,Timelimit,MaxRSS,ReqMem,NodeList -P \
        | column -t -s'|'
else
    echo "(sacct not available on this host)"
fi
echo

echo "=== Log errors (non-warning) ==="
found_err=0
for err in "$PROJECT_ROOT/logs/fep_jorgensen.${job_id}"_*.err; do
    [[ -f "$err" ]] || continue
    if grep -vE 'SyntaxWarning|Warning on use of the timeseries|PyMBAR will use' "$err" | grep -qiE 'error|traceback|exception|killed|oom|cuda|failed'; then
        echo "--- $(basename "$err") ---"
        grep -vE 'SyntaxWarning|Warning on use of the timeseries|PyMBAR will use' "$err" \
            | grep -iE 'error|traceback|exception|killed|oom|cuda|failed' || true
        found_err=1
    fi
done
if [[ "$found_err" -eq 0 ]]; then
    echo "none (only OpenMM SyntaxWarnings — ignore)"
fi
echo

echo "=== Window outputs ($WINDOWS) ==="
if [[ ! -d "$WINDOWS" ]]; then
    echo "Missing windows dir: $WINDOWS" >&2
    exit 1
fi

complete=0
partial=0
missing=0
for state in $(seq 0 10); do
    tag="$(printf '%02d' "$state")"
    csv="$WINDOWS/state_${tag}_energies.csv"
    chk="$WINDOWS/state_${tag}.chk"
    if [[ ! -f "$csv" ]]; then
        echo "state ${tag}: MISSING csv"
        missing=$((missing + 1))
        continue
    fi
    samples=$(( $(wc -l < "$csv") - 1 ))
    chk_note=""
    [[ -f "$chk" ]] && chk_note=" chk=yes"
    if [[ "$samples" -eq "$TARGET_SAMPLES" ]]; then
        echo "state ${tag}: OK  ${samples}/${TARGET_SAMPLES} samples${chk_note}"
        complete=$((complete + 1))
    else
        echo "state ${tag}: PARTIAL ${samples}/${TARGET_SAMPLES} samples${chk_note}"
        partial=$((partial + 1))
    fi
done
echo
echo "Summary: ${complete}/11 complete, ${partial} partial, ${missing} missing"
if [[ "$complete" -eq 11 ]]; then
    echo "Ready for rsync + local: python -m scripts.fep_jorgensen.analyze --target V106A"
elif grep -q OUT_OF_MEMORY <(sacct -j "$job_id" --format=State -P 2>/dev/null || true); then
    echo "Hint: resubmit failed tasks with SHERLOCK_MEM=32G (or 64G)"
fi
