#!/bin/bash
#
# Resume the G190E 20 ns-equilibration campaign: audit equil, then chain
# extract -> switch.
#
# WHY THIS EXISTS
#   On 2026-08-28 equil array element 3 died on a transient CUDA fault
#   (cudaErrorLaunchFailure after 9200 steps). Because the stages are chained
#   `afterok`, that made the dependency permanently unsatisfiable and the queued
#   extract/switch jobs (41149158, 41149159) had to be cancelled. They must be
#   resubmitted once all 12 equil units are present -- 11 from array 41149155
#   plus the repair job 41155898 (holo rep_02, lambda 1).
#
#   This script does the audit, refuses to proceed unless equil is complete, and
#   submits the two remaining stages with per-batch TASK_ID_FILEs.
#
# USAGE
#   bash ops/slurm/fep/resume_g190e_chain.sh            # audit + submit
#   AUDIT_ONLY=1 bash ops/slurm/fep/resume_g190e_chain.sh
#
# See RUNBOOK_G190E_SEM.md for the full campaign and the read-outs that matter.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"
export PYTHONPATH="${PROJECT_ROOT:-$PWD}/src:${PYTHONPATH:-}"

MANIFEST="${MANIFEST:-results/analysis/fep_pmx/neq_g190e_equil20_manifest.csv}"
EXCLUDE_NODES="${EXCLUDE_NODES:-sh03-12n12,sh02-16n06,sh03-13n01,sh03-12n01}"
export MANIFEST EXCLUDE_NODES

echo "=========================================="
echo "G190E chain resume"
echo "=========================================="
echo "Manifest: $MANIFEST"
echo ""

if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: manifest not found. Did you git pull?" >&2; exit 1
fi

# --- 1. is anything still running? -------------------------------------------
ACTIVE="$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -c "pmx_neq" || true)"
if [ "${ACTIVE:-0}" -gt 0 ]; then
    echo "NOTE: $ACTIVE pmx_neq job(s) still active:"
    squeue -u "$USER" -o "%.14i %.16j %.8T %.10M %R" | grep -E "pmx_neq|JOBID" || true
    echo ""
fi

# --- 2. equil must be 12/12 --------------------------------------------------
echo "--- audit ---"
# Capture the FULL audit, not just the '===' lines. Filtering here with a
# pipeline meant that when the audit printed nothing matching (or failed
# outright), grep returned non-zero, `set -e` killed the script, and the real
# error was never shown. Never let a filter decide whether an error is visible.
# NOTE: audit_neq_panel.py exits NON-ZERO whenever any unit is incomplete, which
# is the normal mid-campaign state -- not an error. Do not let `set -e` or
# pipefail treat that as a crash; judge completeness from the summary lines.
AUDIT_RAW="$(python -m nnrti.fep.audit_neq_panel --manifest "$MANIFEST" 2>&1 || true)"
AUDIT="$(echo "$AUDIT_RAW" | grep '===' || true)"
if [ -z "$AUDIT" ]; then
    echo "ERROR: audit produced no '===' summary lines. Full output:" >&2
    echo "$AUDIT_RAW" >&2
    exit 1
fi
echo "$AUDIT"
echo ""

EQUIL_LINE="$(echo "$AUDIT" | grep -i 'EQUIL' || true)"
if ! echo "$EQUIL_LINE" | grep -q "12/12 ok"; then
    echo "STOP: equil is not complete -- $EQUIL_LINE" >&2
    echo "" >&2
    echo "If jobs are still RUNNING, wait. If an element FAILED, find it with:" >&2
    echo "  sacct -j <arrayjobid> --format=JobID%20,State,ExitCode,Elapsed,NodeList%20 | grep -v COMPLETED" >&2
    echo "then repair it exactly as on 2026-08-28: build a one-row manifest for the" >&2
    echo "failed panel_task_id, submit STAGE=equil against it with its own" >&2
    echo "TASK_ID_FILE and EXCLUDE_NODES. Recipe in RUNBOOK_G190E_SEM.md." >&2
    exit 1
fi
echo "equil complete."

# Guard: a truncated equil.gro counts as "present" to run_neq_task and would be
# skipped. All 12 must share an atom count.
echo ""
echo "--- equil.gro sanity (line counts must match within each replicate) ---"
for ph in holo apo; do for r in 01 02 03; do for l in 0 1; do
  d="results/analysis/fep_pmx/legs/wt_to_G190E/$ph/rep_$r/neq/eq_lambda$l"
  printf "  %-4s rep%s l%s: " "$ph" "$r" "$l"
  if [ -f "$d/equil.gro" ]; then wc -l < "$d/equil.gro"; else echo "MISSING"; fi
done; done; done

if [ "${AUDIT_ONLY:-0}" = "1" ]; then
    echo ""; echo "AUDIT_ONLY set -- stopping before submission."; exit 0
fi

# --- 3. chain extract -> switch ----------------------------------------------
echo ""
echo "--- submitting extract -> switch ---"
EXTRACT=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_extract_g190e_e20_task_ids.txt \
          STAGE=extract bash ops/slurm/fep/submit_p0_neq.sh | tail -1)
echo "extract=$EXTRACT"
SWITCH=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_switch_g190e_e20_task_ids.txt \
         STAGE=switch DEPENDENCY=afterok:$EXTRACT bash ops/slurm/fep/submit_p0_neq.sh | tail -1)
echo "switch=$SWITCH"

echo ""
for j in $EXTRACT $SWITCH; do
    echo "== $j =="; scontrol show job "$j" | grep -ioE 'Dependency=[^ ]*' || true
done

echo ""
echo "=========================================="
echo "Submitted. Expect 12 extract (CPU, minutes) then 60 switch (~8 h each)."
echo ""
echo "Sanity checks that matter:"
echo "  * switch must be 60 tasks in ONE chunk. 120 would exceed the gpu QOS"
echo "    cap of 100 submitted jobs and be rejected."
echo "  * each submit must print a chunk file with the g190e_e20 qualifier."
echo ""
echo "Monitor:  python -m nnrti.fep.audit_neq_panel --manifest \$MANIFEST | grep '==='"
echo ""
echo "When SWITCH is 60/60, analyse ON SHERLOCK and WITHOUT --force:"
echo "  python -m nnrti.fep.combine_neq --targets G190E --replicates 3"
echo "  python -m nnrti.fep.qc_neq --legs wt_to_G190E --replicates 3"
echo ""
echo "Then REBUILD THE FULL PANEL -- combine_neq --targets X rewrites panel_ddg.csv"
echo "with only those targets. Full genotype list is in RUNBOOK_G190E_SEM.md."
echo ""
echo "Read-out: sigma_DDG (spread of the three per-rep ddG), NOT the SEM and NOT"
echo "the per-phase spread. 3.06 -> under ~1.7 means n=3 clears SEM<1 alone."
echo "Secondary: sd(holo-apo hysteresis), currently 7.39; toward ~2 means the"
echo "mechanism worked even though overlap stays ~0."
echo "=========================================="
