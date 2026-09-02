#!/bin/bash
#
# Extend WT apo MD from 10 ns -> 100 ns on Sherlock.
#
# Context: the WT apo runs (results/md_runs/apo/wt/rep_{01,02,03}) were only
# taken to 10 ns (5,000,000 steps) because MM/GBSA did not need more. WT apo is
# the shared endpoint for every single-mutation FEP leg, so to seed the FEP
# endpoints from a converged 100 ns ensemble (instead of the 5 ns hybrid
# re-equilibration) we first extend it to 100 ns.
#
# The 10 ns runs are marked status="ok", so a plain rerun would skip them:
# we force the rerun but keep SKIP_IF_AT_TARGET so anything already at 100 ns is
# left alone, and resume from the existing checkpoint (runs the *remaining*
# 45M steps, i.e. to 100 ns total -- see src/nnrti/md/openmm/md_protocol.py).
#
# Usage (on Sherlock, after Batch A finishes):
#   ./scripts/fep_pmx/submit_wt_apo_md.sh            # 3 reps, up to 6 concurrent
#   ./scripts/fep_pmx/submit_wt_apo_md.sh 3 6        # explicit batch/concurrency
#
# Each SLURM job is 12 h; 90 ns of extension will not finish in one job, so
# rerun this same command after each batch completes -- it resumes each rep's
# checkpoint and skips reps that have reached 100 ns.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

export MUTATION_ALLOWLIST="${MUTATION_ALLOWLIST:-wt}"
export MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-100.0}"
export MD_FORCE_RERUN="${MD_FORCE_RERUN:-1}"          # override status=ok on the 10 ns runs
export SKIP_IF_AT_TARGET="${SKIP_IF_AT_TARGET:-1}"     # leave reps already at 100 ns alone
export MD_RESUME_FROM_CHECKPOINT="${MD_RESUME_FROM_CHECKPOINT:-1}"

exec bash "$PROJECT_ROOT/scripts/sherlock/submit_apo_md_batched.sh" "${1:-3}" "${2:-6}"
