#!/usr/bin/env bash
# Stage 2 - equilibrium MD (GPU cluster, ~2 weeks for 60 jobs).
#
# This is a pointer, not a runner: MD is submitted as a Slurm array from the
# cluster checkout. See ops/slurm/cluster/submit_md_batched.sh for the array and
# docs/ for the runbooks.
#
#   ssh <cluster>
#   cd "$SCRATCH/nnrti-mechanisms-work"
#   export PROJECT_ROOT="$PWD"
#   ./ops/slurm/cluster/submit_md_batched.sh

set -euo pipefail
sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
echo
echo "Per-task entry point: python -m nnrti.md.sherlock.run_md_job --help"
