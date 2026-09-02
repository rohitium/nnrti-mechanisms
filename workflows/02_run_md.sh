#!/usr/bin/env bash
# Stage 2 - equilibrium MD (GPU cluster, ~2 weeks for 60 jobs).
#
# This is a pointer, not a runner: MD is submitted as a Slurm array from the
# cluster checkout. See scripts/sherlock/submit_md_batched.sh for the array and
# docs/ for the runbooks.
#
#   ssh <cluster>
#   cd nnrti-mechanisms && git pull
#   ./scripts/sherlock/submit_md_batched.sh

set -euo pipefail
cat "$(dirname "$0")/02_run_md.sh" | sed -n '2,12p' | sed 's/^# \{0,1\}//'
echo
echo "Per-task entry point: python -m nnrti.md.sherlock.run_md_job --help"
