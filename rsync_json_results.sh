#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   SHERLOCK_USER=rsatija bash rsync_json_results.sh

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija"
  exit 1
fi

REMOTE="${SHERLOCK_USER}@login.sherlock.stanford.edu:/scratch/users/${SHERLOCK_USER}/nnrti-mechanisms/results/fep_runs/"
LOCAL="results/fep_runs/"

mkdir -p "${LOCAL}"

rsync -avz \
  --include='*/' \
  --include='*.json' \
  --exclude='*' \
  "${REMOTE}" \
  "${LOCAL}"

echo "JSON sync complete. You can now run: python -m src.main --collect-results"
