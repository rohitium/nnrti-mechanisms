#!/usr/bin/env bash
set -euo pipefail

# Pull completed results from Sherlock and run the full analysis pipeline.
#
# Usage:
#   SHERLOCK_USER=rsatija bash scripts/rsync_and_analyze.sh
#
# Optional:
#   MD_PRODUCTION_NS=10.0   target ns for COMPLETE_ONLY filtering (default 10.0)

SHERLOCK_USER="${SHERLOCK_USER:-}"
if [[ -z "${SHERLOCK_USER}" ]]; then
  echo "Set SHERLOCK_USER first, e.g.: export SHERLOCK_USER=rsatija"
  exit 1
fi

export SHERLOCK_USER
export MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-10.0}"

echo "=== Pulling holo results ==="
COMPLETE_ONLY=1 bash scripts/rsync_results.sh pull

echo "=== Pulling apo results ==="
COMPLETE_ONLY=1 bash scripts/rsync_apo.sh pull

echo "=== Running holo analysis ==="
bash scripts/run_analysis.sh

echo "=== Running apo analysis ==="
bash scripts/run_apo_analysis.sh

echo "Done."
