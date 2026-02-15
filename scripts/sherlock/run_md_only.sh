#!/bin/bash
#
# Submit ONLY MD jobs on Sherlock (script-based path).
# Thin wrapper around submit_md_batched.sh.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_SIZE="${1:-6}"
MAX_CONCURRENT="${2:-12}"

echo "Submitting MD jobs only via submit_md_batched.sh..."
echo "  batch_size=${BATCH_SIZE} max_concurrent=${MAX_CONCURRENT}"
echo

bash "${SCRIPT_DIR}/submit_md_batched.sh" "${BATCH_SIZE}" "${MAX_CONCURRENT}"
