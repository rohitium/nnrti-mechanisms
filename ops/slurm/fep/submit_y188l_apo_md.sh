#!/bin/bash
#
# Submit Y188L apo MD (P0 blocker) on Sherlock.
#
# Usage:
#   MUTATION_ALLOWLIST=y188l ./ops/slurm/fep/submit_y188l_apo_md.sh
#
# Or locally (GPU/OpenMM):
#   conda activate nnrti-openmm   # or env with OpenMM CUDA/CPU
#   python -m nnrti.fep.run_apo_md --mutations Y188L --production-ns 100

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

export MUTATION_ALLOWLIST="${MUTATION_ALLOWLIST:-y188l}"
export MD_PRODUCTION_NS="${MD_PRODUCTION_NS:-100.0}"
export MD_FORCE_RERUN="${MD_FORCE_RERUN:-0}"
export SKIP_IF_AT_TARGET="${SKIP_IF_AT_TARGET:-1}"

exec bash "$PROJECT_ROOT/ops/slurm/cluster/submit_apo_md_batched.sh" "${1:-3}" "${2:-6}"
