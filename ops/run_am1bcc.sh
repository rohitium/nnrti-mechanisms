#!/bin/bash
# Wait for the unrestrained relaxation run to free the CPU, then re-score the
# charge-model subset with AM1-BCC under the main-text protocol.
set -u
cd /Users/rohitpro/Career/00_Github/nnrti-mechanisms
while pgrep -f "compute_mmgbsa_safe.*relaxsens" >/dev/null; do sleep 120; done
echo "=== relaxation run finished, starting AM1-BCC $(date +%H:%M:%S) ==="
PYTHONPATH=src ~/miniconda3/envs/nnrti-prep/bin/python -m nnrti.cli.compute_mmgbsa_safe \
    --manifest manifests/md_manifest_charge_subset.csv \
    --output results/.checkpoints/.checkpoint_mmgbsa_am1bcc.csv \
    --snapshots 100 --discard-fraction 0.25 --frame-sampling even \
    --snapshot-relaxation unrestrained --relaxation-iterations 2000 \
    --charge-method am1bcc --workers 4 --force
echo "=== AM1-BCC exit=$? $(date +%H:%M:%S) ==="
