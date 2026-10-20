#!/bin/bash
# Supplementary Table 5: MM/GBSA sensitivity to the snapshot relaxation protocol.
# Same subset, same sampling as the published run; only --snapshot-relaxation varies.
set -u
FAIL=0
cd /Users/rohitpro/Career/00_Github/nnrti-mechanisms
PY=~/miniconda3/envs/nnrti-prep/bin/python
MAN=manifests/md_manifest_relaxation_subset.csv
for MODE in none h_relax unrestrained; do
  OUT="results/.checkpoints/.checkpoint_mmgbsa_relaxsens_${MODE}.csv"
  echo "=== $MODE -> $OUT  ($(date +%H:%M:%S)) ==="
  PYTHONPATH=src $PY -m nnrti.cli.compute_mmgbsa_safe \
      --manifest "$MAN" --output "$OUT" \
      --snapshots 100 --discard-fraction 0.25 --frame-sampling even \
      --snapshot-relaxation "$MODE" --relaxation-iterations 2000 \
      --workers 4 --force
  RC=$?; [ $RC -ne 0 ] && FAIL=1; echo "=== $MODE exit=$RC ($(date +%H:%M:%S)) ==="
done
[ $FAIL -eq 0 ] && echo "ALL DONE OK" || echo "ALL DONE WITH FAILURES"
