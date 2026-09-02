#!/usr/bin/env bash
# Stage 1 - build the 19 RT-DOR systems (laptop, minutes).
#
# Starts from PDB 4NCG, introduces each p66 substitution in silico, rebuilds DOR
# from its SDF template and writes manifests/md_manifest.csv, which stages 2-4
# all read.

set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-$HOME/miniconda3/envs/nnrti-prep/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

"$PYTHON" -m nnrti.structure_prep.preparation "$@"
echo "Wrote manifests/md_manifest.csv"
