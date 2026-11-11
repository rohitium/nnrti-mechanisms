#!/usr/bin/env bash
# Stage 3 - non-equilibrium alchemical FEP (GPU cluster, ~3 weeks).
#
# 19 alchemical legs x {holo, apo} x 3 replicates x 100 forward + 100 reverse
# switches. Submitted from the cluster checkout:
#
#   ssh <cluster>
#   cd "$SCRATCH/nnrti-mechanisms-work"
#   export PROJECT_ROOT="$PWD"
#   # Prepare hybrids locally first; transfer inputs before cluster submission.
#   ./ops/slurm/fep/prepare_p0_neq.sh
#   ./ops/slurm/fep/submit_p0_neq_pipeline.sh
#
# Then, once the arrays drain, integrate the work values into the DDG panel:
#
#   python -m nnrti.fep.combine_neq
#
# Pipeline reference: ops/slurm/fep/README.md

set -euo pipefail
sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
