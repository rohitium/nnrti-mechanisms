#!/usr/bin/env bash
# Stage 3 - non-equilibrium alchemical FEP (GPU cluster, ~3 weeks).
#
# 19 alchemical legs x {holo, apo} x 3 replicates x 100 forward + 100 reverse
# switches. Submitted from the cluster checkout:
#
#   ssh <cluster>
#   cd nnrti-mechanisms && git pull
#   ./scripts/fep_pmx/prepare_p0_hybrids.sh
#   ./scripts/fep_pmx/prepare_p0_neq.sh
#   ./scripts/fep_pmx/submit_p0_neq_pipeline.sh
#
# Then, once the arrays drain, integrate the work values into the DDG panel:
#
#   python -m nnrti.fep.combine_neq
#
# Runbooks: scripts/fep_pmx/OPERATIONS.md, scripts/fep_pmx/RUNBOOK_G190E_SEM.md

set -euo pipefail
sed -n '2,20p' "$(dirname "$0")/03_run_fep.sh" | sed 's/^# \{0,1\}//'
