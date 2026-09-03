#!/usr/bin/env bash
# Stage 4 - trajectory analysis (laptop).
#
# Consumes the stripped analysis DCDs under results/md_runs/ and produces the
# derived tables that stage 5 turns into manuscript artifacts. Everything here
# is deterministic given the trajectories; nothing needs a GPU.
#
#   ./workflows/04_analysis.sh
#
# MDAnalysis and MDTraj are required, so this needs the nnrti-prep environment.

set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:-src}"

run() { echo; echo "=== $* ==="; "$@"; }

# Convergence: DOR pose RMSD and DOR-RT centre-of-mass distance (Supp. Fig. 1).
run "$PYTHON" -m nnrti.cli.compute_md_convergence

# Interface geometry: the observables behind Table 3 and Figure 3.
run "$PYTHON" -m nnrti.cli.compute_mechanism_coordinates
run "$PYTHON" -m nnrti.cli.compute_dor_moiety_contacts \
    --mutations WT F227C G190A G190E G190S K103N K103N+M230L K103N+P225H \
                L100I+K103N V106A V106A+F227L V106A+L234I V106A+P225H \
                V106I V106I+F227C V106M Y181C Y188L Y318F A98G+F227C

# NNIBP pocket volume (the V(NNIBP) column of Table 3).
run "$PYTHON" -m nnrti.cli.compute_modern_md_suite

# MM/GBSA interface energies. 100 evenly spaced snapshots, no frame filtering,
# minimisation to convergence, double precision -- see Methods.
run "$PYTHON" -m nnrti.cli.compute_mmgbsa_safe --frame-sampling even --snapshots 100
run "$PYTHON" -m nnrti.cli.rebuild_binding_energy_sources

echo; echo "Stage 4 complete. Run workflows/05_manuscript_artifacts.sh next."
