#!/usr/bin/env bash
# Stage 5 - regenerate every numbered table and figure in the manuscript.
#
#   ./workflows/05_manuscript_artifacts.sh
#
# Reads only committed/deposited derived data, so it runs in about five minutes
# on a laptop. manuscript/ARTIFACTS.md maps each output back to its figure or
# table number.

set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-$HOME/miniconda3/envs/nnrti-prep/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

run() { echo; echo "=== $* ==="; "$@"; }

# --- Tables -----------------------------------------------------------------
run "$PYTHON" -m nnrti.cli.plot_dor_susceptibility_bars      # Table 1
run "$PYTHON" -m nnrti.cli.build_table_2                     # Table 2
run "$PYTHON" -m nnrti.cli.build_supplementary_table_4       # Table 3 + Supp. Table 4
run "$PYTHON" -m nnrti.cli.build_supplementary_table_3       # Supp. Table 3

# --- Figures ----------------------------------------------------------------
run "$PYTHON" -m nnrti.cli.plot_dor_schematic                # Figure 1B
run "$PYTHON" -m nnrti.fep.plot_protocol_schematic           # Figure 2A-C
run "$PYTHON" -m nnrti.cli.plot_panel_by_resistance_category # Figure 2D
run "$PYTHON" -m nnrti.cli.plot_mechanism_panel              # Figure 3
run "$PYTHON" -m nnrti.cli.plot_convergence_panel            # Supp. Figure 1
run "$PYTHON" -m nnrti.cli.plot_fep_work_distributions       # Supp. Figure 2

# --- Statistics quoted in the text ------------------------------------------
run "$PYTHON" -m nnrti.cli.classification_performance

echo; echo "Done. See manuscript/ARTIFACTS.md for the artifact -> file map."
