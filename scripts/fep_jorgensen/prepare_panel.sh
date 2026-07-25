#!/bin/bash
#
# Prepare all manuscript FEP legs with Perses hybrid topology (holo + apo).
# Skips legs that already have both holo/ and apo/ hybrid_system.xml.
#
# Usage (local Mac, nnrti-prep after setup_perses_env.sh):
#   bash scripts/fep_jorgensen/prepare_panel.sh
#   bash scripts/fep_jorgensen/prepare_panel.sh --force   # redo all legs
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

OUTPUT_DIR="${FEP_OUTPUT_DIR:-results/analysis/fep_jorgensen}"
REPLICATE="${FEP_REPLICATE:-1}"
FORCE=0
if [[ "${1:-}" == "--force" ]]; then
    FORCE=1
fi

PYTHON="${FEP_PYTHON:-python}"
if ! command -v "$PYTHON" >/dev/null 2>&1; then
    PYTHON=python3
fi

run_prepare() {
    local mutation="$1"
    local start_label="$2"
    local end_label="$3"
    local input_pdb="$4"
    local leg_id
    leg_id="$(PYTHONPATH=. "$PYTHON" - <<PY
from scripts.fep_jorgensen.mutations import MutationLeg
print(MutationLeg("$start_label", "$end_label", "$mutation").leg_id)
PY
)"
    local holo_xml="$OUTPUT_DIR/legs/${leg_id}/holo/hybrid_system.xml"
    local apo_xml="$OUTPUT_DIR/legs/${leg_id}/apo/hybrid_system.xml"
    if [[ "$FORCE" -eq 0 && -f "$holo_xml" && -f "$apo_xml" ]]; then
        echo "SKIP ${leg_id} (exists: $holo_xml and $apo_xml)"
        return 0
    fi
    echo "PREP ${leg_id} (${start_label} -> ${end_label}, mutation ${mutation}, holo+apo)"
    PYTHONNOUSERSITE=1 PYTHONPATH=. "$PYTHON" -m scripts.fep_jorgensen.prepare \
        --backend perses \
        --phase all \
        --mutation "$mutation" \
        --start-label "$start_label" \
        --end-label "$end_label" \
        --input-complex-pdb "$input_pdb" \
        --output-dir "$OUTPUT_DIR"
}

# WT -> single mutants (run singles before compound legs for easier monitoring)
run_prepare F227C WT F227C results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare G190A WT G190A results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare G190E WT G190E results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare G190S WT G190S results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare K103N WT K103N results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare V106A WT V106A results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare V106I WT V106I results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare V106M WT V106M results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare Y181C WT Y181C results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare Y188L WT Y188L results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
run_prepare Y318F WT Y318F results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb

# compound legs (background MD structures already on disk)
run_prepare A98G F227C A98G+F227C results/md_runs/F227C/rep_01/assets/F227C_md_rep01_start.pdb
run_prepare M230L K103N K103N+M230L results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
run_prepare P225H K103N K103N+P225H results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
run_prepare L100I K103N L100I+K103N results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
run_prepare F227L V106A V106A+F227L results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
run_prepare L234I V106A V106A+L234I results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
run_prepare P225H V106A V106A+P225H results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
run_prepare F227C V106I V106I+F227C results/md_runs/V106I/rep_01/assets/V106I_md_rep01_start.pdb

echo "Panel prep complete."
