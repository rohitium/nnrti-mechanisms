#!/usr/bin/env bash
# Prepare missing Perses holo/apo hybrids in parallel (local Mac).
#
# Usage:
#   MAX_JOBS=6 bash scripts/fep_jorgensen/prepare_panel_parallel.sh
#   MAX_JOBS=4 bash scripts/fep_jorgensen/prepare_panel_parallel.sh --holo-only
#   MAX_JOBS=6 bash scripts/fep_jorgensen/prepare_panel_parallel.sh --apo-only
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

OUTPUT_DIR="${FEP_OUTPUT_DIR:-results/analysis/fep_jorgensen}"
MAX_JOBS="${MAX_JOBS:-6}"
PYTHON="${FEP_PYTHON:-/Users/rohitpro/miniconda3/envs/nnrti-prep/bin/python}"
LOG_DIR="${FEP_PREP_LOG_DIR:-logs/fep_jorgensen_prepare_parallel}"
MODE="all"
for arg in "$@"; do
    case "$arg" in
        --holo-only) MODE="holo" ;;
        --apo-only) MODE="apo" ;;
        *) echo "Unknown arg: $arg" >&2; exit 1 ;;
    esac
done

mkdir -p "$LOG_DIR"

run_one() {
    local phase="$1"
    local mutation="$2"
    local start_label="$3"
    local end_label="$4"
    local input_pdb="$5"
    local leg_id
    leg_id="$(PYTHONNOUSERSITE=1 PYTHONPATH=. "$PYTHON" - <<PY
from scripts.fep_jorgensen.mutations import MutationLeg
print(MutationLeg("$start_label", "$end_label", "$mutation").leg_id)
PY
)"
    local marker="$OUTPUT_DIR/legs/${leg_id}/${phase}/hybrid_system.xml"
    if [[ -f "$marker" ]]; then
        echo "SKIP ${leg_id} ${phase} (exists)"
        return 0
    fi
    local log="$LOG_DIR/${leg_id}_${phase}.log"
    echo "START ${leg_id} ${phase} -> $log"
    PYTHONNOUSERSITE=1 PYTHONPATH=. "$PYTHON" -m scripts.fep_jorgensen.prepare \
        --backend perses \
        --phase "$phase" \
        --mutation "$mutation" \
        --start-label "$start_label" \
        --end-label "$end_label" \
        --input-complex-pdb "$input_pdb" \
        --output-dir "$OUTPUT_DIR" \
        >"$log" 2>&1
    echo "DONE ${leg_id} ${phase}"
}

export -f run_one
export OUTPUT_DIR PYTHON PROJECT_ROOT

TASKS="$LOG_DIR/tasks.tsv"
: >"$TASKS"

append_tasks() {
    local phase="$1"
    while IFS=$'\t' read -r mutation start end pdb; do
        [[ -n "$mutation" ]] || continue
        echo -e "${phase}\t${mutation}\t${start}\t${end}\t${pdb}" >>"$TASKS"
    done
}

if [[ "$MODE" == "all" || "$MODE" == "holo" ]]; then
    append_tasks holo <<'EOF'
L100I	K103N	L100I+K103N	results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
F227L	V106A	V106A+F227L	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
L234I	V106A	V106A+L234I	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
P225H	V106A	V106A+P225H	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
F227C	V106I	V106I+F227C	results/md_runs/V106I/rep_01/assets/V106I_md_rep01_start.pdb
EOF
fi

if [[ "$MODE" == "all" || "$MODE" == "apo" ]]; then
    append_tasks apo <<'EOF'
F227C	WT	F227C	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
G190A	WT	G190A	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
G190E	WT	G190E	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
G190S	WT	G190S	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
K103N	WT	K103N	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
V106I	WT	V106I	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
V106M	WT	V106M	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
Y181C	WT	Y181C	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
Y188L	WT	Y188L	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
Y318F	WT	Y318F	results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
A98G	F227C	A98G+F227C	results/md_runs/F227C/rep_01/assets/F227C_md_rep01_start.pdb
M230L	K103N	K103N+M230L	results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
P225H	K103N	K103N+P225H	results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
L100I	K103N	L100I+K103N	results/md_runs/K103N/rep_01/assets/K103N_md_rep01_start.pdb
F227L	V106A	V106A+F227L	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
L234I	V106A	V106A+L234I	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
P225H	V106A	V106A+P225H	results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb
F227C	V106I	V106I+F227C	results/md_runs/V106I/rep_01/assets/V106I_md_rep01_start.pdb
EOF
fi

echo "Tasks: $(wc -l <"$TASKS")  MAX_JOBS=$MAX_JOBS  MODE=$MODE"
echo "Logs: $LOG_DIR"

if command -v parallel >/dev/null 2>&1; then
    parallel --jobs "$MAX_JOBS" --colsep '\t' --halt soon,fail=1 run_one {1} {2} {3} {4} {5} :::: "$TASKS"
else
    while IFS=$'\t' read -r phase mutation start end pdb; do
        while (( $(jobs -pr | wc -l | tr -d ' ') >= MAX_JOBS )); do
            sleep 2
        done
        run_one "$phase" "$mutation" "$start" "$end" "$pdb" &
    done <"$TASKS"
    wait
fi

echo "Parallel panel prep finished."
