#!/bin/bash
#
# Submit exact Perses point-mutation FEP thermodynamic-cycle jobs.
#
# Example pilot:
#   MUTATION_ALLOWLIST=V106A,Y181C,Y188L REPLICATES=1 \
#   FEP_N_CYCLES=500 FEP_N_STATES=11 ./scripts/sherlock/submit_perses_point_mutation_fep.sh
#
# Production example:
#   MUTATION_ALLOWLIST=V106A,V106I,Y181C,Y188L,G190A,G190E,Y318F \
#   REPLICATES=1,2,3 FEP_N_CYCLES=5000 ./scripts/sherlock/submit_perses_point_mutation_fep.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"

SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_GRES="${SHERLOCK_GRES:-gpu:1}"
SHERLOCK_TIME="${SHERLOCK_TIME:-48:00:00}"
SHERLOCK_MEM="${SHERLOCK_MEM:-32G}"
SHERLOCK_QOS="${SHERLOCK_QOS:-}"
CONDA_ENV="${CONDA_ENV:-nnrti-fep}"

MUTATION_ALLOWLIST="${MUTATION_ALLOWLIST:-V106A,V106I,Y181C,Y188L,G190A,G190E,Y318F}"
REPLICATES="${REPLICATES:-1,2,3}"
FEP_OUTPUT_DIR="${FEP_OUTPUT_DIR:-results/analysis/perses_point_mutation_fep}"
FEP_N_STATES="${FEP_N_STATES:-11}"
FEP_N_CYCLES="${FEP_N_CYCLES:-5000}"
FEP_STEPS_PER_CYCLE="${FEP_STEPS_PER_CYCLE:-250}"
FEP_PLATFORM="${FEP_PLATFORM:-CUDA}"
FEP_CHAIN_ID="${FEP_CHAIN_ID:-A}"

if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found. Run this script on Sherlock login node." >&2
    exit 1
fi

mkdir -p logs
IFS=',' read -r -a MUTATIONS <<< "$MUTATION_ALLOWLIST"
IFS=',' read -r -a REPS <<< "$REPLICATES"

for mutation in "${MUTATIONS[@]}"; do
    mutation="$(echo "$mutation" | xargs)"
    [ -n "$mutation" ] || continue
    for rep in "${REPS[@]}"; do
        rep="$(echo "$rep" | xargs)"
        [ -n "$rep" ] || continue
        wt_pdb="results/md_runs/wt/rep_$(printf '%02d' "$rep")/assets/wt_md_rep$(printf '%02d' "$rep")_start.pdb"
        job_name="pfep_${mutation}_${rep}"
        sbatch \
            --job-name="$job_name" \
            --partition="$SHERLOCK_PARTITION" \
            --gres="$SHERLOCK_GRES" \
            --time="$SHERLOCK_TIME" \
            --mem="$SHERLOCK_MEM" \
            ${SHERLOCK_QOS:+--qos="$SHERLOCK_QOS"} \
            --output="logs/${job_name}.%j.out" \
            --error="logs/${job_name}.%j.err" \
            --wrap="source ~/.bashrc && conda activate '$CONDA_ENV' && cd '$PROJECT_ROOT' && PYTHONPATH=. python -m src.analysis.cli.run_perses_point_mutation_cycle --mutation '$mutation' --chain-id '$FEP_CHAIN_ID' --wt-complex-pdb '$wt_pdb' --output-dir '$FEP_OUTPUT_DIR' --n-states '$FEP_N_STATES' --n-cycles '$FEP_N_CYCLES' --steps-per-cycle '$FEP_STEPS_PER_CYCLE' --platform '$FEP_PLATFORM' --skip-endstate-validation"
    done
done

