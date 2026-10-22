#!/bin/bash
#
# Densely sampled apo collapse runs WITH THE CAVITY FILLED, as a Slurm array.
#
# The companion array (submit_collapse.sh) runs the systems as prepared, with
# the NNIBP left empty, and is the control. These runs fill the cavity with
# water first, settle it with the protein restrained so the pocket is held open
# while excess water leaves, then release and sample densely. Comparing the two
# separates genuine closing from relaxation into a vacuum.
#
# Pocket volumes are NOT computed here. Coordinates come back and are measured
# locally with the same routine as every other number in the paper.
#
# Usage:
#   bash ops/slurm/apo_collapse/submit_collapse.sh [PICOSECONDS]
#
set -euo pipefail

PS_TOTAL="${1:-200}"
SETTLE_PS="${2:-100}"
GH="${GROUP_HOME:?GROUP_HOME not set}"
PROJECT="${GH}/nnrti-mechanisms"
ASSETS="${PROJECT}/md_runs/apo"
OUT="${SCRATCH:?SCRATCH not set}/apo_collapse_solvated"
LOGS="${OUT}/logs"
MANIFEST="${OUT}/manifest.txt"

mkdir -p "${OUT}" "${LOGS}"

# one line per genotype/replicate that actually has its start assets
: > "${MANIFEST}"
for d in "${ASSETS}"/*/; do
    g=$(basename "$d")
    for rep in 1 2 3; do
        r=$(printf "%02d" "${rep}")
        if [[ -f "${d}/rep_${r}/assets/${g}_apo_md_rep${r}_start.pdb" ]]; then
            echo "${g} ${rep}" >> "${MANIFEST}"
        fi
    done
done
N=$(wc -l < "${MANIFEST}")
echo "manifest: ${N} runs -> ${MANIFEST}"
echo "production ${PS_TOTAL} ps after ${SETTLE_PS} ps restrained settling"
[[ "${N}" -gt 0 ]] || { echo "no runs found under ${ASSETS}" >&2; exit 1; }

sbatch --parsable \
    --job-name=apo_solv \
    --partition=gpu \
    --gres=gpu:1 \
    --time=02:00:00 \
    --mem=16G \
    --cpus-per-task=4 \
    --array=1-"${N}"%12 \
    --output="${LOGS}/collapse_%A_%a.log" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail
module load chemistry py-openmm/8.1.1_py312

read -r G REP <<< "\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
echo "task \${SLURM_ARRAY_TASK_ID}: \${G} rep\${REP} on \$(hostname)"

python3 "${PROJECT}/ops/apo_collapse/collapse_run_solvated.py" \
    "\${G}" "\${REP}" "${ASSETS}" "${OUT}" "${PS_TOTAL}" "${SETTLE_PS}"
SBATCH_EOF
