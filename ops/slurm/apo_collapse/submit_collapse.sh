#!/bin/bash
#
# Densely sampled apo collapse runs for the whole panel, as a Slurm array.
#
# The production apo trajectories save a frame every 50-500 ps and the pocket
# is already shut at the first one, so the collapse cannot be drawn from them.
# These runs save a frame every picosecond. A 5 ps test showed the NNIBP going
# from 170 A^3 to ~80 A^3 within 2 ps, so 200 ps is mostly baseline.
#
# Pocket volumes are NOT computed here. Coordinates come back and are measured
# locally with the same routine as every other number in the paper.
#
# Usage:
#   bash ops/slurm/apo_collapse/submit_collapse.sh [PICOSECONDS]
#
set -euo pipefail

PS_TOTAL="${1:-200}"
GH="${GROUP_HOME:?GROUP_HOME not set}"
PROJECT="${GH}/nnrti-mechanisms"
ASSETS="${PROJECT}/md_runs/apo"
OUT="${SCRATCH:?SCRATCH not set}/apo_collapse"
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
[[ "${N}" -gt 0 ]] || { echo "no runs found under ${ASSETS}" >&2; exit 1; }

sbatch --parsable \
    --job-name=apo_collapse \
    --partition=gpu \
    --gres=gpu:1 \
    --time=01:00:00 \
    --mem=16G \
    --cpus-per-task=4 \
    --array=1-"${N}"%12 \
    --output="${LOGS}/collapse_%A_%a.log" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail
module load chemistry py-openmm/8.1.1_py312

read -r G REP <<< "\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
echo "task \${SLURM_ARRAY_TASK_ID}: \${G} rep\${REP} on \$(hostname)"

python3 "${PROJECT}/ops/apo_collapse/collapse_run.py" \
    "\${G}" "\${REP}" "${ASSETS}" "${OUT}" "${PS_TOTAL}"
SBATCH_EOF
