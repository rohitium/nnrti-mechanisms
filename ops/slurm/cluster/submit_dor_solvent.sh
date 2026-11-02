#!/bin/bash
#
# Triplicate 100 ns MD of doravirine alone in explicit solvent.
#
# The free-ligand leg of the binding thermodynamic cycle: the reference against
# which the bound torsions measure conformational restriction, and the solution
# term in dE_strain. See src/nnrti/cli/run_dor_solvent.py for why every ligand
# parameter is pinned to the production values.
#
# The system is ~2,500 atoms against ~201,000 for the solvated complex, so each
# replicate gets its own GPU rather than sharing one: without MPS, CUDA
# time-slices contexts and three processes on one device buy little over one.
#
# Usage:
#   bash ops/slurm/cluster/submit_dor_solvent.sh
#   REPLICATES="2 3" SHERLOCK_TIME=08:00:00 bash ops/slurm/cluster/submit_dor_solvent.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
cd "$PROJECT_ROOT"

RUNS_ROOT="${DOR_SOLVENT_ROOT:-/scratch/users/${USER}/nnrti-mechanisms-git/results/md_runs/dor_solvent}"
# Optional. The ligand atom names the torsion analysis needs are pinned in
# run_dor_solvent.PRODUCTION_ATOM_NAMES; set this to a holo
# *_analysis_topology.pdb to read them from a run instead and cross-check.
# Not defaulted to one: on the cluster those paths are symlinks into a
# visualization directory that the purge has already broken.
NAME_REFERENCE="${NAME_REFERENCE:-}"
LOG_DIR="${PROJECT_ROOT}/logs/dor_solvent"
mkdir -p "${LOG_DIR}"

REPLICATES="${REPLICATES:-1 2 3}"
PRODUCTION_NS="${PRODUCTION_NS:-100.0}"
ANALYSIS_FRAMES="${ANALYSIS_FRAMES:-2000}"
# Generous: throughput for a system this small on the gpu partition's mixed
# hardware has not been measured yet. Checkpoints land every 5,000 steps
# (0.01 ns), so a TIMEOUT costs essentially nothing and the job resumes.
SHERLOCK_TIME="${SHERLOCK_TIME:-16:00:00}"
SHERLOCK_PARTITION="${SHERLOCK_PARTITION:-gpu}"
SHERLOCK_MEM="${SHERLOCK_MEM:-16G}"

NAME_REFERENCE_ARG=()
if [ -n "${NAME_REFERENCE}" ]; then
    if [ ! -f "${NAME_REFERENCE}" ]; then
        echo "NAME_REFERENCE set but not found: ${NAME_REFERENCE}" >&2
        exit 1
    fi
    NAME_REFERENCE_ARG=(--name-reference "${NAME_REFERENCE}")
fi

# Prep once, here, rather than in the jobs: three jobs starting together would
# all see the assets missing and race to write the same two files. It is a few
# seconds of CPU on ~2,500 atoms, which is light enough for a login node.
if [ ! -f "${RUNS_ROOT}/prepared/dor_solvent_system.xml" ]; then
    echo "[prep] building solvated DOR assets in ${RUNS_ROOT}/prepared"
    module reset >/dev/null 2>&1 || true
    module load chemistry py-openmm/8.1.1_py312
    PYTHONPATH="${GROUP_HOME}/nnrti-mechanisms-repo/src:${PYTHONPATH:-}" \
        python3 -m nnrti.cli.run_dor_solvent \
            --root "${RUNS_ROOT}" \
            prep "${NAME_REFERENCE_ARG[@]+"${NAME_REFERENCE_ARG[@]}"}"
else
    echo "[prep] reusing ${RUNS_ROOT}/prepared"
fi

for REP in ${REPLICATES}; do
    REP2=$(printf "%02d" "${REP}")

    sbatch \
        --parsable \
        --job-name="dorsol_${REP2}" \
        --partition="${SHERLOCK_PARTITION}" \
        --time="${SHERLOCK_TIME}" \
        --mem="${SHERLOCK_MEM}" \
        -G 1 \
        --output="${LOG_DIR}/dor_solvent_rep${REP2}_%j.log" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail

module reset >/dev/null 2>&1 || true
module load chemistry py-openmm/8.1.1_py312

cd ${PROJECT_ROOT}

# py-openmm puts openmm on PYTHONPATH, not in the interpreter's site-packages,
# so this must KEEP what module load set. The heredoc is unquoted so REP and
# PROJECT_ROOT resolve at submit time; \${PYTHONPATH:-} is escaped so the JOB
# expands it, not the login node where it is empty. Getting that wrong exported
# ".../src:" over openmm's path and killed 55 of 55 jobs once already.
export PYTHONPATH="$GROUP_HOME/nnrti-mechanisms-repo/src:\${PYTHONPATH:-}"

python3 -c "import openmm, nnrti" || {
    echo "openmm/nnrti not importable on \$(hostname)" >&2
    echo "  PYTHONPATH=\${PYTHONPATH:-<empty>}" >&2
    exit 1
}

nvidia-smi --query-gpu=name,memory.total,memory.used,compute_mode --format=csv,noheader 2>&1 | sed "s/^/[gpu] /" || true

# Assets are built once at submit time. The job never writes them: re-running
# prep would re-minimize and move the start geometry out from under a
# checkpoint written against the old assets.
if [ ! -f "${RUNS_ROOT}/prepared/dor_solvent_system.xml" ]; then
    echo "prepared assets missing under ${RUNS_ROOT}/prepared" >&2
    exit 1
fi

python3 -X faulthandler -m nnrti.cli.run_dor_solvent \\
    --root "${RUNS_ROOT}" \\
    run \\
    --replicate ${REP} \\
    --production-ns ${PRODUCTION_NS} \\
    --analysis-frames ${ANALYSIS_FRAMES}
SBATCH_EOF

    echo "submitted rep ${REP2}"
done
