#!/bin/bash
# Pose-independence test (Reviewer 1, comment 2): four WT RT-DOR systems that
# differ only in the DOR starting pose. 2 ns triage, one replicate each.
#
# See results/md_runs/pose_test/HANDOFF.md for provenance and how to read the
# result. Short jobs by design, so Slurm can backfill them around the apo runs.
set -euo pipefail

ROOT="${PROJECT_ROOT:-/scratch/users/rsatija/nnrti-mechanisms-git}"
SRC="${NNRTI_SRC:-$GROUP_HOME/nnrti-mechanisms/src}"
TIME="${SHERLOCK_TIME:-03:00:00}"
NS="${PRODUCTION_NS:-2}"
LOGS="$ROOT/logs/pose_test"; mkdir -p "$LOGS"
cd "$ROOT"

for POSE in pose_alt_01 pose_alt_02 pose_alt_05 pose_tilt_02; do
    D="results/md_runs/pose_test/$POSE/rep_01"
    if squeue -u "$USER" -h -o "%j" | grep -qx "pt_$POSE"; then
        echo "  $POSE already queued, skipping"; continue
    fi
    sbatch -J "pt_$POSE" -p gpu -G 1 -c 4 --mem=16G -t "$TIME" \
           -o "$LOGS/${POSE}_%j.log" <<SBATCH_EOF
#!/bin/bash
set -euo pipefail
module reset >/dev/null 2>&1 || true
module load chemistry py-openmm/8.1.1_py312
cd ${ROOT}
# Escaped: an unquoted heredoc would expand this on the login node, where it is
# empty, erasing py-openmm's path. That killed 55 apo jobs.
export PYTHONPATH="${SRC}:\${PYTHONPATH:-}"
python3 -c "import openmm, nnrti" || {
    echo "openmm/nnrti not importable on \$(hostname)" >&2
    echo "  PYTHONPATH=\${PYTHONPATH:-<empty>}" >&2
    exit 1
}
nvidia-smi --query-gpu=name,memory.used --format=csv,noheader 2>&1 | sed "s/^/[gpu] /" || true

python3 -X faulthandler -m nnrti.md.sherlock.run_md_job \
    --mutation "${POSE}" --replicate 1 --task-id 0 \
    --system-xml    "${D}/assets/${POSE}_md_rep01_system.xml" \
    --topology-pdb  "${D}/assets/${POSE}_md_rep01_start.pdb" \
    --minimized-pdb "${D}/${POSE}_minimized_rep01.pdb" \
    --output-json   "${D}/${POSE}_rep01.json" \
    --ligand-sdf "data/ligands/dor.sdf" --ligand-resname "2KW" \
    --heating-ps 25 --production-ns ${NS} \
    --report-interval 2000 --checkpoint-interval 5000
SBATCH_EOF
    echo "  submitted $POSE"
    sleep 0.5
done
