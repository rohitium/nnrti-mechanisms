#!/bin/bash
#
# Interactive GPU smoke test: Y188L apo MD (0.01 ns) before batch submit.
#
# Usage:
#   cd $PROJECT_ROOT
#   git pull origin main
#   bash ops/slurm/fep/salloc_apo_gpu.sh          # get gpu node
#   bash ops/slurm/fep/test_y188l_apo_gpu.sh      # run inside allocation
#
# Optional:
#   REP=02 bash ops/slurm/fep/test_y188l_apo_gpu.sh
#   MD_PRODUCTION_NS=0.05 bash ops/slurm/fep/test_y188l_apo_gpu.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
cd "$PROJECT_ROOT"
export PYTHONPATH="${PROJECT_ROOT:-$PWD}/src:${PYTHONPATH:-}"

source ops/slurm/cluster/load_openmm_module.sh

MUTATION="y188l"
REP="${REP:-01}"
REP_INT=$((10#$REP))
PRODUCTION_NS="${MD_PRODUCTION_NS:-0.01}"

PARENT="results/md_runs/apo/${MUTATION}/rep_${REP}"
ASSETS="${PARENT}/assets"
SYSTEM_XML="${ASSETS}/${MUTATION}_apo_md_rep${REP}_system.xml"
TOPOLOGY_PDB="${ASSETS}/${MUTATION}_apo_md_rep${REP}_start.pdb"
HOLO_REP="results/md_runs/Y188L/rep_${REP}"
MINIMIZED_PDB="${HOLO_REP}/Y188L_minimized_rep${REP}.pdb"
if [[ ! -f "$MINIMIZED_PDB" ]]; then
    # run_prepared_md ignores this; only stored in output JSON metadata
    MINIMIZED_PDB="$TOPOLOGY_PDB"
fi
OUTPUT_JSON="${PARENT}/${MUTATION}_apo_rep${REP}_TEST.json"

echo "=========================================="
echo "Y188L apo GPU smoke test"
echo "=========================================="
echo "Project:       $PROJECT_ROOT"
echo "Replicate:     $REP"
echo "Production ns: $PRODUCTION_NS"
echo "System XML:    $SYSTEM_XML"
echo "Topology PDB:  $TOPOLOGY_PDB"
echo "Output JSON:   $OUTPUT_JSON (test only)"
echo ""

python3 - <<'PY'
import openmm
from openmm import Platform
print("OpenMM:", openmm.__version__)
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f"  platform {i}: {p.getName()}", end="")
    try:
        print(f" (speed {p.getSpeed()})")
    except Exception:
        print()
PY

MISSING=0
for f in "$SYSTEM_XML" "$TOPOLOGY_PDB"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing $f" >&2
        MISSING=1
    else
        echo "✓ $f"
    fi
done
if [[ -f "${HOLO_REP}/Y188L_minimized_rep${REP}.pdb" ]]; then
    echo "✓ ${HOLO_REP}/Y188L_minimized_rep${REP}.pdb"
else
    echo "⚠ holo minimized PDB missing; using apo topology for JSON metadata only"
fi
if [[ "$MISSING" -eq 1 ]]; then
    echo ""
    echo "Sync apo assets to Sherlock if needed, e.g.:"
    echo "  rsync -av results/md_runs/apo/y188l/ \\"
    echo "    <user>@<cluster>:\$PROJECT_ROOT/results/md_runs/apo/y188l/"
    exit 1
fi

echo ""
echo "Running apo MD smoke test..."
python3 -m nnrti.md.sherlock.run_md_job \
    --mutation "y188l" \
    --replicate "$REP_INT" \
    --task-id 0 \
    --system-xml "$SYSTEM_XML" \
    --topology-pdb "$TOPOLOGY_PDB" \
    --minimized-pdb "$MINIMIZED_PDB" \
    --output-json "$OUTPUT_JSON" \
    --ligand-sdf "" \
    --ligand-resname "" \
    --production-ns "$PRODUCTION_NS" \
    --no-resume \
    --force

echo ""
echo "=========================================="
echo "✓ Smoke test finished"
echo "=========================================="
cat "$OUTPUT_JSON"
echo ""
echo "If status=ok, submit batch:"
echo "  bash ops/slurm/fep/submit_y188l_apo_md.sh"
echo ""
echo "Remove test artifacts (optional):"
echo "  rm -f ${OUTPUT_JSON} ${PARENT}/${MUTATION}_rep${REP}_*"
