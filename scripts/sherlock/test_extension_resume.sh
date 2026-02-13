#!/bin/bash
#
# Smoke test for checkpoint-based trajectory extension on Sherlock.
#
# This validates the exact path we need for 2 ns -> 10 ns:
# - existing run with checkpoint
# - force rerun with a slightly higher production target
# - confirm resume happened and target steps were reached
#
# Usage:
#   # Run inside an sh_dev allocation/session
#   # Optional args: [mutation] [rep]
#   bash scripts/sherlock/test_extension_resume.sh [mutation] [rep]
#
# Optional env vars:
#   TARGET_NS                     (optional override; if unset, uses current_ns + 0.01)
#   BOOTSTRAP_NS                  (default: 0.02, only used when no source checkpoint exists)
#   BOOTSTRAP_CHECKPOINT_INTERVAL (default: 1000 steps, only used in bootstrap run)
#   KEEP_SMOKE_DIR                (default: 0; set to 1 to keep temp outputs)
#   OPENMM_PLATFORM               (optional; CUDA/CPU/OpenCL)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

TARGET_NS="${TARGET_NS:-}"
BOOTSTRAP_NS="${BOOTSTRAP_NS:-0.02}"
BOOTSTRAP_CHECKPOINT_INTERVAL="${BOOTSTRAP_CHECKPOINT_INTERVAL:-1000}"
KEEP_SMOKE_DIR="${KEEP_SMOKE_DIR:-0}"
MUTATION="${1:-}"
REP="${2:-}"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required for this smoke test." >&2
  exit 1
fi

source /etc/profile.d/modules.sh 2>/dev/null || true
module load chemistry py-openmm/8.1.1_py312

pick_candidate() {
  local chosen=""
  local best_score=-1
  local j
  for j in $(find results/md_runs -name '*_rep*.json' | sort); do
    local d base chk steps status score has_chk
    d=$(dirname "$j")
    base=$(basename "$j" .json)
    chk="$d/${base}_md.chk"

    # Keep only canonical replicate JSON names (e.g., wt_rep01.json)
    if [[ ! "$base" =~ _rep[0-9][0-9]$ ]]; then
      continue
    fi

    steps=$(jq -r '.md_production_steps_completed // .md_production_steps // 0' "$j" 2>/dev/null || echo "0")
    status=$(jq -r '.status // ""' "$j" 2>/dev/null || echo "")

    # Score candidates: prefer status=ok; then existing checkpoint; then larger step counts.
    score=0
    if [ "$status" = "ok" ]; then
      score=$((score + 1000000000))
    fi
    has_chk=0
    if [ -f "$chk" ]; then
      has_chk=1
      score=$((score + 100000000))
    fi
    if [ "${steps}" -gt 0 ] 2>/dev/null; then
      score=$((score + steps))
    fi

    if [ "$score" -gt "$best_score" ]; then
      best_score="$score"
      chosen="$j"
    fi
  done
  echo "$chosen"
}

if [ -n "$MUTATION" ] && [ -n "$REP" ]; then
  CANDIDATE_JSON="results/md_runs/${MUTATION}/rep_${REP}/${MUTATION}_rep${REP}.json"
else
  CANDIDATE_JSON="$(pick_candidate)"
fi

if [ -z "$CANDIDATE_JSON" ] || [ ! -f "$CANDIDATE_JSON" ]; then
  echo "ERROR: Could not find a suitable replicate JSON candidate." >&2
  echo "Debug: json count=$(find results/md_runs -name '*_rep*.json' | wc -l | tr -d ' '), chk count=$(find results/md_runs -name '*_md.chk' | wc -l | tr -d ' ')" >&2
  echo "Tip: pass explicit args: bash $0 <mutation> <rep>" >&2
  exit 1
fi

REP_DIR="$(dirname "$CANDIDATE_JSON")"
MUTATION="$(basename "$(dirname "$REP_DIR")")"
REP_DIR_BASENAME="$(basename "$REP_DIR")"
REP="${REP_DIR_BASENAME#rep_}"
REP_INT=$((10#$REP))

ASSETS_DIR="${REP_DIR}/assets"
SYSTEM_XML="${ASSETS_DIR}/${MUTATION}_md_rep${REP}_system.xml"
TOPOLOGY_PDB="${ASSETS_DIR}/${MUTATION}_md_rep${REP}_start.pdb"
MINIMIZED_PDB="${REP_DIR}/${MUTATION}_minimized_rep${REP}.pdb"
SOURCE_CHK="${REP_DIR}/${MUTATION}_rep${REP}_md.chk"

for f in "$SYSTEM_XML" "$TOPOLOGY_PDB" "$MINIMIZED_PDB"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f" >&2
    exit 1
  fi
done

SMOKE_DIR="$(mktemp -d "/tmp/nnrti_ext_smoke_${MUTATION}_rep${REP}_XXXXXX")"
BOOTSTRAP_JSON="${SMOKE_DIR}/${MUTATION}_rep${REP}_bootstrap.json"
OUTPUT_JSON="${SMOKE_DIR}/${MUTATION}_rep${REP}_smoketest.json"
TEST_CHK="${SMOKE_DIR}/${MUTATION}_rep${REP}_md.chk"
CURRENT_STEPS=0
RESUME_SOURCE="bootstrap"

if [ -f "$SOURCE_CHK" ]; then
  cp "$SOURCE_CHK" "$TEST_CHK"
  CURRENT_STEPS="$(jq -r '.md_production_steps_completed // .md_production_steps // 0' "$CANDIDATE_JSON" 2>/dev/null || echo "0")"
  RESUME_SOURCE="existing"
else
  echo "No source checkpoint found for candidate; bootstrapping a short run to create one..."
  python3 -m src.md.sherlock.run_md_job \
    --mutation "$MUTATION" \
    --replicate "$REP_INT" \
    --task-id 999998 \
    --system-xml "$SYSTEM_XML" \
    --topology-pdb "$TOPOLOGY_PDB" \
    --minimized-pdb "$MINIMIZED_PDB" \
    --output-json "$BOOTSTRAP_JSON" \
    --production-ns "$BOOTSTRAP_NS" \
    --checkpoint-interval "$BOOTSTRAP_CHECKPOINT_INTERVAL" \
    --no-resume \
    --force

  if [ ! -f "$TEST_CHK" ]; then
    echo "ERROR: Bootstrap run did not create checkpoint: $TEST_CHK" >&2
    exit 1
  fi
  CURRENT_STEPS="$(jq -r '.md_production_steps_completed // .md_production_steps // 0' "$BOOTSTRAP_JSON" 2>/dev/null || echo "0")"
fi

if ! [[ "$CURRENT_STEPS" =~ ^[0-9]+$ ]]; then
  CURRENT_STEPS=0
fi

if [ -n "$TARGET_NS" ]; then
  TARGET_STEPS="$(python3 - <<PY
ns = float("${TARGET_NS}")
print(max(1, int(round((ns * 1_000_000.0) / 2.0))))
PY
)"
else
  # Default smoke increment = +0.01 ns = +5000 steps at 2 fs.
  TARGET_STEPS=$((CURRENT_STEPS + 5000))
fi

TARGET_NS_EFFECTIVE="$(python3 - <<PY
steps = int("${TARGET_STEPS}")
print(f"{(steps * 2.0) / 1_000_000.0:.6f}")
PY
)"

echo "=========================================="
echo "Checkpoint Extension Smoke Test"
echo "=========================================="
echo "Candidate JSON: $CANDIDATE_JSON"
echo "Mutation:       $MUTATION"
echo "Replicate:      $REP"
echo "Resume source:  $RESUME_SOURCE checkpoint"
echo "Current steps:  $CURRENT_STEPS"
echo "Target ns:      $TARGET_NS_EFFECTIVE"
echo "Target steps:   $TARGET_STEPS"
echo "Smoke dir:      $SMOKE_DIR"
echo ""

python3 -m src.md.sherlock.run_md_job \
  --mutation "$MUTATION" \
  --replicate "$REP_INT" \
  --task-id 999999 \
  --system-xml "$SYSTEM_XML" \
  --topology-pdb "$TOPOLOGY_PDB" \
  --minimized-pdb "$MINIMIZED_PDB" \
  --output-json "$OUTPUT_JSON" \
  --production-ns "$TARGET_NS_EFFECTIVE" \
  --resume \
  --force

if [ ! -f "$OUTPUT_JSON" ]; then
  echo "ERROR: Smoke test output JSON not created." >&2
  exit 1
fi

STATUS="$(jq -r '.status // ""' "$OUTPUT_JSON")"
RESUMED="$(jq -r '.resumed_from_checkpoint // false' "$OUTPUT_JSON")"
COMPLETED_STEPS="$(jq -r '.md_production_steps_completed // .md_production_steps // 0' "$OUTPUT_JSON")"

echo ""
echo "Result summary:"
echo "  status:                  $STATUS"
echo "  resumed_from_checkpoint: $RESUMED"
echo "  production_steps_done:   $COMPLETED_STEPS"

if [ "$STATUS" != "ok" ]; then
  echo "ERROR: Smoke test status is not ok." >&2
  exit 1
fi

if [ "$RESUMED" != "true" ]; then
  echo "ERROR: Run did not resume from checkpoint." >&2
  exit 1
fi

if [ "$COMPLETED_STEPS" -lt "$TARGET_STEPS" ]; then
  echo "ERROR: Completed steps ($COMPLETED_STEPS) < target steps ($TARGET_STEPS)." >&2
  exit 1
fi

echo ""
echo "PASS: checkpoint extension path is working."

if [ "$KEEP_SMOKE_DIR" = "1" ]; then
  echo "Preserving smoke outputs: $SMOKE_DIR"
else
  rm -rf "$SMOKE_DIR"
  echo "Cleaned smoke outputs."
fi
