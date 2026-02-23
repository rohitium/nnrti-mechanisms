#!/usr/bin/env bash
# Run the apo trajectory analysis pipeline.
#
# Prerequisites:
#   1. Apo systems prepared:  python -m src.md.dor_md_pipeline_apo
#   2. Apo MD completed on Sherlock and synced locally.
#
# This script:
#   1. Fixes PBC imaging artifacts in apo *_analysis.dcd files.
#   2. Computes NNBP tunnel gate distances (Hypothesis 2: tunnel opening).
#      - All apo mutations, output to results/apo_nnbp_tunnel_{dynamics,summary}.csv
#   3. Computes DCCM allosteric coupling scores (Hypothesis 1: RT processivity).
#      - Key mutations for H1: WT, F227C, V106A, A98G+F227C, V106I+F227C
#      - Output to results/apo_dccm_allosteric_coupling.csv
#   4. (Combined) Re-runs tunnel dynamics on the HOLO manifest for the same
#      mutations so apo vs holo gate distances can be compared directly.

PYTHON=~/miniconda3/envs/nnrti-prep/bin/python
APO_MANIFEST="results/apo_md_manifest.csv"
HOLO_MANIFEST="results/md_manifest.csv"

set -euo pipefail

# ── Step 1: PBC correction ────────────────────────────────────────────────────
echo "=== Step 1: PBC correction (apo trajectories) ==="
$PYTHON -m src.analysis.cli.fix_pbc_trajectories \
    --root results/apo_md_runs \
    --in-place

# ── Step 2: Tunnel dynamics (apo) ────────────────────────────────────────────
echo ""
echo "=== Step 2: NNBP tunnel gate distances — apo ==="
$PYTHON -m src.analysis.cli.compute_nnbp_tunnel_dynamics \
    --manifest "${APO_MANIFEST}" \
    --output-csv results/apo_nnbp_tunnel_dynamics.csv \
    --summary-csv results/apo_nnbp_tunnel_summary.csv \
    --plots-dir results/plots/apo_nnbp_tunnel

# ── Step 3: Tunnel dynamics (holo, same mutations for direct comparison) ─────
# Only needed if not already computed from run_analysis.sh.
echo ""
echo "=== Step 3: NNBP tunnel gate distances — holo (comparison baseline) ==="
$PYTHON -m src.analysis.cli.compute_nnbp_tunnel_dynamics \
    --manifest "${HOLO_MANIFEST}" \
    --mutations WT F227C V106A "V106A+P225H" "K103N+M230L" "A98G+F227C" "V106I+F227C" \
    --output-csv results/holo_nnbp_tunnel_dynamics.csv \
    --summary-csv results/holo_nnbp_tunnel_summary.csv \
    --plots-dir results/plots/holo_nnbp_tunnel

# ── Step 4: DCCM (apo) — Hypothesis 1 mutations ──────────────────────────────
echo ""
echo "=== Step 4: DCCM allosteric coupling — apo (H1 mutations) ==="
$PYTHON -m src.analysis.cli.compute_dccm \
    --manifest "${APO_MANIFEST}" \
    --mutations WT F227C V106A "A98G+F227C" "V106I+F227C" \
    --output-csv results/apo_dccm_allosteric_coupling.csv \
    --plots-dir results/plots/apo_dccm \
    --save-plots

# ── Step 5: DCCM (holo) — same mutations for direct comparison ───────────────
echo ""
echo "=== Step 5: DCCM allosteric coupling — holo (H1 comparison baseline) ==="
$PYTHON -m src.analysis.cli.compute_dccm \
    --manifest "${HOLO_MANIFEST}" \
    --mutations WT F227C V106A "A98G+F227C" "V106I+F227C" \
    --output-csv results/holo_dccm_allosteric_coupling.csv \
    --plots-dir results/plots/holo_dccm \
    --save-plots

echo ""
echo "=== Apo analysis complete ==="
echo ""
echo "Key outputs:"
echo "  results/apo_nnbp_tunnel_summary.csv     — gate mean±std per replicate (apo)"
echo "  results/holo_nnbp_tunnel_summary.csv    — gate mean±std per replicate (holo)"
echo "  results/apo_dccm_allosteric_coupling.csv  — NNBP↔domain coupling (apo)"
echo "  results/holo_dccm_allosteric_coupling.csv — NNBP↔domain coupling (holo)"
echo "  results/plots/apo_nnbp_tunnel/          — apo gate distance plots"
echo "  results/plots/apo_dccm/                 — apo DCCM heatmaps"
echo ""
echo "To compare apo vs holo, load the two summary CSVs and inspect the"
echo "gate_K101_Y188_CA and NNBP_fingers coupling columns."
