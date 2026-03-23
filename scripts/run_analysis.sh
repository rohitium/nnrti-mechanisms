#!/usr/bin/env bash
# Run the full holo+apo analysis pipeline using the nnrti-prep conda environment.
# MDAnalysis-dependent scripts (drm_distances, dor_key_contacts, apo/holo comparisons)
# require this env.

PYTHON=~/miniconda3/envs/nnrti-prep/bin/python

set -euo pipefail

run() {
  echo ""
  echo "=== $* ==="
  "$@"
}

HOLO_TABLES_DIR="results/tables/holo"
APO_TABLES_DIR="results/tables/apo"
PLOTS_PNG_DIR="results/plots/png"
KEY_CONTACTS_BY_MUT_DIR="${HOLO_TABLES_DIR}/dor_key_contacts_timeseries_by_mutation"

mkdir -p "${HOLO_TABLES_DIR}" "${APO_TABLES_DIR}" "${PLOTS_PNG_DIR}" "${KEY_CONTACTS_BY_MUT_DIR}"

run "$PYTHON" -m src.analysis.cli.fix_pbc_trajectories --root results/md_runs --in-place
run "$PYTHON" -m src.analysis.cli.audit_pbc_trajectories \
  --root results/md_runs \
  --output-csv results/tables/pbc_audit.csv
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step collect --force
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step metrics --force
run "$PYTHON" -m src.analysis.cli.compute_mmgbsa_safe --snapshots 100 --timestep-fs 2.0 --workers 2
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step plots --force
run "$PYTHON" -m src.analysis.cli.plot_all_mutation_drm_distances \
  --plots-dir "${PLOTS_PNG_DIR}/drm_distances" \
  --output-csv "${HOLO_TABLES_DIR}/drm_sidechain_distance_timeseries_all_mutations.csv" \
  --interesting-csv "${HOLO_TABLES_DIR}/drm_sidechain_distance_interesting_traces.csv"
run "$PYTHON" -m src.analysis.cli.plot_all_mutation_dor_key_contacts \
  --plots-dir "${PLOTS_PNG_DIR}/dor_key_contacts" \
  --output-csv "${HOLO_TABLES_DIR}/dor_key_contacts_timeseries_all_mutations.csv" \
  --contact-defs-csv "${HOLO_TABLES_DIR}/dor_key_contact_definitions_4ncg.csv"

# Split all-mutations CSV into per-mutation CSVs required by the heatmap script.
run "$PYTHON" - <<'PYEOF'
import pandas as pd
from pathlib import Path
src = Path("results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv")
out_dir = Path("results/tables/holo/dor_key_contacts_timeseries_by_mutation")
out_dir.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(src)
for mut, grp in df.groupby("mutation"):
    safe = str(mut).replace("+", "_")
    grp.to_csv(out_dir / f"{safe}_dor_key_contacts_timeseries.csv", index=False)
PYEOF

run "$PYTHON" -m src.analysis.cli.curate_interesting_drm_traces \
  --input-csv "${HOLO_TABLES_DIR}/drm_sidechain_distance_timeseries_all_mutations.csv" \
  --output-csv "${HOLO_TABLES_DIR}/drm_sidechain_distance_interesting_traces.csv" \
  --plots-dir "${PLOTS_PNG_DIR}/drm_distances"
run "$PYTHON" -m src.analysis.cli.plot_interesting_drm_distance_traces \
  --interesting "${HOLO_TABLES_DIR}/drm_sidechain_distance_interesting_traces.csv" \
  --timeseries "${HOLO_TABLES_DIR}/drm_sidechain_distance_timeseries_all_mutations.csv" \
  --output "${PLOTS_PNG_DIR}/interesting_drm_distance_traces.png"
run "$PYTHON" -m src.analysis.cli.plot_mmgbsa_component_signatures \
  --ddg-full "${HOLO_TABLES_DIR}/ddg_full.csv" \
  --output "${PLOTS_PNG_DIR}/manuscript_global_signatures.png"
run "$PYTHON" -m src.analysis.cli.plot_key_contact_occupancy_heatmap \
  --timeseries-dir "${KEY_CONTACTS_BY_MUT_DIR}" \
  --contact-defs "${HOLO_TABLES_DIR}/dor_key_contact_definitions_4ncg.csv" \
  --output "${PLOTS_PNG_DIR}/dor_key_contact_occupancy_heatmap.png" \
  --corr-output "${PLOTS_PNG_DIR}/dor_key_contact_selected_vs_fold_reduction.png" \
  --corr-all-output "${PLOTS_PNG_DIR}/dor_key_contact_all_vs_fold_reduction.png"
run "$PYTHON" -m src.analysis.cli.plot_pocket_volume_distributions \
  --profiles "${HOLO_TABLES_DIR}/pocket_volume_profiles.csv" \
  --output "${PLOTS_PNG_DIR}/pocket_volume_distribution_by_mutation.png"

# Apo/holo comparative analyses (run only when apo manifest is available).
if [[ -f "manifests/apo_md_manifest.csv" ]]; then
  run bash scripts/run_apo_analysis.sh
  run "$PYTHON" -m src.analysis.cli.compute_nnbp_pocket_volume
  run "$PYTHON" -m src.analysis.cli.compute_t290_i63_distance

  # Remove obsolete Y181 chi2 outputs to avoid conflicting interpretations.
  rm -f results/apo_y181_chi2.csv \
        results/tables/apo/apo_y181_chi2.csv \
        results/plots/apo_y181_chi2_timeseries.png \
        results/plots/apo_y181_chi2_distribution.png \
        results/plots/apo_y181_chi2_vs_fold.png \
        results/plots/png/apo/apo_y181_chi2_timeseries.png \
        results/plots/png/apo/apo_y181_chi2_distribution.png \
        results/plots/png/apo/apo_y181_chi2_vs_fold.png
else
  echo ""
  echo "=== Skipping apo analyses (manifests/apo_md_manifest.csv not found) ==="
fi

# Consolidate any root-level CSV/PNG outputs into the structured results layout.
run "$PYTHON" - <<'PYEOF'
from pathlib import Path
import shutil

tables_holo = Path("results/tables/holo")
tables_holo.mkdir(parents=True, exist_ok=True)
png_root = Path("results/plots/png")
png_root.mkdir(parents=True, exist_ok=True)

for csv in Path("results").glob("*.csv"):
    target = tables_holo / csv.name
    if target.exists():
        target.unlink()
    shutil.move(str(csv), str(target))

src_plots = Path("results/plots")
for png in src_plots.rglob("*.png"):
    if png.is_dir() or "results/plots/png/" in str(png):
        continue
    rel_parent = png.parent.relative_to(src_plots)
    dest_dir = png_root / rel_parent
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / png.name
    if target.exists():
        target.unlink()
    shutil.move(str(png), str(target))
PYEOF

echo ""
echo "=== run_analysis.sh complete ==="
