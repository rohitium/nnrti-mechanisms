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

run "$PYTHON" -m src.analysis.cli.fix_pbc_trajectories --root results/md_runs --in-place
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step collect --force
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step metrics --force
run "$PYTHON" -m src.analysis.cli.compute_mmgbsa_safe --snapshots 100 --sample-window-ns 1.0 --timestep-fs 2.0 --workers 2
run "$PYTHON" -m src.analysis.cli.analyze_incremental --step plots --force
run "$PYTHON" -m src.analysis.cli.plot_all_mutation_drm_distances
run "$PYTHON" -m src.analysis.cli.plot_all_mutation_dor_key_contacts

# Split all-mutations CSV into per-mutation CSVs required by the heatmap script.
run "$PYTHON" - <<'PYEOF'
import pandas as pd
from pathlib import Path
src = Path("results/dor_key_contacts_timeseries_all_mutations.csv")
out_dir = Path("results/dor_key_contacts_timeseries_by_mutation")
out_dir.mkdir(exist_ok=True)
df = pd.read_csv(src)
for mut, grp in df.groupby("mutation"):
    safe = str(mut).replace("+", "_")
    grp.to_csv(out_dir / f"{safe}_dor_key_contacts_timeseries.csv", index=False)
PYEOF

run "$PYTHON" -m src.analysis.cli.curate_interesting_drm_traces
run "$PYTHON" -m src.analysis.cli.plot_interesting_drm_distance_traces
run "$PYTHON" -m src.analysis.cli.plot_mmgbsa_component_signatures
run "$PYTHON" -m src.analysis.cli.plot_key_contact_occupancy_heatmap
run "$PYTHON" -m src.analysis.cli.plot_pocket_volume_distributions

# Apo/holo comparative analyses (run only when apo manifest is available).
if [[ -f "results/apo_md_manifest.csv" ]]; then
  run bash scripts/run_apo_analysis.sh
  run "$PYTHON" scripts/compute_nnbp_pocket_volume.py
  run "$PYTHON" scripts/compute_t290_i63_distance.py

  # Remove obsolete Y181 chi2 outputs to avoid conflicting interpretations.
  rm -f results/apo_y181_chi2.csv \
        results/plots/apo_y181_chi2_timeseries.png \
        results/plots/apo_y181_chi2_distribution.png \
        results/plots/apo_y181_chi2_vs_fold.png
else
  echo ""
  echo "=== Skipping apo analyses (results/apo_md_manifest.csv not found) ==="
fi

echo ""
echo "=== run_analysis.sh complete ==="
