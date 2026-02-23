#!/usr/bin/env bash
# Run the full analysis pipeline using the nnrti-prep conda environment.
# MDAnalysis-dependent scripts (drm_distances, dor_key_contacts) require this env.

PYTHON=~/miniconda3/envs/nnrti-prep/bin/python

set -euo pipefail

$PYTHON -m src.analysis.cli.fix_pbc_trajectories --root results/md_runs --in-place
$PYTHON -m src.analysis.cli.analyze_incremental --step collect --force
$PYTHON -m src.analysis.cli.analyze_incremental --step metrics --force
$PYTHON -m src.analysis.cli.compute_mmgbsa_safe --snapshots 100 --sample-window-ns 1.0 --timestep-fs 2.0 --workers 2
$PYTHON -m src.analysis.cli.analyze_incremental --step plots --force
$PYTHON -m src.analysis.cli.plot_all_mutation_drm_distances
$PYTHON -m src.analysis.cli.plot_all_mutation_dor_key_contacts

# Split all-mutations CSV into per-mutation CSVs required by the heatmap script.
$PYTHON - <<'PYEOF'
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

$PYTHON -m src.analysis.cli.curate_interesting_drm_traces
$PYTHON -m src.analysis.cli.plot_interesting_drm_distance_traces
$PYTHON -m src.analysis.cli.plot_mmgbsa_component_signatures
$PYTHON -m src.analysis.cli.plot_key_contact_occupancy_heatmap
$PYTHON -m src.analysis.cli.plot_pocket_volume_distributions
