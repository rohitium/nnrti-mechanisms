"""Snakemake script: collect results, run MM/GBSA + structural metrics + aggregation."""

import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
logging.getLogger("MDAnalysis.analysis").setLevel(logging.WARNING)

from src.md.manifest import MDTask, save_manifest
from src.analysis.result_collector import run_result_collection
from src.analysis.susceptibility import load_dor_susceptibilities

root = Path(".").resolve()

# Build fold-reduction lookup from susceptibility data.
xlsx = Path(snakemake.input.xlsx)  # noqa: F821
dor_df = load_dor_susceptibilities(xlsx, default_chain="A")
fold_lookup = dict(
    zip(dor_df["mutation"].tolist(), dor_df["dor_fold_reduction"].tolist())
)

# Build manifest from completed JSON files.
tasks = []
for i, json_path in enumerate(snakemake.input.jsons):  # noqa: F821
    jp = Path(json_path)
    data = json.loads(jp.read_text())

    mutation = data.get("mutation", "")
    fold_reduction = data.get("fold_reduction")
    if fold_reduction is None and mutation != "WT":
        fold_reduction = fold_lookup.get(mutation)

    tasks.append(
        MDTask(
            task_id=i,
            structure=data.get("structure", "DOR"),
            mutation=mutation,
            safe_label=data.get("safe_label", ""),
            replicate=int(data.get("replicate", 0)),
            minimized_pdb=data.get("minimized_pdb", ""),
            ligand_sdf=data.get("ligand_sdf", str(root / "data" / "ligands" / "dor.sdf")),
            ligand_resname=data.get("ligand_resname", snakemake.params.ligand_resname),  # noqa: F821
            fold_reduction=fold_reduction,
            output_json=str(jp),
            prepared_topology_pdb=data.get("prepared_topology_pdb", ""),
            prepared_system_xml=data.get("prepared_system_xml", ""),
        )
    )

# Write temporary manifest for the collection pipeline.
manifest_path = root / "results" / ".snakemake_manifest.csv"
save_manifest(tasks, manifest_path)

# Run the existing collection pipeline.
run_result_collection(
    manifest_path=manifest_path,
    md_results_dir=root / "results" / "md_runs",
    output_dir=root / "results",
    ligand_resname=snakemake.params.ligand_resname,  # noqa: F821
    compute_structural=True,
    metric_frame_stride=max(1, snakemake.params.metric_frame_stride),  # noqa: F821
    metric_max_frames=max(1, snakemake.params.metric_max_frames),  # noqa: F821
    mmgbsa_snapshots=max(5, snakemake.params.mmgbsa_snapshots),  # noqa: F821
    mmgbsa_discard_fraction=max(0.0, min(0.9, snakemake.params.mmgbsa_discard_fraction)),  # noqa: F821
)
