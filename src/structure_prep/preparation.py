from __future__ import annotations

import logging
import shutil
from pathlib import Path

from ..analysis.susceptibility import load_dor_susceptibilities
from ..md.manifest import MDTask, save_manifest
from ..md.openmm.md_protocol import MDProtocolConfig, prepare_md_assets
from ..utils import ensure_dirs, load_chain_subunits, load_residue_mappings, sanitize_label
from ..utils.mutations import deterministic_seed
from .config import dor_4ncg_spec
from .mutation.mutagenesis import apply_mutations
from .mutation.numbering import detect_numbering_scheme
from .mutation.steps import build_mutation_steps


def prepare_local_openmm_only_for_cluster(
    root: Path,
    susceptibility_xlsx: Path,
    prepared_dir: Path,
    manifest_path: Path,
    results_dir: Path,
    replicates: int = 3,
    jitter_seed_base: int | None = None,
    jitter_angstrom: float = 0.1,
    md_config: MDProtocolConfig | None = None,
    selected_mutations: set[str] | None = None,
) -> list[MDTask]:
    """Prepare WT/mutant systems for Sherlock explicit MD (no alchemical protocol)."""
    from ..md.openmm.minimizer import minimize_system
    from ..md.openmm.require import require_module
    from ..md.openmm.structure import minimize_with_restraints

    cfg = md_config or MDProtocolConfig()

    def _norm_mutation(label: str) -> str:
        return "+".join(
            token.strip().upper().replace(" ", "")
            for token in str(label).replace(",", "+").split("+")
            if token.strip()
        )

    spec = dor_4ncg_spec(root)
    dor_df = load_dor_susceptibilities(susceptibility_xlsx, default_chain="A")
    chain_map = load_chain_subunits(spec.structure.cif_path)
    residue_maps = load_residue_mappings(spec.structure.cif_path)
    numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)

    selected = {_norm_mutation(m) for m in selected_mutations} if selected_mutations else None
    if selected is not None:
        dor_df = dor_df[dor_df["mutation"].map(_norm_mutation).isin(selected)].copy()
        if dor_df.empty:
            raise ValueError(f"No rows matched selected mutations: {', '.join(sorted(selected))}")

    ensure_dirs([prepared_dir, manifest_path.parent, results_dir])

    wt_cif = prepared_dir / "wt_4ncg.cif"
    if not wt_cif.exists():
        shutil.copy2(spec.structure.cif_path, wt_cif)

    structures: list[dict] = [
        {"mutation": "WT", "safe_label": "wt", "cif_path": wt_cif, "fold_reduction": None}
    ]
    for _, row in dor_df.iterrows():
        mutation = row["mutation"]
        safe_label = sanitize_label(mutation)
        chain_list = [c.strip().upper() for c in str(row["chain"]).split("+") if c.strip()]
        mutation_steps, _ = build_mutation_steps(
            mutation_label=mutation,
            chain_list=chain_list,
            residue_maps=residue_maps,
            numbering_scheme=numbering,
        )
        mut_cif = prepared_dir / f"mut_{safe_label}.cif"
        if not mut_cif.exists():
            apply_mutations(spec.structure.cif_path, mutation_steps, mut_cif)
        structures.append(
            {
                "mutation": mutation,
                "safe_label": safe_label,
                "cif_path": mut_cif,
                "fold_reduction": float(row["dor_fold_reduction"]),
            }
        )

    logging.info("Preparing MD assets for %d structures (%d replicates)", len(structures), replicates)

    tasks: list[MDTask] = []
    task_id = 0
    app = require_module("openmm.app")

    for struct in structures:
        safe_label = struct["safe_label"]
        for replicate in range(1, replicates + 1):
            run_dir = results_dir / safe_label / f"rep_{replicate:02d}"
            assets_dir = run_dir / "assets"
            ensure_dirs([run_dir, assets_dir])

            seed = deterministic_seed(jitter_seed_base, safe_label, replicate)

            min_pdb = run_dir / f"{safe_label}_minimized_rep{replicate:02d}.pdb"
            if not min_pdb.exists():
                topology, positions, forcefield = minimize_with_restraints(
                    cif_path=Path(struct["cif_path"]),
                    ligand_resname=spec.structure.ligand_resname,
                    ligand_sdf=spec.structure.ligand_sdf,
                    restraint_radius_angstrom=spec.restraint_radius_angstrom,
                    restraint_k_kj_mol_nm2=spec.restraint_k_kj_mol_nm2,
                    output_path=min_pdb,
                    jitter_seed=seed,
                    jitter_angstrom=jitter_angstrom,
                )
                _, positions = minimize_system(
                    topology,
                    positions,
                    forcefield,
                    restraint_indices=[],
                    restraint_k_kj_mol_nm2=0.0,
                )
                with open(min_pdb, "w") as handle:
                    app.PDBFile.writeFile(topology, positions, handle)

            system_xml = assets_dir / f"{safe_label}_md_rep{replicate:02d}_system.xml"
            topology_pdb = assets_dir / f"{safe_label}_md_rep{replicate:02d}_start.pdb"

            if not (system_xml.exists() and topology_pdb.exists()):
                prepare_md_assets(
                    minimized_pdb_path=min_pdb,
                    ligand_resname=spec.structure.ligand_resname,
                    ligand_sdf=spec.structure.ligand_sdf,
                    topology_pdb_path=topology_pdb,
                    system_xml_path=system_xml,
                    config=cfg,
                )

            output_json = run_dir / f"{safe_label}_rep{replicate:02d}.json"
            tasks.append(
                MDTask(
                    task_id=task_id,
                    structure="DOR",
                    mutation=struct["mutation"],
                    safe_label=safe_label,
                    replicate=replicate,
                    minimized_pdb=str(min_pdb.resolve()),
                    ligand_sdf=str(spec.structure.ligand_sdf.resolve()),
                    ligand_resname=spec.structure.ligand_resname,
                    fold_reduction=struct["fold_reduction"],
                    output_json=str(output_json),
                    input_cif=str(Path(struct["cif_path"]).resolve()),
                    jitter_seed=seed,
                    jitter_angstrom=jitter_angstrom,
                    restraint_radius=spec.restraint_radius_angstrom,
                    restraint_k=spec.restraint_k_kj_mol_nm2,
                    prepared_topology_pdb=str(topology_pdb.resolve()),
                    prepared_system_xml=str(system_xml.resolve()),
                )
            )
            task_id += 1

    save_manifest(tasks, manifest_path)
    logging.info("Wrote MD manifest with %d tasks to %s", len(tasks), manifest_path)
    return tasks
