from __future__ import annotations

import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from .cluster.manifest import FEPTask, save_manifest
from .config import dor_4ncg_spec
from .mutation.rows import build_rows
from .mutation.steps import build_mutation_steps
from .mutation.mutagenesis import apply_mutations
from .numbering import detect_numbering_scheme
from .openmm.alchemy import AlchemicalConfig
from .structure_prep import prepare_structure_local
from .susceptibility_io import load_dor_susceptibilities
from .utils import (
    ensure_dirs,
    load_chain_subunits,
    load_residue_mappings,
    sanitize_label,
)


def prepare_dor_local_structures(
    root: Path,
    susceptibility_xlsx: Path,
    prepared_dir: Path,
    manifest_csv: Path,
) -> pd.DataFrame:
    """Prepare DOR mutation CIF structures for cluster runs (legacy mode)."""
    spec = dor_4ncg_spec(root)
    dor_df = load_dor_susceptibilities(susceptibility_xlsx, default_chain="A")
    chain_map = load_chain_subunits(spec.structure.cif_path)
    residue_maps = load_residue_mappings(spec.structure.cif_path)
    numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)

    ensure_dirs([prepared_dir, manifest_csv.parent])
    wt_cif = prepared_dir / "wt_4ncg.cif"
    if not wt_cif.exists():
        shutil.copy2(spec.structure.cif_path, wt_cif)

    manifest_rows: list[dict] = []
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
        manifest_rows.append(
            {
                "structure": "DOR",
                "mutation": mutation,
                "safe_label": safe_label,
                "chain": "+".join(chain_list),
                "prepared_cif": str(mut_cif.resolve()),
                "wt_cif": str(wt_cif.resolve()),
                "dor_fold_reduction": float(row["dor_fold_reduction"]),
                "mutation_order": int(row["order"]),
            }
        )

    manifest = pd.DataFrame(manifest_rows).sort_values("mutation_order").reset_index(
        drop=True
    )
    manifest.to_csv(manifest_csv, index=False)
    return manifest


def prepare_local_with_fep_manifest(
    root: Path,
    susceptibility_xlsx: Path,
    prepared_dir: Path,
    fep_manifest_path: Path,
    structural_metrics_path: Path,
    fep_results_dir: Path,
    replicates: int = 3,
    jitter_seed_base: int | None = None,
    jitter_angstrom: float = 0.1,
) -> tuple[pd.DataFrame, list[FEPTask]]:
    """Prepare structures locally and generate FEP manifest for cluster runs.

    This function:
    1. Loads mutation data from susceptibility spreadsheet
    2. Creates mutant CIF files using PDBFixer
    3. Minimizes all structures (WT + mutants) with restraints
    4. Computes structural metrics (contacts, H-bonds, pocket volume)
    5. Generates FEP manifest with task assignments for cluster execution

    Args:
        root: Project root directory.
        susceptibility_xlsx: Path to DOR susceptibility workbook.
        prepared_dir: Directory for prepared/minimized PDB files.
        fep_manifest_path: Path to write FEP manifest CSV.
        structural_metrics_path: Path to write structural metrics CSV.
        fep_results_dir: Directory where FEP results will be written.
        replicates: Number of replicates per mutation.
        jitter_seed_base: Base seed for coordinate jitter.
        jitter_angstrom: Coordinate jitter magnitude.

    Returns:
        Tuple of (structural_metrics_df, fep_tasks).
    """
    spec = dor_4ncg_spec(root)
    dor_df = load_dor_susceptibilities(susceptibility_xlsx, default_chain="A")
    chain_map = load_chain_subunits(spec.structure.cif_path)
    residue_maps = load_residue_mappings(spec.structure.cif_path)
    numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)

    ensure_dirs([prepared_dir, fep_manifest_path.parent, structural_metrics_path.parent])

    wt_cif = prepared_dir / "wt_4ncg.cif"
    if not wt_cif.exists():
        shutil.copy2(spec.structure.cif_path, wt_cif)

    structures: list[dict] = [
        {
            "mutation": "WT",
            "safe_label": "wt",
            "cif_path": wt_cif,
            "fold_reduction": None,
        }
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
        structures.append({
            "mutation": mutation,
            "safe_label": safe_label,
            "cif_path": mut_cif,
            "fold_reduction": float(row["dor_fold_reduction"]),
        })

    logging.info("Prepared %d structures (1 WT + %d mutants)", len(structures), len(structures) - 1)

    metrics_rows = []
    fep_tasks = []
    task_id = 0

    for struct in structures:
        for replicate in range(1, replicates + 1):
            safe_label = struct["safe_label"]
            mut_dir = prepared_dir / safe_label / f"rep_{replicate:02d}"
            ensure_dirs([mut_dir])

            min_pdb = mut_dir / f"{safe_label}_minimized_rep{replicate:02d}.pdb"

            seed = None
            if jitter_seed_base is not None:
                seed = jitter_seed_base + hash((safe_label, replicate)) % 100000

            logging.info(
                "Minimizing %s replicate %d", struct["mutation"], replicate
            )
            contacts, pocket_vol = prepare_structure_local(
                cif_path=struct["cif_path"],
                ligand_resname=spec.structure.ligand_resname,
                ligand_sdf=spec.structure.ligand_sdf,
                restraint_radius=spec.restraint_radius_angstrom,
                restraint_k=spec.restraint_k_kj_mol_nm2,
                output_path=min_pdb,
                jitter_seed=seed,
                jitter_angstrom=jitter_angstrom,
            )

            metrics_rows.append({
                "structure": "DOR",
                "mutation": struct["mutation"],
                "safe_label": safe_label,
                "replicate": replicate,
                "contact_count": contacts.contact_count,
                "hbond_count": contacts.hbond_count,
                "pocket_volume_proxy": pocket_vol,
            })

            for leg in ["complex", "solvent"]:
                output_json = (
                    fep_results_dir / safe_label / f"rep_{replicate:02d}" /
                    f"{safe_label}_{leg}_rep{replicate:02d}.json"
                )
                fep_tasks.append(FEPTask(
                    task_id=task_id,
                    structure="DOR",
                    mutation=struct["mutation"],
                    safe_label=safe_label,
                    replicate=replicate,
                    leg=leg,
                    minimized_pdb=str(min_pdb.resolve()),
                    ligand_sdf=str(spec.structure.ligand_sdf.resolve()),
                    ligand_resname=spec.structure.ligand_resname,
                    fold_reduction=struct["fold_reduction"],
                    output_json=str(output_json),
                ))
                task_id += 1

    structural_metrics = pd.DataFrame(metrics_rows)
    structural_metrics.to_csv(structural_metrics_path, index=False)
    logging.info("Wrote structural metrics to %s", structural_metrics_path)

    save_manifest(fep_tasks, fep_manifest_path)
    logging.info(
        "Wrote FEP manifest with %d tasks to %s", len(fep_tasks), fep_manifest_path
    )

    return structural_metrics, fep_tasks


def run_dor_alchemical_manifest(
    root: Path,
    manifest_csv: Path,
    output_dir: Path,
    replicates: int,
    jitter_seed_base: int | None,
    jitter_angstrom: float,
    alchemy_config: AlchemicalConfig,
    mutation_index: int | None = None,
) -> pd.DataFrame:
    """Run DOR alchemical workflow from a prepared manifest (legacy mode)."""
    from .structure_prep import compute_alchemical_binding_metric, prepare_structure

    spec = dor_4ncg_spec(root)
    manifest = pd.read_csv(manifest_csv)
    if mutation_index is not None:
        manifest = manifest.iloc[[mutation_index]]

    ensure_dirs([output_dir])
    records = []
    for replicate in range(1, replicates + 1):
        wt_dir = output_dir / "wt" / f"rep_{replicate:02d}"
        ensure_dirs([wt_dir])
        wt_min_pdb = wt_dir / f"wt_minimized_rep{replicate:02d}.pdb"
        wt_seed = None
        if jitter_seed_base is not None:
            wt_seed = jitter_seed_base + replicate * 100000 + 1
        _, wt_contacts, wt_pocket = prepare_structure(
            cif_path=Path(manifest.iloc[0]["wt_cif"]),
            ligand_resname=spec.structure.ligand_resname,
            ligand_sdf=spec.structure.ligand_sdf,
            restraint_radius=spec.restraint_radius_angstrom,
            restraint_k=spec.restraint_k_kj_mol_nm2,
            output_path=wt_min_pdb,
            jitter_seed=wt_seed,
            jitter_angstrom=jitter_angstrom,
        )
        wt_binding = compute_alchemical_binding_metric(
            minimized_pdb_path=wt_min_pdb,
            ligand_resname=spec.structure.ligand_resname,
            ligand_sdf=spec.structure.ligand_sdf,
            config=alchemy_config,
            output_json=wt_dir / f"wt_alchemy_rep{replicate:02d}.json",
            metadata={"structure": "DOR", "state": "WT", "replicate": replicate},
        )
        wt_metrics = {
            "binding_delta_g_kj_mol": wt_binding,
            "binding_proxy_kj_mol": float("nan"),
            "contact_count": wt_contacts.contact_count,
            "hbond_count": wt_contacts.hbond_count,
            "pocket_volume_proxy": wt_pocket,
        }

        mut_results = []
        for _, row in manifest.iterrows():
            mutation = str(row["mutation"])
            safe_label = str(row["safe_label"])
            mut_dir = output_dir / safe_label / f"rep_{replicate:02d}"
            ensure_dirs([mut_dir])
            mut_min_pdb = mut_dir / f"mut_minimized_{safe_label}_rep{replicate:02d}.pdb"
            _, mut_contacts, mut_pocket = prepare_structure(
                cif_path=Path(str(row["prepared_cif"])),
                ligand_resname=spec.structure.ligand_resname,
                ligand_sdf=spec.structure.ligand_sdf,
                restraint_radius=spec.restraint_radius_angstrom,
                restraint_k=spec.restraint_k_kj_mol_nm2,
                output_path=mut_min_pdb,
                jitter_seed=None,
                jitter_angstrom=jitter_angstrom,
            )
            mut_binding = compute_alchemical_binding_metric(
                minimized_pdb_path=mut_min_pdb,
                ligand_resname=spec.structure.ligand_resname,
                ligand_sdf=spec.structure.ligand_sdf,
                config=alchemy_config,
                output_json=mut_dir / f"mut_alchemy_{safe_label}_rep{replicate:02d}.json",
                metadata={
                    "structure": "DOR",
                    "state": "MUT",
                    "mutation": mutation,
                    "replicate": replicate,
                },
            )
            mut_results.append(
                {
                    "base": {
                        "structure": "DOR",
                        "mutation": mutation,
                        "mutation_order": int(row["mutation_order"]),
                        "category": None,
                        "notes": None,
                        "chain": str(row["chain"]),
                        "subunit": "p66",
                        "replicate": replicate,
                        "dor_fold_reduction": float(row["dor_fold_reduction"]),
                    },
                    "mut_metrics": {
                        "binding_delta_g_kj_mol": mut_binding,
                        "binding_proxy_kj_mol": float("nan"),
                        "contact_count": mut_contacts.contact_count,
                        "hbond_count": mut_contacts.hbond_count,
                        "pocket_volume_proxy": mut_pocket,
                    },
                }
            )
        records.extend(build_rows(mut_results, wt_metrics))
    df = pd.DataFrame(records)
    return df


def summarize_dor_correlations(metrics_df: pd.DataFrame) -> pd.DataFrame:
    """Compute correlations between metrics and fold reduction."""
    muts = (
        metrics_df[metrics_df["state"] == "MUT"][["mutation", "dor_fold_reduction"]]
        .drop_duplicates()
    )
    out = []
    for metric in ["binding_delta_g_kj_mol", "contact_count", "hbond_count", "pocket_volume_proxy"]:
        m = metrics_df[metrics_df["metric"] == metric]
        piv = m.pivot_table(
            index=["mutation", "replicate", "dor_fold_reduction"],
            columns="state",
            values="value",
            aggfunc="first",
        ).reset_index()
        piv = piv.dropna(subset=["WT", "MUT"])
        piv["delta"] = piv["MUT"] - piv["WT"]
        by_mut = piv.groupby("mutation", as_index=False).agg(
            delta_mean=("delta", "mean"),
            dor_fold_reduction=("dor_fold_reduction", "first"),
        )
        merged = by_mut.merge(muts, on=["mutation", "dor_fold_reduction"], how="left")
        if merged.empty:
            continue
        pearson = merged["delta_mean"].corr(merged["dor_fold_reduction"], method="pearson")
        spearman = merged["delta_mean"].corr(
            merged["dor_fold_reduction"], method="spearman"
        )
        out.append(
            {
                "metric": metric,
                "pearson_r": float(pearson) if not np.isnan(pearson) else float("nan"),
                "spearman_rho": float(spearman)
                if not np.isnan(spearman)
                else float("nan"),
                "n_mutations": int(len(merged)),
            }
        )
    return pd.DataFrame(out)
