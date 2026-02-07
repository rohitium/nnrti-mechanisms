from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .manifest import load_manifest


def collect_fep_results(
    manifest_path: Path,
    fep_results_dir: Path,
) -> pd.DataFrame:
    """Collect all FEP leg results from JSON files.

    Args:
        manifest_path: Path to the FEP manifest CSV.
        fep_results_dir: Directory containing FEP result JSON files.

    Returns:
        DataFrame with columns: task_id, structure, mutation, safe_label,
        replicate, leg, delta_g_kj_mol, fold_reduction, minimized_pdb.
    """
    tasks = load_manifest(manifest_path)
    rows = []

    for task in tasks:
        json_path = Path(task.output_json)
        if not json_path.exists():
            logging.warning("Missing result for task %d: %s", task.task_id, json_path)
            continue

        try:
            with open(json_path) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            logging.error("Invalid JSON for task %d: %s", task.task_id, e)
            continue

        rows.append({
            "task_id": task.task_id,
            "structure": task.structure,
            "mutation": task.mutation,
            "safe_label": task.safe_label,
            "replicate": task.replicate,
            "leg": task.leg,
            "delta_g_kj_mol": data.get("delta_g_kj_mol"),
            "fold_reduction": task.fold_reduction,
            "minimized_pdb": task.minimized_pdb or data.get("minimized_pdb", ""),
        })

    return pd.DataFrame(rows)


def compute_structural_metrics(fep_df: pd.DataFrame, ligand_resname: str) -> pd.DataFrame:
    """Compute structural metrics from minimized PDB files.

    Args:
        fep_df: DataFrame from collect_fep_results().
        ligand_resname: Ligand residue name for contact calculations.

    Returns:
        DataFrame with structural metrics per (mutation, replicate).
    """
    from ..analysis_metrics import compute_contacts, pocket_volume_proxy

    # Get unique (mutation, replicate, minimized_pdb) combinations
    # (both legs share the same minimized structure)
    unique_structs = fep_df[["mutation", "safe_label", "replicate", "minimized_pdb", "fold_reduction"]].drop_duplicates()

    rows = []
    for _, row in unique_structs.iterrows():
        pdb_path = row["minimized_pdb"]
        if not pdb_path or not Path(pdb_path).exists():
            logging.warning(
                "Missing minimized PDB for %s rep%d: %s",
                row["mutation"], row["replicate"], pdb_path
            )
            continue

        try:
            contacts = compute_contacts(Path(pdb_path), ligand_resname=ligand_resname)
            pocket_vol = pocket_volume_proxy(Path(pdb_path), ligand_resname=ligand_resname)

            rows.append({
                "structure": "DOR",
                "mutation": row["mutation"],
                "safe_label": row["safe_label"],
                "replicate": row["replicate"],
                "contact_count": contacts.contact_count,
                "hbond_count": contacts.hbond_count,
                "pocket_volume_proxy": pocket_vol,
                "fold_reduction": row["fold_reduction"],
            })
        except Exception as e:
            logging.error(
                "Failed to compute metrics for %s rep%d: %s",
                row["mutation"], row["replicate"], e
            )

    return pd.DataFrame(rows)


def compute_binding_ddg(fep_df: pd.DataFrame) -> pd.DataFrame:
    """Compute binding ΔG and ΔΔG from collected FEP leg results.

    Args:
        fep_df: DataFrame from collect_fep_results().

    Returns:
        DataFrame with columns: structure, mutation, replicate, complex_dg,
        solvent_dg, binding_dg, wt_binding_dg, ddg, fold_reduction.
    """
    # Keep fold_reduction out of the pivot index to avoid dropping WT rows
    # (WT has fold_reduction=NaN, and pivot_table drops NaN index values).
    fold_map = fep_df.drop_duplicates(
        subset=["structure", "mutation", "safe_label", "replicate"]
    )[["structure", "mutation", "safe_label", "replicate", "fold_reduction"]]

    pivot = fep_df.pivot_table(
        index=["structure", "mutation", "safe_label", "replicate"],
        columns="leg",
        values="delta_g_kj_mol",
        aggfunc="first",
    ).reset_index()

    pivot = pivot.merge(fold_map, on=["structure", "mutation", "safe_label", "replicate"], how="left")

    if "complex" not in pivot.columns or "solvent" not in pivot.columns:
        logging.error("Missing complex or solvent leg in results")
        return pd.DataFrame()

    pivot["binding_dg"] = pivot["complex"] - pivot["solvent"]

    wt_rows = pivot[pivot["mutation"] == "WT"]
    if wt_rows.empty:
        logging.warning("No WT results found - cannot compute ΔΔG")
        pivot["wt_binding_dg"] = np.nan
        pivot["ddg"] = np.nan
    else:
        wt_by_rep = wt_rows.set_index(["structure", "replicate"])["binding_dg"]

        def get_wt_dg(row):
            key = (row["structure"], row["replicate"])
            return wt_by_rep.get(key, np.nan)

        pivot["wt_binding_dg"] = pivot.apply(get_wt_dg, axis=1)
        pivot["ddg"] = pivot["binding_dg"] - pivot["wt_binding_dg"]

    result = pivot.rename(columns={
        "complex": "complex_dg",
        "solvent": "solvent_dg",
    })

    return result[[
        "structure", "mutation", "safe_label", "replicate",
        "complex_dg", "solvent_dg", "binding_dg",
        "wt_binding_dg", "ddg", "fold_reduction",
    ]]


def merge_with_structural_metrics(
    ddg_df: pd.DataFrame,
    structural_metrics_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge ΔΔG results with structural metrics.

    Args:
        ddg_df: DataFrame from compute_binding_ddg().
        structural_metrics_df: DataFrame from compute_structural_metrics().

    Returns:
        Merged DataFrame with ΔΔG and structural metrics.
    """
    if structural_metrics_df.empty:
        return ddg_df

    merged = ddg_df.merge(
        structural_metrics_df,
        on=["structure", "mutation", "safe_label", "replicate"],
        how="left",
        suffixes=("", "_struct"),
    )

    # Handle duplicate fold_reduction column
    if "fold_reduction_struct" in merged.columns:
        merged = merged.drop(columns=["fold_reduction_struct"])

    return merged


def compute_correlations(ddg_df: pd.DataFrame) -> pd.DataFrame:
    """Compute correlations between ΔΔG and fold_reduction.

    Args:
        ddg_df: DataFrame from compute_binding_ddg().

    Returns:
        DataFrame with correlation statistics.
    """
    mut_df = ddg_df[ddg_df["mutation"] != "WT"].dropna(subset=["ddg", "fold_reduction"])

    if mut_df.empty:
        logging.warning("No valid mutation data for correlation analysis")
        return pd.DataFrame()

    by_mutation = mut_df.groupby("mutation", as_index=False).agg(
        ddg_mean=("ddg", "mean"),
        ddg_std=("ddg", "std"),
        fold_reduction=("fold_reduction", "first"),
        n_replicates=("replicate", "nunique"),
    )

    results = []

    pearson = by_mutation["ddg_mean"].corr(
        by_mutation["fold_reduction"], method="pearson"
    )
    spearman = by_mutation["ddg_mean"].corr(
        by_mutation["fold_reduction"], method="spearman"
    )

    results.append({
        "metric": "ddg_vs_fold_reduction",
        "pearson_r": float(pearson) if not np.isnan(pearson) else np.nan,
        "spearman_rho": float(spearman) if not np.isnan(spearman) else np.nan,
        "n_mutations": len(by_mutation),
    })

    log_corr = by_mutation["ddg_mean"].corr(
        np.log10(by_mutation["fold_reduction"] + 1), method="pearson"
    )
    results.append({
        "metric": "ddg_vs_log10_fold_reduction",
        "pearson_r": float(log_corr) if not np.isnan(log_corr) else np.nan,
        "spearman_rho": np.nan,
        "n_mutations": len(by_mutation),
    })

    return pd.DataFrame(results)


def summarize_ddg_by_mutation(ddg_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize ΔΔG across replicates for each mutation.

    Args:
        ddg_df: DataFrame from compute_binding_ddg().

    Returns:
        DataFrame with mean and std ΔΔG per mutation.
    """
    summary = ddg_df.groupby(
        ["structure", "mutation", "safe_label", "fold_reduction"], as_index=False
    ).agg(
        ddg_mean=("ddg", "mean"),
        ddg_std=("ddg", "std"),
        ddg_sem=("ddg", lambda x: x.std() / np.sqrt(len(x))),
        binding_dg_mean=("binding_dg", "mean"),
        n_replicates=("replicate", "nunique"),
    )

    return summary.sort_values("ddg_mean", ascending=False).reset_index(drop=True)


def run_result_collection(
    manifest_path: Path,
    fep_results_dir: Path,
    output_dir: Path,
    ligand_resname: str = "2KW",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the full result collection and analysis pipeline.

    Args:
        manifest_path: Path to FEP manifest CSV.
        fep_results_dir: Directory containing FEP result JSONs.
        output_dir: Directory for output files.
        ligand_resname: Ligand residue name for structural metrics.

    Returns:
        Tuple of (ddg_summary, correlation_analysis, full_ddg_df).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Collecting FEP results from %s", fep_results_dir)
    fep_df = collect_fep_results(manifest_path, fep_results_dir)

    if fep_df.empty:
        raise ValueError("No FEP results found")

    logging.info("Computing binding ΔG and ΔΔG")
    ddg_df = compute_binding_ddg(fep_df)

    logging.info("Computing structural metrics from minimized structures")
    struct_df = compute_structural_metrics(fep_df, ligand_resname)
    if not struct_df.empty:
        struct_df.to_csv(output_dir / "structural_metrics.csv", index=False)
        ddg_df = merge_with_structural_metrics(ddg_df, struct_df)

    ddg_df.to_csv(output_dir / "ddg_full.csv", index=False)

    logging.info("Summarizing ΔΔG by mutation")
    ddg_summary = summarize_ddg_by_mutation(ddg_df)
    ddg_summary.to_csv(output_dir / "ddg_summary.csv", index=False)

    logging.info("Computing correlations")
    correlations = compute_correlations(ddg_df)
    correlations.to_csv(output_dir / "correlation_analysis.csv", index=False)

    return ddg_summary, correlations, ddg_df
