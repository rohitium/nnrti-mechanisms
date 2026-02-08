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
    """Collect all FEP leg results from JSON files."""
    del fep_results_dir
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
        except json.JSONDecodeError as exc:
            logging.error("Invalid JSON for task %d: %s", task.task_id, exc)
            continue

        rows.append(
            {
                "task_id": task.task_id,
                "structure": task.structure,
                "mutation": task.mutation,
                "safe_label": task.safe_label,
                "replicate": task.replicate,
                "leg": task.leg,
                "delta_g_kj_mol": data.get("delta_g_kj_mol"),
                "fold_reduction": task.fold_reduction,
                "minimized_pdb": task.minimized_pdb or data.get("minimized_pdb", ""),
                "topology_pdb": data.get("topology_pdb", ""),
                "trajectory_dcd": data.get("trajectory_dcd", ""),
                "physical_trajectory_dcd": data.get("physical_trajectory_dcd", ""),
            }
        )

    return pd.DataFrame(rows)


def collect_lambda_profiles(manifest_path: Path) -> pd.DataFrame:
    """Collect per-window free-energy profile data for protocol sanity checks."""
    tasks = load_manifest(manifest_path)
    rows: list[dict] = []

    for task in tasks:
        json_path = Path(task.output_json)
        if not json_path.exists():
            continue

        try:
            with open(json_path) as f:
                data = json.load(f)
        except Exception:
            continue

        pair = data.get("pair_delta_g_kj_mol") or []
        protocol = data.get("lambda_protocol") or []
        if not pair or len(protocol) != len(pair) + 1:
            continue

        cumulative = 0.0
        for i, pair_dg in enumerate(pair):
            lam0_e, lam0_s = protocol[i]
            lam1_e, lam1_s = protocol[i + 1]
            pair_dg = float(pair_dg)
            cumulative += pair_dg
            rows.append(
                {
                    "structure": task.structure,
                    "mutation": task.mutation,
                    "safe_label": task.safe_label,
                    "replicate": task.replicate,
                    "leg": task.leg,
                    "window_index": i,
                    "lambda0_electrostatics": float(lam0_e),
                    "lambda0_sterics": float(lam0_s),
                    "lambda1_electrostatics": float(lam1_e),
                    "lambda1_sterics": float(lam1_s),
                    "pair_delta_g_kj_mol": pair_dg,
                    "cumulative_delta_g_kj_mol": cumulative,
                    "fold_reduction": task.fold_reduction,
                }
            )

    return pd.DataFrame(rows)


def summarize_lambda_profiles(lambda_df: pd.DataFrame) -> pd.DataFrame:
    if lambda_df.empty:
        return pd.DataFrame()
    return (
        lambda_df.groupby(["mutation", "leg", "window_index"], as_index=False)
        .agg(
            pair_delta_g_mean=("pair_delta_g_kj_mol", "mean"),
            pair_delta_g_std=("pair_delta_g_kj_mol", "std"),
            cumulative_delta_g_mean=("cumulative_delta_g_kj_mol", "mean"),
            cumulative_delta_g_std=("cumulative_delta_g_kj_mol", "std"),
            n_replicates=("replicate", "nunique"),
        )
        .sort_values(["mutation", "leg", "window_index"])
        .reset_index(drop=True)
    )


def compute_structural_metrics(
    fep_df: pd.DataFrame,
    ligand_resname: str,
    frame_stride: int = 5,
    max_frames: int = 200,
) -> pd.DataFrame:
    """Compute ensemble-averaged structural metrics from complex-leg trajectories."""
    from ..analysis_metrics import (
        compute_contacts,
        compute_ensemble_metrics,
        pocket_volume_proxy,
    )

    complex_rows = fep_df[fep_df["leg"] == "complex"].copy()
    if complex_rows.empty:
        return pd.DataFrame()

    rows = []
    for _, row in complex_rows.iterrows():
        mutation = row["mutation"]
        replicate = int(row["replicate"])
        topo_path = Path(str(row.get("topology_pdb") or row.get("minimized_pdb") or ""))
        physical_dcd = Path(str(row.get("physical_trajectory_dcd") or ""))

        if not topo_path.exists():
            logging.warning(
                "Missing topology for %s rep%d: %s",
                mutation,
                replicate,
                topo_path,
            )
            continue

        try:
            if physical_dcd.exists():
                ens = compute_ensemble_metrics(
                    topology_pdb_path=topo_path,
                    trajectory_dcd_path=physical_dcd,
                    ligand_resname=ligand_resname,
                    frame_stride=frame_stride,
                    max_frames=max_frames,
                )
                rows.append(
                    {
                        "structure": row["structure"],
                        "mutation": mutation,
                        "safe_label": row["safe_label"],
                        "replicate": replicate,
                        "contact_count": ens.contact_count_mean,
                        "contact_count_std": ens.contact_count_std,
                        "hbond_count": ens.hbond_count_mean,
                        "hbond_count_std": ens.hbond_count_std,
                        "pocket_volume_proxy": ens.pocket_volume_proxy_mean,
                        "pocket_volume_proxy_std": ens.pocket_volume_proxy_std,
                        "metric_n_frames": ens.n_frames,
                        "metric_source": "trajectory",
                        "fold_reduction": row["fold_reduction"],
                    }
                )
            else:
                contacts = compute_contacts(topo_path, ligand_resname=ligand_resname)
                pocket = pocket_volume_proxy(topo_path, ligand_resname=ligand_resname)
                rows.append(
                    {
                        "structure": row["structure"],
                        "mutation": mutation,
                        "safe_label": row["safe_label"],
                        "replicate": replicate,
                        "contact_count": float(contacts.contact_count),
                        "contact_count_std": np.nan,
                        "hbond_count": float(contacts.hbond_count or 0.0),
                        "hbond_count_std": np.nan,
                        "pocket_volume_proxy": float(pocket),
                        "pocket_volume_proxy_std": np.nan,
                        "metric_n_frames": 1,
                        "metric_source": "single_structure",
                        "fold_reduction": row["fold_reduction"],
                    }
                )
        except Exception as exc:
            logging.error(
                "Failed structural metrics for %s rep%d: %s",
                mutation,
                replicate,
                exc,
            )

    return pd.DataFrame(rows)


def compute_binding_ddg(fep_df: pd.DataFrame) -> pd.DataFrame:
    """Compute binding ΔG and ΔΔG from collected FEP leg results."""
    fold_map = fep_df.drop_duplicates(
        subset=["structure", "mutation", "safe_label", "replicate"]
    )[["structure", "mutation", "safe_label", "replicate", "fold_reduction"]]

    pivot = fep_df.pivot_table(
        index=["structure", "mutation", "safe_label", "replicate"],
        columns="leg",
        values="delta_g_kj_mol",
        aggfunc="first",
    ).reset_index()
    pivot = pivot.merge(
        fold_map,
        on=["structure", "mutation", "safe_label", "replicate"],
        how="left",
    )

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
        pivot["wt_binding_dg"] = pivot.apply(
            lambda row: wt_by_rep.get((row["structure"], row["replicate"]), np.nan),
            axis=1,
        )
        pivot["ddg"] = pivot["binding_dg"] - pivot["wt_binding_dg"]

    result = pivot.rename(columns={"complex": "complex_dg", "solvent": "solvent_dg"})
    return result[
        [
            "structure",
            "mutation",
            "safe_label",
            "replicate",
            "complex_dg",
            "solvent_dg",
            "binding_dg",
            "wt_binding_dg",
            "ddg",
            "fold_reduction",
        ]
    ]


def merge_with_structural_metrics(
    ddg_df: pd.DataFrame,
    structural_metrics_df: pd.DataFrame,
) -> pd.DataFrame:
    if structural_metrics_df.empty:
        return ddg_df

    merged = ddg_df.merge(
        structural_metrics_df,
        on=["structure", "mutation", "safe_label", "replicate"],
        how="left",
        suffixes=("", "_struct"),
    )
    if "fold_reduction_struct" in merged.columns:
        merged = merged.drop(columns=["fold_reduction_struct"])
    return merged


def _corr_rows(metric_name: str, x: pd.Series, y: pd.Series) -> list[dict]:
    mask = np.isfinite(x.values) & np.isfinite(y.values)
    if mask.sum() < 3:
        return []
    x_valid = x[mask]
    y_valid = y[mask]
    return [
        {
            "metric": metric_name,
            "pearson_r": float(x_valid.corr(y_valid, method="pearson")),
            "spearman_rho": float(x_valid.corr(y_valid, method="spearman")),
            "n_mutations": int(mask.sum()),
        },
        {
            "metric": f"{metric_name}_log10_fold",
            "pearson_r": float(x_valid.corr(np.log10(y_valid + 1), method="pearson")),
            "spearman_rho": np.nan,
            "n_mutations": int(mask.sum()),
        },
    ]


def compute_correlations(ddg_df: pd.DataFrame) -> pd.DataFrame:
    """Compute metric-vs-susceptibility correlations."""
    mut_df = ddg_df[ddg_df["mutation"] != "WT"].copy()
    if mut_df.empty:
        logging.warning("No mutation rows available for correlation analysis")
        return pd.DataFrame()

    wt_df = ddg_df[ddg_df["mutation"] == "WT"].set_index(["structure", "replicate"])
    metrics = ["ddg", "contact_count", "hbond_count", "pocket_volume_proxy"]

    for metric in metrics:
        if metric not in mut_df.columns:
            continue
        if metric in wt_df.columns:
            wt_lookup = wt_df[metric]
            mut_df[f"{metric}_delta"] = mut_df.apply(
                lambda row: row[metric]
                - wt_lookup.get((row["structure"], row["replicate"]), np.nan),
                axis=1,
            )

    agg_cols: dict[str, tuple[str, str]] = {
        "fold_reduction": ("fold_reduction", "first"),
        "n_replicates": ("replicate", "nunique"),
    }
    for metric in metrics:
        if metric in mut_df.columns:
            agg_cols[f"{metric}_mean"] = (metric, "mean")
        delta_col = f"{metric}_delta"
        if delta_col in mut_df.columns:
            agg_cols[f"{delta_col}_mean"] = (delta_col, "mean")

    by_mut = mut_df.groupby("mutation", as_index=False).agg(**agg_cols)
    if by_mut.empty:
        logging.warning("No valid mutation aggregates for correlation analysis")
        return pd.DataFrame()

    results: list[dict] = []
    y = by_mut["fold_reduction"]
    for col in by_mut.columns:
        if col.endswith("_mean") and col not in {"fold_reduction", "n_replicates"}:
            results.extend(_corr_rows(col.replace("_mean", ""), by_mut[col], y))

    return pd.DataFrame(results)


def summarize_ddg_by_mutation(ddg_df: pd.DataFrame) -> pd.DataFrame:
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
    """Run the full result collection and analysis pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Collecting FEP results from %s", fep_results_dir)
    fep_df = collect_fep_results(manifest_path, fep_results_dir)
    if fep_df.empty:
        raise ValueError("No FEP results found")

    logging.info("Collecting lambda-window diagnostics")
    lambda_df = collect_lambda_profiles(manifest_path)
    if not lambda_df.empty:
        lambda_df.to_csv(output_dir / "lambda_window_profiles.csv", index=False)
        summarize_lambda_profiles(lambda_df).to_csv(
            output_dir / "lambda_window_summary.csv",
            index=False,
        )

    logging.info("Computing binding ΔG and ΔΔG")
    ddg_df = compute_binding_ddg(fep_df)

    logging.info("Computing ensemble structural metrics")
    struct_df = compute_structural_metrics(fep_df, ligand_resname)
    if not struct_df.empty:
        struct_df.to_csv(output_dir / "structural_metrics.csv", index=False)
        ddg_df = merge_with_structural_metrics(ddg_df, struct_df)

    ddg_df.to_csv(output_dir / "ddg_full.csv", index=False)

    logging.info("Summarizing ΔΔG by mutation")
    ddg_summary = summarize_ddg_by_mutation(ddg_df)
    ddg_summary.to_csv(output_dir / "ddg_summary.csv", index=False)

    logging.info("Computing susceptibility correlations")
    correlations = compute_correlations(ddg_df)
    correlations.to_csv(output_dir / "correlation_analysis.csv", index=False)

    return ddg_summary, correlations, ddg_df
