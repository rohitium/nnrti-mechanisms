#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .plot_binding_energy_summary import _place_greedy_annotations
from .plot_dor_susceptibility_bars import CATEGORY_COLORS, _category_for_mutation
from .plot_triplet_contact_story import _load_replicate_meta
from .run_custom_mechanism_panel_models import _compute_custom_replicate_means


FEATURE_SPECS = [
    ("ser105_dor_distance_angstrom", "SER105-DOR Distance", "SER105-DOR Distance (Å)"),
    ("residue_min_distance_TYR188_angstrom", "TYR188-DOR Distance", "TYR188-DOR Distance (Å)"),
    ("ligand_pose_rmsd_angstrom", "Ligand Pose RMSD", "Ligand Pose RMSD (Å)"),
]

DISPLAY_TEST_SET_LABEL = "Test set"
DISPLAY_CATEGORY_COLORS = {
    "Negative control": CATEGORY_COLORS["Negative control"],
    "Positive control": CATEGORY_COLORS["Positive control"],
    DISPLAY_TEST_SET_LABEL: CATEGORY_COLORS["Uncertain/limited data"],
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot selected mechanism features against DOR fold-change.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument(
        "--frame-features-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--mechanism-panel-csv",
        type=Path,
        default=Path("results/analysis/custom_mechanism_selected_model/tables/mechanism_panel_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/custom_mechanism_selected_model"),
    )
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=5)
    return parser.parse_args()


def _sem(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size <= 1:
        return 0.0
    return float(np.nanstd(arr, ddof=1) / np.sqrt(arr.size))


def _display_category_for_control_category(control_category: str, mutation: str) -> str:
    label = str(mutation).strip().upper()
    category = str(control_category).strip().lower()
    if label == "WT" or category == "wt_reference":
        return "WT"
    if category == "negative_control":
        return "Negative control"
    if category == "positive_control":
        return "Positive control"
    if category == "uncertain_limited":
        return DISPLAY_TEST_SET_LABEL
    return _category_for_mutation(str(mutation))


def _build_feature_summary(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel = pd.read_csv(args.mechanism_panel_csv)
    panel_category_map = panel.set_index(panel["mutation"].astype(str))["control_category"].astype(str).to_dict()
    mutations = set(panel["mutation"].astype(str).tolist())
    metas = _load_replicate_meta(args.manifest, needed_mutations=mutations)

    rep_custom = _compute_custom_replicate_means(
        metas=metas,
        ligand_resname=str(args.ligand_resname),
        resid_offset=int(args.resid_offset),
        frame_stride=int(args.frame_stride),
    )[["mutation", "replicate", "ser105_dor_distance_angstrom_mean"]].rename(
        columns={"ser105_dor_distance_angstrom_mean": "ser105_dor_distance_angstrom"}
    )

    frame_df = pd.read_csv(args.frame_features_csv, usecols=["mutation", "replicate", "fold_reduction", "residue_min_distance_TYR188_angstrom", "ligand_pose_rmsd_angstrom"])
    rep_frame = (
        frame_df.groupby(["mutation", "replicate"], as_index=False)
        .agg(
            fold_reduction=("fold_reduction", "first"),
            residue_min_distance_TYR188_angstrom=("residue_min_distance_TYR188_angstrom", "mean"),
            ligand_pose_rmsd_angstrom=("ligand_pose_rmsd_angstrom", "mean"),
        )
    )

    rep = rep_frame.merge(rep_custom, on=["mutation", "replicate"], how="left")
    rep = rep[rep["mutation"].astype(str).isin(mutations)].copy()
    rep["control_category"] = rep["mutation"].astype(str).map(panel_category_map)
    rep["category"] = [
        _display_category_for_control_category(control_category=row.get("control_category", ""), mutation=row.get("mutation", ""))
        for _, row in rep.iterrows()
    ]

    agg_map: dict[str, tuple[str, str]] = {
        "fold_reduction": ("fold_reduction", "first"),
        "category": ("category", "first"),
        "n_replicates": ("replicate", "nunique"),
    }
    for column, _label, _ylabel in FEATURE_SPECS:
        agg_map[f"{column}_mean"] = (column, "mean")
        agg_map[f"{column}_sem"] = (column, _sem)

    summary = rep.groupby("mutation", as_index=False).agg(**agg_map).reset_index(drop=True)
    return rep, summary


def _plot_feature(summary: pd.DataFrame, column: str, label: str, ylabel: str, output_png: Path) -> dict[str, float]:
    import matplotlib.pyplot as plt

    df = summary[summary["mutation"].astype(str).str.upper() != "WT"].copy()
    df["display_category"] = df["category"].astype(str)
    df["log10_fold_reduction"] = np.log10(pd.to_numeric(df["fold_reduction"], errors="coerce"))
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_fold_reduction"]).reset_index(drop=True)
    if df.empty:
        return {"r_squared": np.nan, "p_value": np.nan}

    fig, ax = plt.subplots(figsize=(14.2, 9.2))
    for category, point_color in [
        ("Negative control", DISPLAY_CATEGORY_COLORS["Negative control"]),
        ("Positive control", DISPLAY_CATEGORY_COLORS["Positive control"]),
        (DISPLAY_TEST_SET_LABEL, DISPLAY_CATEGORY_COLORS[DISPLAY_TEST_SET_LABEL]),
    ]:
        subset = df[df["display_category"] == category].copy()
        if subset.empty:
            continue
        ax.errorbar(
            subset["fold_reduction"],
            subset[f"{column}_mean"],
            yerr=subset[f"{column}_sem"],
            fmt="o",
            ms=7.5,
            lw=0,
            elinewidth=1.15,
            capsize=2.5,
            color=point_color,
            markerfacecolor=point_color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            alpha=0.95,
            label=category,
            zorder=3,
        )

    valid = df[["log10_fold_reduction", f"{column}_mean"]].dropna()
    if len(valid) >= 3:
        slope, intercept, r_value, p_value, _stderr = stats.linregress(valid["log10_fold_reduction"], valid[f"{column}_mean"])
        x_grid = np.geomspace(float(df["fold_reduction"].min()) * 0.9, float(df["fold_reduction"].max()) * 1.1, 300)
        y_grid = slope * np.log10(x_grid) + intercept
        ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.5, zorder=2)
        annotation = f"R^2 = {r_value**2:.3f}\np = {p_value:.3f}"
    else:
        r_value = p_value = np.nan
        annotation = "R^2 = NA\np = NA"

    ax.set_xscale("log")
    ax.grid(alpha=0.25)
    ax.set_xlabel("Fold Reduction")
    ax.set_ylabel(ylabel)
    legend_loc = "lower left" if column == "ser105_dor_distance_angstrom" else "lower right"
    ax.legend(loc=legend_loc, fontsize=10, frameon=True, framealpha=0.9, facecolor="white", edgecolor="#cccccc")
    _place_greedy_annotations(
        ax,
        df["fold_reduction"].to_numpy(dtype=float),
        df[f"{column}_mean"].to_numpy(dtype=float),
        df["mutation"].astype(str).tolist(),
    )
    ax.text(
        0.02,
        0.98,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan, "p_value": float(p_value) if np.isfinite(p_value) else np.nan}


def main() -> int:
    args = _parse_args()
    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    rep, summary = _build_feature_summary(args)
    rep.to_csv(out_tables / "selected_mechanism_feature_replicate_means.csv", index=False)
    summary.to_csv(out_tables / "selected_mechanism_feature_mutation_summary.csv", index=False)

    stats_rows: list[dict[str, object]] = []
    for column, label, ylabel in FEATURE_SPECS:
        slug = column.replace("_angstrom", "").replace("_mean", "")
        stats_dict = _plot_feature(summary, column, label, ylabel, out_plots / f"{slug}_vs_fold_change.png")
        stats_rows.append({"feature": column, "label": label, **stats_dict})

    pd.DataFrame(stats_rows).to_csv(out_tables / "selected_mechanism_feature_vs_fold_stats.csv", index=False)
    (out_config / "selected_mechanism_features_vs_fold_run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "frame_features_csv": str(args.frame_features_csv),
                "mechanism_panel_csv": str(args.mechanism_panel_csv),
                "output_dir": str(args.output_dir),
                "feature_columns": [x[0] for x in FEATURE_SPECS],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
