#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .plot_binding_energy_summary import _place_greedy_annotations
from .plot_dor_susceptibility_bars import CATEGORY_COLORS


FEATURE_SPECS = [
    ("residue_min_distance_SER105_angstrom_mean", "SER105-DOR distance", "ser105_dor_distance"),
    ("residue_min_distance_TYR188_angstrom_mean", "TYR188-DOR distance", "tyr188_dor_distance"),
    ("ligand_pose_rmsd_angstrom_mean", "DOR pose RMSD", "dor_pose_rmsd"),
]

DISPLAY_CATEGORIES = {
    "negative_control": "Negative control",
    "positive_control": "Positive control",
    "uncertain_phenotype": "Uncertain Phenotype",
    "wt_reference": "WT",
}

AXIS_LABEL_SIZE = 22
TICK_LABEL_SIZE = 15
LEGEND_LABEL_SIZE = 16
POINT_LABEL_SIZE = 13
STATS_LABEL_SIZE = 20


def _sem(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size <= 1:
        return 0.0
    return float(np.nanstd(arr, ddof=1) / np.sqrt(arr.size))


def _summarize_features(rep_df: pd.DataFrame) -> pd.DataFrame:
    grouped = rep_df.groupby("mutation", as_index=False).agg(
        control_category=("control_category", "first"),
        target_fold_change=("target_fold_change", "first"),
        target_binary_class=("target_binary_class", "first"),
        n_replicates=("replicate", "nunique"),
    )
    for feature, _label, _slug in FEATURE_SPECS:
        if feature not in rep_df.columns:
            continue
        by_mut = rep_df.groupby("mutation")[feature]
        grouped[f"{feature}_mean"] = by_mut.mean().to_numpy(dtype=float)
        grouped[f"{feature}_sem"] = by_mut.apply(_sem).to_numpy(dtype=float)
    category_order = {"negative_control": 0, "wt_reference": 1, "uncertain_phenotype": 2, "positive_control": 3}
    grouped["_category_order"] = grouped["control_category"].map(category_order).fillna(9).astype(int)
    grouped = grouped.sort_values(["_category_order", "target_fold_change", "mutation"], kind="stable").drop(columns="_category_order")
    return grouped.reset_index(drop=True)


def _plot_feature(summary: pd.DataFrame, feature: str, label: str, output_png: Path) -> dict[str, object]:
    import matplotlib.pyplot as plt

    df = summary[summary["control_category"].isin(["negative_control", "positive_control", "wt_reference"])].copy()
    df = df[pd.to_numeric(df["target_fold_change"], errors="coerce").notna()].reset_index(drop=True)
    df["target_fold_change"] = df["target_fold_change"].astype(float)
    df["log10_fold_change"] = np.log10(df["target_fold_change"].clip(lower=1e-6))
    mean_col = f"{feature}_mean"
    sem_col = f"{feature}_sem"

    fig, ax = plt.subplots(figsize=(14.2, 9.2))
    category_order = ["wt_reference", "negative_control", "positive_control"]
    for category in category_order:
        subset = df[df["control_category"].astype(str) == category].copy()
        if subset.empty:
            continue
        display = DISPLAY_CATEGORIES.get(category, category)
        color = "#333333" if category == "wt_reference" else CATEGORY_COLORS.get(display, "#333333")
        ax.errorbar(
            subset["target_fold_change"],
            subset[mean_col],
            yerr=subset[sem_col],
            fmt="o",
            ms=7.5,
            lw=0,
            elinewidth=1.15,
            capsize=2.5,
            color=color,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            alpha=0.95,
            label=display,
            zorder=3,
        )

    valid = df[["log10_fold_change", mean_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) >= 3:
        slope, intercept, r_value, p_value, _stderr = stats.linregress(valid["log10_fold_change"], valid[mean_col])
        pearson_r, pearson_p = stats.pearsonr(valid["log10_fold_change"], valid[mean_col])
        x_grid = np.geomspace(float(df["target_fold_change"].min()) * 0.9, float(df["target_fold_change"].max()) * 1.1, 300)
        y_grid = slope * np.log10(x_grid) + intercept
        ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.5, zorder=2)
        annotation = f"R\u00b2 = {r_value**2:.3f}\np = {p_value:.3f}"
    else:
        slope = intercept = r_value = p_value = pearson_r = pearson_p = np.nan
        annotation = "R\u00b2 = NA\np = NA"

    fixed_offsets = {
        "V106I+F227C": (-8, 8, "right"),
        "V106A+F227L": (-8, 8, "right"),
        "A98G+F227C": (-8, 8, "right"),
        "Y188L": (-8, 8, "right"),
    }

    ax.set_xscale("log")
    ax.grid(alpha=0.25)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.set_xlabel("Fold-change", fontsize=AXIS_LABEL_SIZE, fontweight="bold")
    ax.set_ylabel(f"{label} (\u00c5)", fontsize=20, fontweight="bold")
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    ax.text(
        0.36,
        0.98,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=STATS_LABEL_SIZE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )
    _place_greedy_annotations(
        ax,
        df["target_fold_change"].to_numpy(dtype=float),
        df[mean_col].to_numpy(dtype=float),
        df["mutation"].astype(str).tolist(),
        fontsize=POINT_LABEL_SIZE,
        fixed_offsets={key: (xoff, yoff) for key, (xoff, yoff, _ha) in fixed_offsets.items()},
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=2.0)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

    return {
        "feature": feature,
        "label": label,
        "n_mutations": int(len(valid)),
        "pearson_r": float(pearson_r) if np.isfinite(pearson_r) else np.nan,
        "pearson_pvalue": float(pearson_p) if np.isfinite(pearson_p) else np.nan,
        "r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan,
        "slope_per_log10_fold": float(slope) if np.isfinite(slope) else np.nan,
        "intercept": float(intercept) if np.isfinite(intercept) else np.nan,
        "output_png": str(output_png),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot logistic-regression feature correlations with DOR fold-change.")
    parser.add_argument(
        "--replicate-feature-csv",
        type=Path,
        default=Path("results/analysis/new_logistic_regression/fixed_ser105_tyr188_pose_rmsd/tables/replicate_level_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/new_logistic_regression/fixed_ser105_tyr188_pose_rmsd"),
    )
    args = parser.parse_args()

    if not args.replicate_feature_csv.exists():
        raise FileNotFoundError(args.replicate_feature_csv)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    for directory in (out_tables, out_plots, out_config):
        directory.mkdir(parents=True, exist_ok=True)

    rep_df = pd.read_csv(args.replicate_feature_csv)
    summary = _summarize_features(rep_df)
    summary.to_csv(out_tables / "logistic_feature_fold_change_summary.csv", index=False)

    stats_rows = []
    for feature, label, slug in FEATURE_SPECS:
        if feature not in rep_df.columns:
            continue
        stats_rows.append(_plot_feature(summary, feature, label, out_plots / f"{slug}_vs_fold_change.png"))
    pd.DataFrame(stats_rows).to_csv(out_tables / "logistic_feature_vs_fold_change_stats.csv", index=False)

    (out_config / "logistic_feature_correlation_config.json").write_text(
        json.dumps(
            {
                "replicate_feature_csv": str(args.replicate_feature_csv),
                "output_dir": str(args.output_dir),
                "features": [feature for feature, _label, _slug in FEATURE_SPECS],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
