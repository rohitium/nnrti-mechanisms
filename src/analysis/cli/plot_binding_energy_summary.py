#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

COMPONENT_SPECS = [
    ("binding_dg", "Total", "#d62828"),
    ("binding_dg_vdw", "vdW", "#1f77b4"),
    ("binding_dg_electrostatic", "Electrostatics", "#2a9d8f"),
    ("binding_dg_gb", "GB (Polar Solvation)", "#f4a261"),
    ("binding_dg_sa", "SA (Nonpolar)", "#6d597a"),
]

DDG_SPECS = [
    ("ddg", "Total ΔΔG", "#d62828"),
    ("ddg_vdw", "ΔΔG vdW", "#1f77b4"),
    ("ddg_electrostatic", "ΔΔG Electrostatics", "#2a9d8f"),
    ("ddg_gb", "ΔΔG GB", "#f4a261"),
    ("ddg_sa", "ΔΔG SA", "#6d597a"),
]


def _mutation_sort_key(mutation: str, fold_lookup: dict[str, float]) -> tuple[float, float, str]:
    if mutation == "WT":
        return (0.0, 0.0, mutation)
    if "+" in mutation:
        return (2.0, float(fold_lookup.get(mutation, np.inf)), mutation)
    return (1.0, float(fold_lookup.get(mutation, np.inf)), mutation)


def _sem(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
    if values.size <= 1:
        return 0.0
    return float(np.nanstd(values, ddof=1) / np.sqrt(values.size))


def _aggregate_summary(df: pd.DataFrame, specs: list[tuple[str, str, str]]) -> pd.DataFrame:
    base = (
        df.groupby("mutation", as_index=False)
        .agg(
            fold_reduction=("fold_reduction", "first"),
            n_replicates=("replicate", "nunique"),
        )
        .reset_index(drop=True)
    )
    for column, _label, _color in specs:
        grouped = df.groupby("mutation")[column]
        base[f"{column}_mean"] = grouped.mean().to_numpy(dtype=float)
        base[f"{column}_std"] = grouped.std(ddof=1).fillna(0.0).to_numpy(dtype=float)
        base[f"{column}_sem"] = grouped.apply(_sem).to_numpy(dtype=float)
    base["is_combo"] = base["mutation"].astype(str).str.contains(r"\+")
    return base


def _plot_by_mutation(summary: pd.DataFrame, specs: list[tuple[str, str, str]], title: str, output_png: Path) -> None:
    fold_lookup = {
        str(row["mutation"]): float(row["fold_reduction"])
        for _, row in summary.iterrows()
        if pd.notna(row["fold_reduction"])
    }
    ordered = summary.copy()
    ordered["sort_key"] = ordered["mutation"].astype(str).map(lambda x: _mutation_sort_key(x, fold_lookup))
    ordered = ordered.sort_values("sort_key").drop(columns=["sort_key"]).reset_index(drop=True)

    x = np.arange(len(ordered), dtype=float)
    fig, axes = plt.subplots(
        nrows=len(specs),
        ncols=1,
        figsize=(max(11.0, 0.75 * len(ordered)), 2.45 * len(specs) + 1.2),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    for ax, (column, label, color) in zip(axes, specs):
        y = ordered[f"{column}_mean"].to_numpy(dtype=float)
        yerr = ordered[f"{column}_sem"].to_numpy(dtype=float)
        ax.scatter(x, y, s=28, color=color, alpha=0.95, zorder=3)
        ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor=color, elinewidth=1.1, capsize=2.5, alpha=0.95, zorder=2)
        ax.grid(axis="y", linestyle=":", alpha=0.32)
        ax.set_ylabel(f"{label}\n(kJ/mol)", fontsize=9)
        if column in {"binding_dg", "ddg"}:
            ax.axhline(0.0, color="#999999", linestyle="--", linewidth=1.0, alpha=0.8)

    axes[-1].set_xticks(x, labels=ordered["mutation"].astype(str).tolist(), rotation=45, ha="right", fontsize=8)
    axes[-1].set_xlabel("Mutation", fontsize=10)
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.995)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _component_fold_stats(summary: pd.DataFrame, specs: list[tuple[str, str, str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    df = summary[summary["mutation"].astype(str) != "WT"].copy()
    df["log10_fold_reduction"] = np.log10(pd.to_numeric(df["fold_reduction"], errors="coerce"))
    for column, label, _color in specs:
        sub = df[["mutation", "fold_reduction", "log10_fold_reduction", f"{column}_mean"]].dropna().copy()
        if len(sub) >= 3:
            pearson_r, pearson_p = stats.pearsonr(sub["log10_fold_reduction"], sub[f"{column}_mean"])
            slope, intercept, r_value, p_value, _stderr = stats.linregress(sub["log10_fold_reduction"], sub[f"{column}_mean"])
        else:
            pearson_r = pearson_p = slope = intercept = r_value = p_value = np.nan
        rows.append(
            {
                "component": column,
                "label": label,
                "n_mutations": int(len(sub)),
                "pearson_r": float(pearson_r) if np.isfinite(pearson_r) else np.nan,
                "pearson_pvalue": float(pearson_p) if np.isfinite(pearson_p) else np.nan,
                "r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan,
                "slope_per_log10_fold": float(slope) if np.isfinite(slope) else np.nan,
                "intercept": float(intercept) if np.isfinite(intercept) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _plot_vs_fold_change(summary: pd.DataFrame, specs: list[tuple[str, str, str]], title: str, output_png: Path) -> None:
    df = summary[summary["mutation"].astype(str) != "WT"].copy()
    df["log10_fold_reduction"] = np.log10(pd.to_numeric(df["fold_reduction"], errors="coerce"))
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_fold_reduction"]).reset_index(drop=True)
    if df.empty:
        return

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.6), constrained_layout=True)
    axes_flat = axes.ravel()
    legend_handles = None
    legend_labels = None

    for ax, (column, label, _color) in zip(axes_flat, specs):
        singles = df[~df["is_combo"]].copy()
        combos = df[df["is_combo"]].copy()

        plot_args = [
            (singles, "#1d3557", "Single DRM", "o"),
            (combos, "#d62828", "Combination DRM", "s"),
        ]
        for subset, color, legend_label, marker in plot_args:
            if subset.empty:
                continue
            handle = ax.errorbar(
                subset["fold_reduction"],
                subset[f"{column}_mean"],
                yerr=subset[f"{column}_sem"],
                fmt=marker,
                ms=6,
                lw=0,
                elinewidth=1.0,
                capsize=2.0,
                color=color,
                alpha=0.9,
                label=legend_label,
            )
            if legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()

        valid = df[["log10_fold_reduction", f"{column}_mean"]].dropna()
        if len(valid) >= 3:
            slope, intercept, r_value, p_value, _stderr = stats.linregress(
                valid["log10_fold_reduction"],
                valid[f"{column}_mean"],
            )
            x_grid = np.geomspace(float(df["fold_reduction"].min()) * 0.9, float(df["fold_reduction"].max()) * 1.1, 200)
            y_grid = slope * np.log10(x_grid) + intercept
            ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.4, zorder=1)
            annotation = f"R^2 = {r_value**2:.3f}\np = {p_value:.3f}"
        else:
            annotation = "R^2 = NA\np = NA"

        ax.set_xscale("log")
        ax.grid(alpha=0.25)
        ax.set_title(label, fontsize=11, fontweight="bold")
        ax.set_xlabel("Fold Reduction")
        ax.set_ylabel("kJ/mol")
        ax.text(
            0.02,
            0.98,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9},
        )

    if len(specs) < len(axes_flat):
        for ax in axes_flat[len(specs) :]:
            ax.axis("off")

    if legend_handles and legend_labels:
        axes_flat[0].legend(legend_handles, legend_labels, loc="best", fontsize=9, frameon=False)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create consolidated MM/GBSA binding-energy summaries.")
    parser.add_argument("--replicate-csv", type=Path, default=Path("results/mmgbsa_replicate_metrics.csv"))
    parser.add_argument("--ddg-csv", type=Path, default=Path("results/ddg_full.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/binding_energy"))
    args = parser.parse_args()

    if not args.replicate_csv.exists():
        raise FileNotFoundError(args.replicate_csv)
    if not args.ddg_csv.exists():
        raise FileNotFoundError(args.ddg_csv)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    rep_df = pd.read_csv(args.replicate_csv)
    ddg_df = pd.read_csv(args.ddg_csv)

    rep_summary = _aggregate_summary(rep_df, COMPONENT_SPECS)
    ddg_summary = _aggregate_summary(ddg_df[ddg_df["mutation"].astype(str) != "WT"].copy(), DDG_SPECS)

    rep_df.to_csv(out_tables / "mmgbsa_replicate_metrics.csv", index=False)
    ddg_df.to_csv(out_tables / "ddg_full.csv", index=False)
    rep_summary.to_csv(out_tables / "mutation_component_summary.csv", index=False)
    ddg_summary.to_csv(out_tables / "mutation_ddg_summary.csv", index=False)
    _component_fold_stats(rep_summary, COMPONENT_SPECS).to_csv(out_tables / "component_vs_fold_change_stats.csv", index=False)
    _component_fold_stats(ddg_summary, DDG_SPECS).to_csv(out_tables / "ddg_vs_fold_change_stats.csv", index=False)

    _plot_by_mutation(rep_summary, COMPONENT_SPECS, "MM/GBSA Components By Mutation", out_plots / "mmgbsa_components_by_mutation.png")
    _plot_by_mutation(ddg_summary, DDG_SPECS, "MM/GBSA Components Relative To WT", out_plots / "mmgbsa_ddg_components_vs_wt.png")
    _plot_vs_fold_change(rep_summary, COMPONENT_SPECS, "MM/GBSA Components Vs Fold Reduction", out_plots / "mmgbsa_components_vs_fold_change.png")
    _plot_vs_fold_change(ddg_summary, DDG_SPECS, "MM/GBSA ΔΔG Components Vs Fold Reduction", out_plots / "mmgbsa_ddg_vs_fold_change.png")

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "replicate_csv": str(args.replicate_csv),
                "ddg_csv": str(args.ddg_csv),
                "component_columns": [x[0] for x in COMPONENT_SPECS],
                "ddg_columns": [x[0] for x in DDG_SPECS],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
