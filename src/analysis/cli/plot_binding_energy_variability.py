#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .plot_binding_energy_summary import (
    AXIS_LABEL_SIZE,
    CATEGORY_COLORS,
    LEGEND_LABEL_SIZE,
    POINT_LABEL_SIZE,
    STATS_LABEL_SIZE,
    TICK_LABEL_SIZE,
)
from .plot_dor_susceptibility_bars import _category_for_mutation


def _plot_ddg_vdw_std_vs_fold_change(summary: pd.DataFrame, output_png: Path) -> dict[str, object]:
    df = summary[summary["mutation"].astype(str) != "WT"].copy()
    df["category"] = df["mutation"].astype(str).map(_category_for_mutation)
    df = df[df["category"].isin(["Negative control", "Positive control"])].copy()
    df["fold_reduction"] = pd.to_numeric(df["fold_reduction"], errors="coerce")
    df["ddg_vdw_std"] = pd.to_numeric(df["ddg_vdw_std"], errors="coerce")
    df["log10_fold_reduction"] = np.log10(df["fold_reduction"])
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_fold_reduction", "ddg_vdw_std"])
    if df.empty:
        raise ValueError("No negative/positive control rows with finite fold-change and ddg_vdw_std values.")

    fig, ax = plt.subplots(figsize=(14.2, 9.2))
    categories = [
        ("Negative control", CATEGORY_COLORS["Negative control"]),
        ("Positive control", CATEGORY_COLORS["Positive control"]),
    ]
    for category, color in categories:
        subset = df[df["category"] == category]
        if subset.empty:
            continue
        ax.scatter(
            subset["fold_reduction"],
            subset["ddg_vdw_std"],
            s=72,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.95,
            label=category,
            zorder=3,
        )

    valid = df[["log10_fold_reduction", "ddg_vdw_std"]].dropna()
    if len(valid) >= 3:
        slope, intercept, r_value, p_value, stderr = stats.linregress(
            valid["log10_fold_reduction"],
            valid["ddg_vdw_std"],
        )
        x_grid = np.geomspace(float(df["fold_reduction"].min()) * 0.85, float(df["fold_reduction"].max()) * 1.6, 300)
        y_grid = slope * np.log10(x_grid) + intercept
        ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.5, zorder=2)
        annotation = f"R\u00b2 = {r_value**2:.3f}\np = {p_value:.3f}"
    else:
        slope = intercept = r_value = p_value = stderr = np.nan
        annotation = "R\u00b2 = NA\np = NA"

    label_offsets = {
        "A98G+F227C": (-10, 10),
        "V106I+F227C": (-10, 14),
        "V106A+F227L": (-10, -14),
        "V106A+P225H": (12, 14),
        "V106A+L234I": (12, 34),
    }
    for _, row in df.iterrows():
        mutation = str(row["mutation"])
        dx, dy = label_offsets.get(mutation, (8, 8))
        ax.annotate(
            mutation,
            xy=(float(row["fold_reduction"]), float(row["ddg_vdw_std"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=POINT_LABEL_SIZE,
            color="#333333",
            alpha=0.92,
            ha="right" if dx < 0 else "left",
            va="top" if dy < 0 else "bottom",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        )

    ax.set_xscale("log")
    ax.set_xlim(float(df["fold_reduction"].min()) * 0.85, float(df["fold_reduction"].max()) * 1.6)
    ax.grid(alpha=0.25)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.set_xlabel("Fold-change", fontsize=AXIS_LABEL_SIZE, fontweight="bold")
    ax.set_ylabel("SD of vdW \u0394\u0394G (kJ/mol)", fontsize=AXIS_LABEL_SIZE, fontweight="bold")
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.62, 0.98),
        fontsize=LEGEND_LABEL_SIZE,
        frameon=True,
        framealpha=0.9,
        facecolor="white",
        edgecolor="#cccccc",
    )
    ax.text(
        0.02,
        0.98,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=STATS_LABEL_SIZE,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=2.0)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

    return {
        "feature": "ddg_vdw_std",
        "label": "SD of vdW ΔΔG",
        "n_mutations": int(len(valid)),
        "pearson_r": float(r_value) if np.isfinite(r_value) else np.nan,
        "pearson_pvalue": float(p_value) if np.isfinite(p_value) else np.nan,
        "r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan,
        "slope_per_log10_fold": float(slope) if np.isfinite(slope) else np.nan,
        "intercept": float(intercept) if np.isfinite(intercept) else np.nan,
        "stderr": float(stderr) if np.isfinite(stderr) else np.nan,
        "output_png": str(output_png),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot binding-energy variability against DOR fold-change.")
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/last20frames/tables/mutation_ddg_summary.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/binding_energy/last20frames"),
    )
    args = parser.parse_args()

    if not args.summary_csv.exists():
        raise FileNotFoundError(args.summary_csv)
    summary = pd.read_csv(args.summary_csv)

    out_plots = args.output_dir / "plots"
    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    stats_row = _plot_ddg_vdw_std_vs_fold_change(
        summary,
        out_plots / "mmgbsa_ddg_vdw_std_vs_fold_change.png",
    )
    pd.DataFrame([stats_row]).to_csv(out_tables / "ddg_vdw_std_vs_fold_change_stats.csv", index=False)
    (out_config / "binding_energy_variability_config.json").write_text(
        json.dumps(
            {
                "summary_csv": str(args.summary_csv),
                "output_dir": str(args.output_dir),
                "plot": stats_row["output_png"],
                "feature": "ddg_vdw_std",
                "included_categories": ["Negative control", "Positive control"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
