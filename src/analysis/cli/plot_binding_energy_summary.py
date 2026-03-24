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


def _place_greedy_annotations(
    ax,
    x_values: np.ndarray,
    y_values: np.ndarray,
    labels: list[str],
    *,
    text_color: str = "#333333",
) -> None:
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    candidate_offsets = [
        (8, 8),
        (8, 18),
        (8, -18),
        (-8, 8),
        (-8, 18),
        (-8, -18),
        (14, 0),
        (-14, 0),
        (14, 22),
        (-14, 22),
        (14, -22),
        (-14, -22),
        (24, 8),
        (-24, 8),
        (24, -8),
        (-24, -8),
        (30, 18),
        (-30, 18),
        (30, -18),
        (-30, -18),
    ]
    placed_boxes: list[tuple[float, float, float, float]] = []

    order = np.argsort(x_values)
    for idx in order:
        x = float(x_values[idx])
        y = float(y_values[idx])
        label = str(labels[idx])
        anchor_px = ax.transData.transform((x, y))
        anchor_disp_x = float(anchor_px[0])
        anchor_disp_y = float(anchor_px[1])

        best = None
        best_score = None
        for dx_pt, dy_pt in candidate_offsets:
            ann = ax.annotate(
                label,
                xy=(x, y),
                xytext=(dx_pt, dy_pt),
                textcoords="offset points",
                fontsize=9,
                color=text_color,
                alpha=0.92,
                ha="left" if dx_pt >= 0 else "right",
                va="bottom" if dy_pt >= 0 else "top",
                arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
                bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
            )
            bbox = ann.get_window_extent(renderer=renderer)
            ann.remove()

            overlaps = 0
            overlap_area = 0.0
            for other in placed_boxes:
                x0 = max(float(bbox.x0), other[0])
                y0 = max(float(bbox.y0), other[1])
                x1 = min(float(bbox.x1), other[2])
                y1 = min(float(bbox.y1), other[3])
                if x1 > x0 and y1 > y0:
                    overlaps += 1
                    overlap_area += (x1 - x0) * (y1 - y0)
            center_x = 0.5 * (float(bbox.x0) + float(bbox.x1))
            center_y = 0.5 * (float(bbox.y0) + float(bbox.y1))
            anchor_penalty = abs(center_x - anchor_disp_x) + abs(center_y - anchor_disp_y)
            score = (overlaps * 1_000_000.0) + (overlap_area * 100.0) + anchor_penalty
            if best_score is None or score < best_score:
                best_score = score
                best = (dx_pt, dy_pt, bbox)

        assert best is not None
        dx_pt, dy_pt, bbox = best
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(dx_pt, dy_pt),
            textcoords="offset points",
            fontsize=9,
            color=text_color,
            alpha=0.92,
            ha="left" if dx_pt >= 0 else "right",
            va="bottom" if dy_pt >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        placed_boxes.append((float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)))


def _plot_ddg_bar_chart(summary: pd.DataFrame, specs: list[tuple[str, str, str]], title: str, output_png: Path) -> None:
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
        facecolors = ["#d62828" if bool(flag) else "#1d3557" for flag in ordered["is_combo"].tolist()]
        ax.bar(x, y, color=facecolors, edgecolor=color, linewidth=0.8, alpha=0.82, zorder=2)
        ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#333333", elinewidth=1.1, capsize=2.5, alpha=0.95, zorder=3)
        ax.grid(axis="y", linestyle=":", alpha=0.32)
        ax.set_ylabel(f"{label}\n(kJ/mol)", fontsize=9)
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


def _plot_component_vs_fold_change(
    summary: pd.DataFrame,
    column: str,
    label: str,
    color: str,
    output_png: Path,
) -> None:
    df = summary[summary["mutation"].astype(str) != "WT"].copy()
    df["log10_fold_reduction"] = np.log10(pd.to_numeric(df["fold_reduction"], errors="coerce"))
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_fold_reduction"]).reset_index(drop=True)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14.2, 9.2))
    singles = df[~df["is_combo"]].copy()
    combos = df[df["is_combo"]].copy()
    for subset, point_color, legend_label, marker in [
        (singles, "#1d3557", "Single DRM", "o"),
        (combos, "#d62828", "Combination DRM", "s"),
    ]:
        if subset.empty:
            continue
        ax.errorbar(
            subset["fold_reduction"],
            subset[f"{column}_mean"],
            yerr=subset[f"{column}_sem"],
            fmt=marker,
            ms=7.5,
            lw=0,
            elinewidth=1.15,
            capsize=2.5,
            color=point_color,
            markerfacecolor=point_color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            alpha=0.95,
            label=legend_label,
            zorder=3,
        )

    valid = df[["log10_fold_reduction", f"{column}_mean"]].dropna()
    if len(valid) >= 3:
        slope, intercept, r_value, p_value, _stderr = stats.linregress(
            valid["log10_fold_reduction"],
            valid[f"{column}_mean"],
        )
        x_grid = np.geomspace(float(df["fold_reduction"].min()) * 0.9, float(df["fold_reduction"].max()) * 1.1, 300)
        y_grid = slope * np.log10(x_grid) + intercept
        ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.5, zorder=2)
        annotation = f"R^2 = {r_value**2:.3f}\np = {p_value:.3f}"
    else:
        annotation = "R^2 = NA\np = NA"

    ax.set_xscale("log")
    ax.grid(alpha=0.25)
    ax.set_title(f"{label} Vs Fold Reduction", fontsize=14, fontweight="bold")
    ax.set_xlabel("Fold Reduction")
    ax.set_ylabel("kJ/mol")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, fontsize=10, frameon=False)
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

    _plot_ddg_bar_chart(ddg_summary, DDG_SPECS, "MM/GBSA Components Relative To WT", out_plots / "mmgbsa_ddg_components_vs_wt.png")
    for column, label, color in COMPONENT_SPECS:
        slug = column.replace("binding_dg_", "").replace("binding_dg", "total")
        _plot_component_vs_fold_change(
            rep_summary,
            column,
            label,
            color,
            out_plots / f"mmgbsa_{slug}_vs_fold_change.png",
        )

    for stale in [
        out_plots / "mmgbsa_components_by_mutation.png",
        out_plots / "mmgbsa_components_vs_fold_change.png",
        out_plots / "mmgbsa_ddg_vs_fold_change.png",
    ]:
        if stale.exists():
            stale.unlink()

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
