#!/usr/bin/env python3
"""Plot WT-referenced contact-occupancy shifts as stacked tick-line tracks."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..susceptibility import load_dor_susceptibilities
from .plot_triplet_contact_story import _normalize_mutation_token


SUSCEPTIBLE = {"V106I", "K103N", "Y181C", "G190A"}
RESISTANT = {
    "V106A",
    "Y188L",
    "Y318F",
    "A98G+F227C",
    "V106A+F227L",
    "V106A+L234I",
    "V106A+P225H",
    "V106I+F227C",
    "K103N+M230L",
}
UNCERTAIN = {"L100I+K103N", "K103N+P225H", "V106M", "G190E", "G190S"}

CATEGORY_ORDER = {"Susceptible": 0, "Uncertain": 1, "Resistant": 2}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot all-mutation WT-referenced occupancy shifts as stacked tick-line tracks."
    )
    parser.add_argument(
        "--matrix-csv",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/"
            "all_mutation_wt_referenced_occupancy_heatmap_wt_contacted_residues_by_region_excluding_f227c_matrix.csv"
        ),
    )
    parser.add_argument(
        "--residue-csv",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/"
            "all_mutation_wt_referenced_occupancy_heatmap_wt_contacted_residues_by_region_excluding_f227c_display_residues.csv"
        ),
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/plots/"
            "all_mutation_wt_referenced_occupancy_tick_lines_wt_contacted_residues_by_region_excluding_f227c.png"
        ),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/"
            "all_mutation_wt_referenced_occupancy_tick_lines_thresholded_matrix.csv"
        ),
    )
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument("--zero-threshold", type=float, default=0.1)
    parser.add_argument("--min-tick-height", type=float, default=0.13)
    parser.add_argument("--max-tick-height", type=float, default=0.42)
    parser.add_argument("--tick-linewidth", type=float, default=4.0)
    parser.add_argument("--dpi", type=int, default=320)
    return parser.parse_args()


def _category_for_mutation(mutation: str) -> str:
    label = str(mutation).strip()
    if label in SUSCEPTIBLE:
        return "Susceptible"
    if label in RESISTANT:
        return "Resistant"
    if label in UNCERTAIN:
        return "Uncertain"
    raise ValueError(f"Mutation is not assigned to susceptible/resistant/uncertain: {mutation}")


def _fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    return {
        _normalize_mutation_token(str(row["mutation"])): float(row["dor_fold_reduction"])
        for _, row in df.iterrows()
    }


def _mutation_sort_key(mutation: str, fold_map: dict[str, float]) -> tuple[int, float, str]:
    return (
        CATEGORY_ORDER[_category_for_mutation(mutation)],
        float(fold_map.get(_normalize_mutation_token(mutation), np.inf)),
        str(mutation),
    )


def _format_mutation_label(mutation: str, fold_map: dict[str, float]) -> str:
    fold = fold_map.get(_normalize_mutation_token(mutation))
    if fold is None or not np.isfinite(fold):
        return str(mutation)
    return f"{mutation} ({fold:g}x)"


def _plot_tick_lines(
    matrix: pd.DataFrame,
    residue_df: pd.DataFrame,
    fold_map: dict[str, float],
    output_png: Path,
    zero_threshold: float,
    min_tick_height: float,
    max_tick_height: float,
    tick_linewidth: float,
    dpi: int,
) -> None:
    import matplotlib.pyplot as plt
    import textwrap
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.transforms import blended_transform_factory

    matrix = matrix.copy()
    ordered_mutations = sorted(matrix.index.tolist(), key=lambda m: _mutation_sort_key(m, fold_map=fold_map))
    matrix = matrix.loc[ordered_mutations, residue_df["label"].tolist()]

    arr = matrix.to_numpy(dtype=float)
    thresholded = arr.copy()
    thresholded[np.abs(thresholded) < float(zero_threshold)] = 0.0

    n_rows, n_cols = thresholded.shape
    x = np.arange(n_cols, dtype=float)
    y_base = np.arange(n_rows, dtype=float)

    fig_w = max(15.0, 4.6 + 0.55 * n_cols)
    fig_h = max(8.2, 1.8 + 0.42 * n_rows)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("coolwarm")
    norm = Normalize(vmin=-1.0, vmax=1.0)

    for row_i, mutation in enumerate(matrix.index.tolist()):
        ax.hlines(y_base[row_i], -0.45, n_cols - 0.55, color="#d8d8d8", linewidth=1.25, zorder=1)
        for col_i, v in enumerate(thresholded[row_i, :]):
            if not np.isfinite(v) or float(v) == 0.0:
                continue
            tick_abs = float(min_tick_height) + (
                (abs(float(v)) - float(zero_threshold))
                / max(1.0 - float(zero_threshold), 1e-9)
                * (float(max_tick_height) - float(min_tick_height))
            )
            tick = np.sign(float(v)) * tick_abs
            ax.vlines(
                x[col_i],
                y_base[row_i],
                y_base[row_i] - tick,
                color=cmap(norm(float(v))),
                linewidth=float(tick_linewidth),
                zorder=4,
            )

    ax.axhline(-0.5, color="#eeeeee", linewidth=1.0)
    for i in range(n_rows - 1):
        cat_a = _category_for_mutation(matrix.index[i])
        cat_b = _category_for_mutation(matrix.index[i + 1])
        if cat_a != cat_b:
            ax.axhline(i + 0.5, color="#555555", linewidth=2.0, zorder=2)

    ax.set_xlim(-0.45, n_cols - 0.45)
    ax.set_ylim(n_rows - 0.45, -1.1)
    ax.set_yticks(y_base)
    ax.set_yticklabels([_format_mutation_label(m, fold_map=fold_map) for m in matrix.index.tolist()], fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(residue_df["label"].tolist(), rotation=45, ha="right", fontsize=11)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)

    if {"pocket_region", "region_color"}.issubset(residue_df.columns):
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        for region, g in residue_df.groupby("pocket_region", sort=False):
            left = int(g.index.min())
            right = int(g.index.max())
            center = (left + right) / 2.0
            color = str(g["region_color"].iloc[0])
            label = "\n".join(textwrap.wrap(str(region), width=17))
            ax.plot(
                [left - 0.45, right + 0.45],
                [1.02, 1.02],
                color=color,
                linewidth=4.0,
                transform=trans,
                clip_on=False,
            )
            ax.text(
                center,
                1.055,
                label,
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
                color=color,
                transform=trans,
                clip_on=False,
            )
            if right < len(residue_df) - 1:
                ax.axvline(right + 0.5, color="#efefef", linewidth=1.0, zorder=0)

    ytrans = blended_transform_factory(ax.transAxes, ax.transData)
    row_categories = pd.Series([_category_for_mutation(m) for m in matrix.index], index=matrix.index)
    for category, rows in row_categories.groupby(row_categories, sort=False):
        positions = [matrix.index.get_loc(idx) for idx in rows.index]
        center = (min(positions) + max(positions)) / 2.0
        ax.text(
            1.015,
            center,
            str(category),
            ha="center",
            va="center",
            rotation=90,
            fontsize=12,
            fontweight="bold",
            color="#555555",
            transform=ytrans,
            clip_on=False,
        )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.028, pad=0.09)
    cbar.set_label("Δ occupancy vs WT", fontsize=12)
    cbar.set_ticks(np.arange(-1.0, 1.01, 0.2))
    cbar.ax.tick_params(labelsize=11)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    matrix = pd.read_csv(args.matrix_csv)
    if "mutation" not in matrix.columns:
        raise ValueError(f"{args.matrix_csv} must contain a mutation column.")
    matrix = matrix.set_index("mutation")
    residue_df = pd.read_csv(args.residue_csv)
    fold_map = _fold_map(args.susceptibility_xlsx)
    required = {"label", "pocket_region", "region_color"}
    missing = required.difference(residue_df.columns)
    if missing:
        raise ValueError(f"{args.residue_csv} is missing required columns: {sorted(missing)}")

    thresholded = matrix.copy()
    thresholded = thresholded.apply(pd.to_numeric, errors="coerce")
    thresholded[thresholded.abs() < float(args.zero_threshold)] = 0.0
    kept_residues = thresholded.columns[(thresholded != 0.0).any(axis=0)].tolist()
    if not kept_residues:
        raise ValueError(
            f"No residue columns have |Δ occupancy| >= {float(args.zero_threshold):.2f}; "
            "lower --zero-threshold or inspect the input matrix."
        )
    matrix = matrix[kept_residues].copy()
    thresholded = thresholded[kept_residues].copy()
    residue_df = residue_df[residue_df["label"].isin(kept_residues)].copy()
    residue_df["_display_order"] = residue_df["label"].map({label: i for i, label in enumerate(kept_residues)})
    residue_df = residue_df.sort_values("_display_order", kind="stable").drop(columns="_display_order").reset_index(drop=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    thresholded.reset_index().rename(columns={"index": "mutation"}).to_csv(args.output_csv, index=False)

    _plot_tick_lines(
        matrix=matrix.apply(pd.to_numeric, errors="coerce"),
        residue_df=residue_df,
        fold_map=fold_map,
        output_png=args.output_png,
        zero_threshold=float(args.zero_threshold),
        min_tick_height=float(args.min_tick_height),
        max_tick_height=float(args.max_tick_height),
        tick_linewidth=float(args.tick_linewidth),
        dpi=int(args.dpi),
    )

    config_path = args.output_png.parent.parent / "config" / "wt_referenced_occupancy_tick_lines_config.json"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(
            {
                "matrix_csv": str(args.matrix_csv),
                "residue_csv": str(args.residue_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_png": str(args.output_png),
                "output_csv": str(args.output_csv),
                "zero_threshold": float(args.zero_threshold),
                "min_tick_height": float(args.min_tick_height),
                "max_tick_height": float(args.max_tick_height),
                "tick_linewidth": float(args.tick_linewidth),
                "n_display_residues": int(len(kept_residues)),
                "display_residues": kept_residues,
                "category_labels": ["Susceptible", "Uncertain", "Resistant"],
            },
            indent=2,
        )
    )
    print(f"Saved {args.output_png}")
    print(f"Saved {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
