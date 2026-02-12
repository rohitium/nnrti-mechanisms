#!/usr/bin/env python3
"""Generate resistance heatmap showing single and combination DRM effects.

Diagonal: Single DRM effects
Off-diagonal: Combination DRM effects (epistasis)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.analysis.susceptibility import load_dor_susceptibilities


def parse_mutation(mut_str: str) -> tuple[str, ...]:
    """Parse mutation string into individual DRMs.

    Examples:
        Y188L -> (Y188L,)
        K103N+P225H -> (K103N, P225H)
        V106A+L234I -> (V106A, L234I)
    """
    if "+" in mut_str:
        return tuple(part.strip() for part in mut_str.split("+"))
    return (mut_str.strip(),)


def build_resistance_matrix(df: pd.DataFrame) -> tuple[np.ndarray, list[str], dict]:
    """Build resistance matrix for heatmap.

    Returns:
        matrix: NxN array of log10(fold reduction), NaN for missing combinations
        drm_order: List of ordered DRM names (for axis labels)
        metadata: Dict with additional info (fold values, mutation counts, etc.)
    """
    # Extract all unique single DRMs
    all_drms = set()
    for _, row in df.iterrows():
        mut = str(row["mutation"])
        if mut == "WT":
            continue
        drms = parse_mutation(mut)
        all_drms.update(drms)

    # Sort DRMs by position number
    import re
    def extract_position(drm: str) -> int:
        match = re.search(r"(\d+)", drm)
        return int(match.group(1)) if match else 999

    drm_order = sorted(all_drms, key=extract_position)
    n = len(drm_order)

    # Initialize matrix with NaN
    matrix = np.full((n, n), np.nan)
    fold_matrix = np.full((n, n), np.nan)

    # Build lookup from DRM tuple to fold reduction
    mut_to_fold = {}
    for _, row in df.iterrows():
        mut = str(row["mutation"])
        if mut == "WT":
            continue
        drms = parse_mutation(mut)
        fold = float(row["fold_reduction"])
        mut_to_fold[drms] = fold

    # Fill matrix
    for i, drm1 in enumerate(drm_order):
        for j, drm2 in enumerate(drm_order):
            if i == j:
                # Diagonal: single DRM
                if (drm1,) in mut_to_fold:
                    fold = mut_to_fold[(drm1,)]
                    fold_matrix[i, j] = fold
                    # Use log10(fold) for color scale
                    if fold > 0:
                        matrix[i, j] = np.log10(fold)
            else:
                # Off-diagonal: check for combination in both orders
                combo1 = tuple(sorted([drm1, drm2]))
                combo2 = (drm1, drm2)
                combo3 = (drm2, drm1)

                for combo in [combo1, combo2, combo3]:
                    if combo in mut_to_fold:
                        fold = mut_to_fold[combo]
                        fold_matrix[i, j] = fold
                        if fold > 0:
                            matrix[i, j] = np.log10(fold)
                        break

    metadata = {
        "fold_matrix": fold_matrix,
        "drm_order": drm_order,
        "n_drms": n,
    }

    return matrix, drm_order, metadata


def plot_resistance_heatmap(
    matrix: np.ndarray,
    drm_order: list[str],
    metadata: dict,
    output_path: Path,
) -> None:
    """Create resistance heatmap visualization."""
    import matplotlib.pyplot as plt

    n = len(drm_order)
    fold_matrix = metadata["fold_matrix"]

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Create custom colormap: white for NaN, blue-yellow-red for values
    from matplotlib.colors import LinearSegmentedColormap
    colors = ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"]
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("resistance", colors, N=n_bins)
    cmap.set_bad(color="#f0f0f0")  # Gray for missing data

    # Plot heatmap
    im = ax.imshow(
        matrix,
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
        vmin=np.nanmin(matrix),
        vmax=np.nanmax(matrix),
    )

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log₁₀(Fold Reduction)", fontsize=12, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    # Set ticks and labels
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(drm_order, fontsize=10, rotation=45, ha="left")
    ax.set_yticklabels(drm_order, fontsize=10)

    # Move x-axis labels to top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    # Add grid
    ax.set_xticks(np.arange(n) - 0.5, minor=True)
    ax.set_yticks(np.arange(n) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", size=0)

    # Annotate cells with fold reduction values
    for i in range(n):
        for j in range(n):
            fold = fold_matrix[i, j]
            if not np.isnan(fold):
                # Determine text color based on background
                log_val = matrix[i, j]
                text_color = "white" if log_val > np.nanmedian(matrix) else "black"

                # Format value
                if fold >= 10:
                    text = f"{fold:.0f}×"
                else:
                    text = f"{fold:.1f}×"

                ax.text(
                    j, i, text,
                    ha="center", va="center",
                    fontsize=9,
                    fontweight="bold",
                    color=text_color,
                )

    # Title and labels
    ax.set_title(
        "Doravirine Resistance Landscape",
        fontsize=14,
        fontweight="bold",
        pad=20,
        y=1.08,  # Move title up to make room for top x-axis labels
    )
    ax.set_xlabel("DRM", fontsize=12, fontweight="bold")
    ax.set_ylabel("DRM", fontsize=12, fontweight="bold")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot resistance heatmap")
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument("--output", type=Path, default=Path("results/plots/resistance_heatmap.png"))
    args = parser.parse_args()

    if not args.susceptibility_xlsx.exists():
        print(f"Error: {args.susceptibility_xlsx} not found")
        return 1

    susc_df = load_dor_susceptibilities(args.susceptibility_xlsx, default_chain="A")
    df = susc_df.rename(columns={"dor_fold_reduction": "fold_reduction"})[
        ["mutation", "fold_reduction"]
    ].copy()

    # Build resistance matrix
    matrix, drm_order, metadata = build_resistance_matrix(df)

    # Plot heatmap
    plot_resistance_heatmap(matrix, drm_order, metadata, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
