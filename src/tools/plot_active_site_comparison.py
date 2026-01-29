#!/usr/bin/env python3
"""
Generate a static bar chart comparison of WT vs mutant active site metrics.

Usage:
    arch -arm64 uv run python -m src.tools.plot_active_site_comparison [--mutation Y181C_K101E]
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path


def main():
    """Main entry point."""
    ROOT = Path(__file__).resolve().parents[2]
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

    import matplotlib.pyplot as plt
    import pandas as pd

    parser = argparse.ArgumentParser(description="Plot active site comparison")
    parser.add_argument("--mutation", default="Y181C_K101E")
    args = parser.parse_args()

    csv_path = ROOT / "results" / "metrics_summary.csv"
    output_dir = ROOT / "results" / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    mutation = args.mutation
    mutation_label = mutation.replace("_", "+")

    # Load metrics
    df = pd.read_csv(csv_path)
    mut_df = df[(df["structure"] == "RPV") & (df["mutation"] == mutation_label) & (df["replicate"] == 1)]

    metrics = {"WT": {}, "MUT": {}}
    for state in ["WT", "MUT"]:
        for _, row in mut_df[mut_df["state"] == state].iterrows():
            metrics[state][row["metric"]] = row["value"]

    # Plot only pocket volume and contacts
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    wt_color, mut_color = "#4CAF50", "#FF5722"

    for ax, (metric, title) in zip(axes, [
        ("pocket_volume_proxy", "Pocket Volume (Å³)"),
        ("contact_count", "Ligand-Protein Contacts"),
    ]):
        wt_val = metrics["WT"].get(metric, 0)
        mut_val = metrics["MUT"].get(metric, 0)
        delta = mut_val - wt_val

        bars = ax.bar(["WT", mutation_label], [wt_val, mut_val], color=[wt_color, mut_color], edgecolor="black")
        for bar, val in zip(bars, [wt_val, mut_val]):
            ax.annotate(f"{val:.1f}", xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                       xytext=(0, 3), textcoords="offset points", ha="center", fontweight="bold")

        delta_color = "red" if (metric == "pocket_volume_proxy" and delta > 0) or (metric == "contact_count" and delta < 0) else "green"
        ax.text(0.5, 0.02, f"Δ = {delta:+.1f}", transform=ax.transAxes, ha="center",
               fontweight="bold", color=delta_color, bbox=dict(boxstyle="round", facecolor="white", edgecolor=delta_color))
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel("Value")

    fig.suptitle(f"RT/DNA/RPV: WT vs {mutation_label}", fontweight="bold")
    plt.tight_layout()

    output_path = output_dir / f"active_site_comparison_{mutation.lower()}.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
