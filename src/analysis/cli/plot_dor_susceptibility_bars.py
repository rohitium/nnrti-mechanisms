#!/usr/bin/env python3
"""Plot a publication-style DOR susceptibility bar chart."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from ..susceptibility import load_dor_susceptibilities

NEGATIVE_CONTROLS = {
    "K103N",
    "Y181C",
    "G190A",
    "V106I",
    "F227C",
}

POSITIVE_CONTROLS = {
    "V106A",
    "Y188L",
    "Y318F",
    "A98G+F227C",
    "V106A+F227L",
    "V106A+L234I",
    "V106A+P225H",
    "V106I+F227C",
}

UNCERTAIN_LIMITED = {
    "L100I+K103N",
    "K103N+P225H",
    "K103N+M230L",
    "V106M",
    "G190E",
    "G190S",
}

CATEGORY_COLORS = {
    "Negative control": "#4c78a8",
    "Positive control": "#e45756",
    "Uncertain/limited data": "#9aa0a6",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot DOR susceptibility bar chart from the curated spreadsheet.")
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/dor_susceptibility_bar_chart"),
    )
    parser.add_argument(
        "--star-mutation",
        type=str,
        default="V106A+L234I",
        help="Mutation to annotate with a star.",
    )
    parser.add_argument(
        "--dagger-mutation",
        type=str,
        default="F227C",
        help="Mutation to annotate with a dagger symbol.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="",
        help="Optional figure title. Leave blank for a clean publication panel.",
    )
    return parser.parse_args()


def _sort_key(label: str) -> tuple[int, str]:
    positions = [int(match) for match in re.findall(r"(\d+)", str(label))]
    if not positions:
        return (10**9, str(label).upper())
    return (min(positions), str(label).upper())


def _category_for_mutation(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation in NEGATIVE_CONTROLS:
        return "Negative control"
    if mutation in POSITIVE_CONTROLS:
        return "Positive control"
    if mutation in UNCERTAIN_LIMITED:
        return "Uncertain/limited data"
    raise ValueError(f"Missing category for mutation: {label}")


def _plot(df, output_png: Path, *, star_mutation: str, dagger_mutation: str, title: str) -> None:
    labels = df["mutation"].astype(str).tolist()
    values = df["dor_fold_reduction"].astype(float).to_numpy(dtype=float)
    colors = df["category"].map(CATEGORY_COLORS).tolist()
    x = np.arange(len(labels), dtype=float)

    ymax = float(np.nanmax(values))
    ypad = max(8.0, 0.12 * ymax)

    fig, ax = plt.subplots(figsize=(10.8, 4.6))
    bars = ax.bar(
        x,
        values,
        width=0.68,
        color=colors,
        edgecolor="#222222",
        linewidth=0.9,
        zorder=3,
    )

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#d9d9d9", linewidth=0.9)
    ax.xaxis.grid(False)

    ax.set_ylabel("Fold-change", fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=50, ha="right", fontsize=9)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(0.0, ymax + ypad)

    tick_step = 10.0 if ymax <= 100 else 20.0
    ax.set_yticks(np.arange(0.0, ymax + ypad + tick_step, tick_step))

    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

    legend_handles = [
        Patch(facecolor=CATEGORY_COLORS["Negative control"], edgecolor="#222222", label="Negative control"),
        Patch(facecolor=CATEGORY_COLORS["Positive control"], edgecolor="#222222", label="Positive control"),
        Patch(facecolor=CATEGORY_COLORS["Uncertain/limited data"], edgecolor="#222222", label="Uncertain/limited data"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=9)

    annotations = {
        str(star_mutation).strip().upper(): "*",
        str(dagger_mutation).strip().upper(): "\u2020",
    }
    for bar, label, value in zip(bars, labels, values):
        symbol = annotations.get(str(label).upper())
        if symbol is None:
            continue
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            float(value) + 0.03 * ymax,
            symbol,
            ha="center",
            va="bottom",
            fontsize=16,
            fontweight="bold",
            color="#111111",
        )

    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    df = (
        load_dor_susceptibilities(args.susceptibility_xlsx)
        .sort_values("mutation", key=lambda s: s.map(_sort_key), kind="stable")
        .reset_index(drop=True)
    )
    df["category"] = df["mutation"].map(_category_for_mutation)
    df.to_csv(out_tables / "dor_susceptibility_values.csv", index=False)

    _plot(
        df,
        out_plots / "dor_susceptibility_bar_chart.png",
        star_mutation=args.star_mutation,
        title=args.title,
        dagger_mutation=args.dagger_mutation,
    )

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "star_mutation": str(args.star_mutation),
                "dagger_mutation": str(args.dagger_mutation),
                "title": str(args.title),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
