#!/usr/bin/env python3
"""WT-referenced heatmap of the Table 3 structural metrics.

Reviewer 1 (comment 4) noted that the structural results are spread across
dense tables and asked for a WT-referenced heatmap to make the patterns legible.

Each cell is the genotype's deviation from wild type expressed in standard
errors of that metric, using the pooled replicate SEM across the panel as the
scale. Standard errors rather than raw units because the metrics have
incompatible units -- angstroms, contact counts, degrees, cubic angstroms --
and a single colour scale over raw values would be dominated by the pocket
volume. Expressing everything on a common uncertainty scale also makes the
colour a reader can trust: a cell beyond about two units is a shift the
replicates support, and a pale cell is a difference within noise.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_structural_heatmap
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


def parse_cell(s: str) -> tuple[float, float]:
    """'2.97 ± 0.01' -> (2.97, 0.01); '*' and blanks -> (nan, nan)."""
    if not isinstance(s, str):
        return (float("nan"), float("nan"))
    m = re.match(r"\s*(-?[\d.]+)\s*±\s*([\d.]+)", s)
    if not m:
        return (float("nan"), float("nan"))
    return (float(m.group(1)), float(m.group(2)))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", type=Path, default=Path("paper/tables/Table-3-structural.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/mechanisms/plots"))
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)

    d = pd.read_csv(args.table)
    metric_cols = [c for c in d.columns if c not in ("Category", "Genotype")]

    means = pd.DataFrame({c: d[c].map(lambda s: parse_cell(s)[0]) for c in metric_cols})
    sems = pd.DataFrame({c: d[c].map(lambda s: parse_cell(s)[1]) for c in metric_cols})
    means.index = d.Genotype
    sems.index = d.Genotype

    wt = means.loc["WT"]
    # Pooled SEM per metric: typical replicate uncertainty across the panel.
    scale = sems.replace(0, np.nan).median()
    z = (means - wt) / scale
    z = z.drop(index="WT")
    cats = d.set_index("Genotype").Category.drop(index="WT")

    order = [g for c in ("Susceptible", "Resistant", "Uncertain") for g in cats[cats == c].index]
    z = z.loc[order]

    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    lim = float(np.nanpercentile(np.abs(z.to_numpy()), 97))
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    im = ax.imshow(z.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(len(metric_cols)))
    ax.set_xticklabels(metric_cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(z.index, fontsize=8)

    # annotate each cell, and mark those beyond 2 SEM
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            v = z.iat[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=6, color="0.4")
                continue
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(v) > lim * 0.6 else "black",
                    fontweight="bold" if abs(v) >= 2 else "normal")

    # category separators
    bounds, seen = [], None
    for i, g in enumerate(z.index):
        c = cats[g]
        if seen is not None and c != seen:
            bounds.append(i - 0.5)
        seen = c
    for b in bounds:
        ax.axhline(b, color="black", lw=1.4)

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("deviation from WT (pooled replicate SEM)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title("WT-referenced structural metrics\nbold cells differ from WT by at least 2 SEM",
                 fontsize=10)
    fig.tight_layout()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(args.output_dir / f"structural_heatmap_wt_referenced.{ext}",
                    dpi=args.dpi, bbox_inches="tight")
    tab = args.output_dir.parent / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    z.round(2).to_csv(tab / "structural_heatmap_wt_referenced.csv")

    print(f"wrote {args.output_dir}/structural_heatmap_wt_referenced.png|pdf")
    big = (z.abs() >= 2).sum().sort_values(ascending=False)
    print("\nmetrics with the most genotypes beyond 2 SEM from WT:")
    print(big.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
