#!/usr/bin/env python3
"""Plot a multi-panel figure of the most "interesting" DRM distance traces.

Selection:
  - Uses results/drm_sidechain_distance_interesting_traces.csv (precomputed scoring)
  - Picks top-N Mutant traces by interesting_score (default 9)

Traces:
  - Uses results/drm_sidechain_distance_timeseries_all_mutations.csv which contains
    both Mutant and WT baseline traces for each mutation-context.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def plot_interesting_traces(
    interesting_csv: Path,
    timeseries_csv: Path,
    output_png: Path,
    top_n: int,
) -> None:
    import matplotlib.pyplot as plt

    interesting = pd.read_csv(interesting_csv)
    if interesting.empty:
        raise ValueError(f"Empty interesting traces CSV: {interesting_csv}")

    interesting = interesting.copy()
    interesting["interesting_score"] = pd.to_numeric(interesting["interesting_score"], errors="coerce")
    interesting["is_interesting"] = interesting.get("is_interesting").astype(bool)

    # Focus on mutant traces that were flagged interesting.
    sel = interesting[(interesting["system"] == "Mutant") & (interesting["is_interesting"])].copy()
    if sel.empty:
        raise ValueError("No rows with system=Mutant and is_interesting=True.")
    sel = sel.sort_values("interesting_score", ascending=False).head(max(1, int(top_n))).reset_index(drop=True)

    ts = pd.read_csv(timeseries_csv)
    if ts.empty:
        raise ValueError(f"Empty timeseries CSV: {timeseries_csv}")

    ts["time_ns"] = pd.to_numeric(ts["time_ns"], errors="coerce")
    ts["distance_angstrom"] = pd.to_numeric(ts["distance_angstrom"], errors="coerce")

    n = len(sel)
    ncols = 3 if n >= 3 else n
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.2 * nrows), constrained_layout=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.ravel()

    for i, row in sel.iterrows():
        mut = str(row["mutation"])
        metric = str(row["metric"])
        rep = int(row["replicate"])

        ax = axes[i]
        sub_mut = ts[(ts["mutation"].astype(str) == mut) & (ts["metric"].astype(str) == metric) & (ts["replicate"] == rep)].copy()
        if sub_mut.empty:
            ax.set_axis_off()
            continue

        # Expect two systems: Mutant and WT (baseline) for the same mutation-context.
        for sys_label, style in [("WT", {"color": "#7f7f7f", "alpha": 0.7, "lw": 1.0}), ("Mutant", {"color": "#d62728", "alpha": 0.95, "lw": 1.6})]:
            g = sub_mut[sub_mut["system"].astype(str) == sys_label].sort_values("time_ns")
            if g.empty:
                continue
            ax.plot(g["time_ns"], g["distance_angstrom"], label=sys_label, **style)

        title = f"{mut} rep{rep} ({metric})"
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (Å)")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8, frameon=False)

        reasons = str(row.get("interesting_reasons", "")).strip()
        if reasons:
            ax.text(
                0.01,
                0.01,
                reasons.replace(";", "\n"),
                transform=ax.transAxes,
                fontsize=7,
                va="bottom",
                ha="left",
                color="#333333",
            )

    # Hide any unused axes.
    for j in range(n, len(axes)):
        axes[j].set_axis_off()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot top-N interesting DRM distance traces.")
    parser.add_argument("--interesting", type=Path, default=Path("results/drm_sidechain_distance_interesting_traces.csv"))
    parser.add_argument("--timeseries", type=Path, default=Path("results/drm_sidechain_distance_timeseries_all_mutations.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/plots/interesting_drm_distance_traces.png"))
    parser.add_argument("--top-n", type=int, default=9)
    args = parser.parse_args()

    for p in [args.interesting, args.timeseries]:
        if not p.exists():
            raise FileNotFoundError(f"Missing: {p}")

    plot_interesting_traces(
        interesting_csv=args.interesting,
        timeseries_csv=args.timeseries,
        output_png=args.output,
        top_n=max(1, int(args.top_n)),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

