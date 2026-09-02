#!/usr/bin/env python3
"""Two-panel molecular-mechanism figure.

One panel per mechanism, each showing the replicate-averaged trajectory of the
coordinate that best separates the genotype group from WT:

  A  Y188L   burial of the DOR chlorocyanophenyl ring
  B  V106A   displacement of DOR toward Ser105, with loss of overall packing

Traces are means over replicates on a common time grid truncated to the
shortest replicate; shading is the SEM across replicates, matching the
convention used for the convergence figures.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

WT_COLOUR = "#B2182B"
GRID_POINTS = 300

PANELS = [
    dict(
        tag="A",
        column="chl_ring_burial",
        # Two corrections: the cutoff is 4.0 A (changed 2026-08-31), and the
        # quantity is a count of atom PAIRS within it, not of atoms -- an RT atom
        # close to three ring atoms contributes three. See
        # paper/contact-cutoff-sensitivity.md.
        ylabel="DOR ring burial\n(RT heavy-atom contacts < 4.0 Å)",
        title="Y188L: loss of aromatic packing",
        legend_loc="upper center",
        series=[("Y188L", "#08519C", 2.4, "-")],
    ),
    dict(
        tag="B",
        column="dor_to_S105_mindist",
        ylabel="DOR displacement\n(distance to Ser105, Å)",
        title="V106A: DOR slips out of position",
        legend_loc="upper center",
        series=[
            ("V106A", "#08519C", 2.2, "-"),
            ("V106A+F227L", "#2171B5", 1.6, "-"),
            ("V106A+L234I", "#4292C6", 1.6, "-"),
            ("V106A+P225H", "#6BAED6", 1.6, "-"),
        ],
    ),
]


def traces(df: pd.DataFrame, column: str, mutations: list[str]) -> dict:
    out = {}
    for mutation in mutations:
        block = df[df["mutation"] == mutation]
        reps = [r for _k, r in block.groupby("replicate") if r[column].notna().any()]
        if len(reps) < 2:
            continue
        t_end = min(float(r["time_ns"].max()) for r in reps)
        grid = np.linspace(0.0, t_end, GRID_POINTS)
        curves = np.vstack(
            [
                np.interp(
                    grid,
                    r["time_ns"].to_numpy(),
                    pd.to_numeric(r[column], errors="coerce")
                    .interpolate(limit_direction="both")
                    .to_numpy(),
                )
                for r in reps
            ]
        )
        n = curves.shape[0]
        out[mutation] = (grid, curves.mean(axis=0), curves.std(axis=0, ddof=1) / np.sqrt(n))
    return out


def main() -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = argparse.ArgumentParser()
    p.add_argument(
        "--csv",
        type=Path,
        default=Path("results/analysis/mechanisms/mechanism_coordinates.csv"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/analysis/mechanisms/plots/mechanism_panel.png"),
    )
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    for spec, ax in zip(PANELS, axes.ravel()):
        wanted = ["WT"] + [s[0] for s in spec["series"]]
        tr = traces(df, spec["column"], wanted)
        for name, colour, lw, ls in spec["series"]:
            if name not in tr:
                continue
            g, m, e = tr[name]
            ax.fill_between(g, m - e, m + e, color=colour, alpha=0.13, lw=0, zorder=2)
            ax.plot(g, m, color=colour, lw=lw, ls=ls, zorder=3, label=name)
        if "WT" in tr:
            g, m, e = tr["WT"]
            ax.fill_between(g, m - e, m + e, color=WT_COLOUR, alpha=0.18, lw=0, zorder=4)
            ax.plot(g, m, color=WT_COLOUR, lw=2.6, zorder=5, label="WT")

        ax.set_xlabel("Time (ns)", fontsize=12)
        ax.set_ylabel(spec["ylabel"], fontsize=12)
        ax.set_title(spec["title"], fontsize=13, pad=10)
        ax.tick_params(labelsize=11)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlim(0, max(g.max() for g, _m, _e in tr.values()))
        # headroom so the legend sits in blank space rather than over the traces
        lo = min(float((m - e).min()) for _g, m, e in tr.values())
        hi = max(float((m + e).max()) for _g, m, e in tr.values())
        span = hi - lo
        ax.set_ylim(lo - 0.06 * span, hi + 0.30 * span)
        ax.legend(fontsize=8.5, frameon=False, loc="upper center", ncol=3,
                  columnspacing=1.2, handlelength=1.8)
        ax.text(
            -0.13,
            1.06,
            spec["tag"],
            transform=ax.transAxes,
            fontsize=17,
            fontweight="bold",
            va="top",
        )

    fig.tight_layout(w_pad=3.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {args.output}")

    # companion table: per-genotype mean +/- SEM over replicates
    rows = []
    for spec in PANELS:
        col = spec["column"]
        rep = df.groupby(["mutation", "replicate"])[col].mean()
        m = rep.groupby("mutation").mean()
        s = rep.groupby("mutation").std() / np.sqrt(3)
        for mutation in ["WT"] + [x[0] for x in spec["series"]]:
            if mutation in m.index and np.isfinite(m[mutation]):
                rows.append(
                    dict(
                        panel=spec["tag"],
                        coordinate=col,
                        mutation=mutation,
                        mean=round(float(m[mutation]), 2),
                        sem=round(float(s[mutation]), 2),
                    )
                )
    table = pd.DataFrame(rows)
    tpath = args.output.parent.parent / "mechanism_summary.csv"
    table.to_csv(tpath, index=False)
    print(f"wrote {tpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
