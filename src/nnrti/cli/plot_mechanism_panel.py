#!/usr/bin/env python3
"""Three-panel molecular-mechanism figure.

One panel per mechanism, each showing the replicate-averaged trajectory of the
coordinate that best separates the genotype group from WT:

  A  Y188L   burial of the DOR chlorocyanophenyl ring
  B  V106A   displacement of DOR toward Ser105, with loss of overall packing
  C  G190    expansion of the NNIBP, graded by the size of the residue 190
             side chain

Traces are means over replicates on a common time grid truncated to the
shortest replicate; shading is the SEM across replicates, matching the
convention used for the convergence figures.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import textwrap

import pandas as pd

WT_COLOUR = "#B2182B"
GRID_POINTS = 300

PANELS = [
    dict(
        tag="A",
        column="chl_ring_burial",
        # A count of atom PAIRS within 4.0 A, not of atoms: an RT atom close to
        # three ring atoms contributes three.
        ylabel="DOR chlorocyanophenyl ring burial\n(RT heavy-atom contacts < 4.0 Å)",
        title="Y188L causes loss of aromatic packing between the DOR "
              "chlorocyanophenyl ring and Tyr188 sidechain",
        legend_loc="upper center",
        series=[("Y188L", "#08519C", 2.4, "-")],
    ),
    dict(
        tag="B",
        column="dor_to_S105_mindist",
        ylabel="DOR displacement from the crystallographic pose\n"
               "(minimum distance to Ser105, Å)",
        title="V106A causes DOR to slip out of its crystal structure pose",
        legend_loc="upper center",
        series=[
            ("V106A", "#08519C", 2.2, "-"),
            ("V106A+F227L", "#2171B5", 1.6, "-"),
            ("V106A+L234I", "#4292C6", 1.6, "-"),
            ("V106A+P225H", "#6BAED6", 1.6, "-"),
        ],
    ),
    dict(
        tag="C",
        # Not in mechanism_coordinates.csv: the pocket volume is computed by a
        # separate pass over the PBC-corrected trajectories.
        csv=Path("results/analysis/pocket_volume/tables/"
                 "holo_volume_per_frame.csv"),
        column="pocket_volume",
        ylabel="NNIBP volume (\u00c5\u00b3)",
        # Only G190E is resolved against WT (+57 +/- 10 A^3); G190A is +5 A^3
        # and is plotted as the contrast, so the title names only G190E.
        title="G190E expands the DOR binding pocket",
        legend_loc="upper center",
        # G190E rep_03 stops at 77.7 ns; carry the trace to the end of the
        # longest replicate rather than truncating the genotype there.
        extend=True,
        # Ordered by side-chain size, and shaded to match.
        series=[
            ("G190A", "#6BAED6", 1.8, "-"),
            ("G190E", "#08519C", 2.2, "-"),
        ],
    ),
]


def traces(df: pd.DataFrame, column: str, mutations: list[str],
           extend: bool = False) -> dict:
    """Replicate mean and SEM of `column` on a common time grid.

    By default the grid stops at the shortest replicate. With `extend`, it
    runs to the longest, and past a replicate's own end that replicate drops
    out of the average: the trace continues on the remaining replicates
    rather than ending, and its SEM band changes width where n changes.
    """
    out = {}
    for mutation in mutations:
        block = df[df["mutation"] == mutation]
        reps = [r for _k, r in block.groupby("replicate") if r[column].notna().any()]
        if len(reps) < 2:
            continue
        ends = [float(r["time_ns"].max()) for r in reps]
        grid = np.linspace(0.0, max(ends) if extend else min(ends), GRID_POINTS)
        rows = []
        for r in reps:
            t = r["time_ns"].to_numpy()
            y = (
                pd.to_numeric(r[column], errors="coerce")
                .interpolate(limit_direction="both")
                .to_numpy()
            )
            curve = np.interp(grid, t, y)
            # np.interp holds the last value flat beyond the data; blank it
            # instead so a finished replicate does not fake a plateau.
            curve[grid > t.max()] = np.nan
            rows.append(curve)
        curves = np.vstack(rows)
        n = np.isfinite(curves).sum(axis=0)
        with np.errstate(invalid="ignore"):
            mean = np.nanmean(curves, axis=0)
            sd = np.nanstd(curves, axis=0, ddof=1)
        sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), 0.0)
        out[mutation] = (grid, mean, sem)
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
    # A panel may read its coordinate from its own file; cache so a second
    # panel on the same file does not re-read it.
    extra: dict[Path, pd.DataFrame] = {}

    def frame_for(spec) -> pd.DataFrame:
        path = spec.get("csv")
        if path is None:
            return df
        if path not in extra:
            extra[path] = pd.read_csv(path)
        return extra[path]

    fig, axes = plt.subplots(1, len(PANELS),
                             figsize=(6.5 * len(PANELS), 5.6))

    for spec, ax in zip(PANELS, np.atleast_1d(axes).ravel()):
        wanted = ["WT"] + [s[0] for s in spec["series"]]
        tr = traces(frame_for(spec), spec["column"], wanted,
                    extend=bool(spec.get("extend")))
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
        # The titles are full sentences now, so they wrap rather than run off
        # the axes; A is two lines and B one, and both panels reserve two so
        # the two plots keep a common top edge.
        ax.set_title(textwrap.fill(spec["title"], 58), fontsize=12.5, pad=10,
                     linespacing=1.25)
        ax.tick_params(labelsize=11)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlim(0, max(g.max() for g, _m, _e in tr.values()))
        # headroom so the legend sits in blank space rather than over the traces
        lo = min(float(np.nanmin(m - e)) for _g, m, e in tr.values())
        hi = max(float(np.nanmax(m + e)) for _g, m, e in tr.values())
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
        rep = frame_for(spec).groupby(["mutation", "replicate"])[col].mean()
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
    tpath = args.output.parent / "mechanism_summary.csv"
    table.to_csv(tpath, index=False)
    print(f"wrote {tpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
