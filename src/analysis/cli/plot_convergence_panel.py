#!/usr/bin/env python3
"""Single-panel convergence figures for the whole genotype panel.

One mean trace per genotype with a shaded across-replicate SEM band, coloured by
DOR fold reduction so the reader can see that the metric carries no resistance
ordering. Fold values come from the authoritative susceptibility spreadsheet,
never from the manifest.

``--metric`` selects the quantity: backbone RMSD, ligand RMSD, or the
ligand-pocket centre-of-mass separation. The first two start at zero by
construction; the COM distance is an absolute separation and does not.

Replicates run to slightly different lengths, so each genotype's trace stops at
its shortest replicate: the mean is always over n = 3, never a changing
denominator.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


#: Per-metric presentation. The inset gets a compact label so it cannot overhang
#: into the main panel.
METRICS = {
    "ca_rmsd": {
        "column": "ca_rmsd_angstrom",
        "ylabel": "C\u03b1 RMSD (\u00c5)",
        "inset_ylabel": "C\u03b1 RMSD (\u00c5)",
        "title": "RT\u2013DOR backbone convergence across the genotype panel",
        "stem": "rmsd_convergence",
        "inset_rect": (0.575, 0.125, 0.40, 0.285),
    },
    "dor_rmsd": {
        "column": "dor_rmsd_angstrom",
        "ylabel": "DOR RMSD (\u00c5)",
        "inset_ylabel": "DOR RMSD (\u00c5)",
        "title": "Doravirine pose convergence across the genotype panel",
        "stem": "dor_rmsd_convergence",
        "inset_rect": (0.575, 0.125, 0.40, 0.285),
    },
    "com_distance": {
        "column": "com_distance_angstrom",
        "ylabel": "DOR\u2013pocket COM distance (\u00c5)",
        "inset_ylabel": "COM dist. (\u00c5)",
        "title": "Doravirine\u2013pocket separation across the genotype panel",
        "stem": "com_distance_convergence",
        "inset_rect": (0.585, 0.125, 0.39, 0.285),
    },
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--rmsd-csv", type=Path,
        default=root / "results/analysis/md_convergence/tables/frame_traces.csv",
        help=(
            "Per-frame traces. The md_convergence table references each replicate to its own "
            "first production frame, so traces start at zero. Note results/rmsd_ca_profiles.csv "
            "instead references the setup topology, which precedes heating and equilibration, so "
            "its traces start 1.5-3.8 A off zero."
        ),
    )
    p.add_argument("--susceptibility-xlsx", type=Path,
                   default=root / "data/DRM-susceptibilities.csv.xlsx")
    p.add_argument("--metric", choices=sorted(METRICS), default="ca_rmsd")
    p.add_argument("--output", type=Path, default=None,
                   help="Defaults to the metric's conventional filename under md_convergence/plots.")
    p.add_argument("--inset-headroom", type=float, default=0.62,
                   help=(
                       "Blank vertical space added below the data, as a fraction of the data "
                       "range, so the inset does not overlap the traces. 0 disables it."
                   ))
    p.add_argument("--inset-rect", type=float, nargs=4, default=None,
                   metavar=("X0", "Y0", "W", "H"),
                   help="Inset position in axes fractions; defaults per metric.")
    p.add_argument("--scoring-window-ns", type=float, default=10.0,
                   help="Width of the terminal window used for MM/GBSA scoring, shaded on the figure.")
    p.add_argument(
        "--exclude", nargs="*", default=["F227C"],
        help=(
            "Genotypes to omit. F227C by default: its replicate 2 RMSD profile is an "
            "alignment/PBC artifact (the complex translates a box length, giving a 107 A "
            "Ca RMSD while Rg and inter-chain separation are unchanged, so the structure is "
            "intact and only the profile is wrong). F227C is also outside the 18-genotype "
            "manuscript panel."
        ),
    )
    p.add_argument("--inset-ns", type=float, default=5.0,
                   help="Span of the magnified inset showing the equilibration rise.")
    p.add_argument("--grid-points", type=int, default=400)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def load_fold_map(xlsx: Path) -> dict[str, float]:
    table = pd.read_excel(xlsx, header=1)
    table.columns = ["mutation", "rpv_fold", "dor_fold"]
    table = table.dropna(subset=["mutation", "dor_fold"])
    keys = table["mutation"].astype(str).str.replace(r",\s*", "+", regex=True).str.strip()
    return dict(zip(keys, table["dor_fold"].astype(float)))


def genotype_traces(rmsd: pd.DataFrame, grid_points: int, column: str):
    """Mean and SEM on a common grid, truncated to each genotype's shortest replicate."""
    out = {}
    if "time_ns" not in rmsd.columns:
        rmsd = rmsd.assign(time_ns=rmsd["time_ps"] / 1000.0)
    for mutation, block in rmsd.groupby("mutation"):
        reps = list(block.groupby("replicate"))
        if len(reps) < 2:
            continue
        t_end = min(r["time_ns"].max() for _rep, r in reps)
        grid = np.linspace(0.0, t_end, grid_points)
        curves = np.vstack([
            np.interp(grid, r["time_ns"].to_numpy(), r[column].to_numpy())
            for _rep, r in reps
        ])
        n = curves.shape[0]
        out[mutation] = (grid, curves.mean(axis=0), curves.std(axis=0, ddof=1) / np.sqrt(n), n, t_end)
    return out


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable

    args = parse_args()
    spec = METRICS[args.metric]
    column, ylabel, title = spec["column"], spec["ylabel"], spec["title"]
    stem, default_rect = spec["stem"], spec["inset_rect"]
    output = args.output or (
        repo_root() / "results/analysis/md_convergence/plots" / f"{stem}.png"
    )
    inset_rect = tuple(args.inset_rect) if args.inset_rect else default_rect

    rmsd = pd.read_csv(args.rmsd_csv)
    if column not in rmsd.columns:
        raise ValueError(f"{args.rmsd_csv} has no column {column!r}")
    if args.exclude:
        rmsd = rmsd[~rmsd["mutation"].isin(args.exclude)]
        print(f"excluded: {', '.join(args.exclude)}")
    folds = load_fold_map(args.susceptibility_xlsx)
    traces = genotype_traces(rmsd, args.grid_points, column)

    mutants = sorted((m for m in traces if m != "WT"), key=lambda m: folds.get(m, np.nan))
    values = np.array([folds[m] for m in mutants if m in folds], dtype=float)
    norm = LogNorm(vmin=values.min(), vmax=values.max())
    # Single-hue sequential ramp: magnitude, light -> dark. Truncated so the
    # palest series is still legible on white.
    ramp = plt.get_cmap("Blues")
    colour = lambda v: ramp(0.35 + 0.65 * norm(v))

    fig, ax = plt.subplots(figsize=(10.0, 6.0))

    t_max = max(v[4] for v in traces.values())

    for mutation in mutants:
        grid, mean, sem, _n, _end = traces[mutation]
        c = colour(folds[mutation]) if mutation in folds else "0.6"
        ax.fill_between(grid, mean - sem, mean + sem, color=c, alpha=0.10, lw=0, zorder=2)
        ax.plot(grid, mean, color=c, lw=1.4, zorder=3)

    if "WT" in traces:
        grid, mean, sem, _n, _end = traces["WT"]
        ax.fill_between(grid, mean - sem, mean + sem, color="#B2182B", alpha=0.18, lw=0, zorder=4)
        ax.plot(grid, mean, color="#B2182B", lw=2.6, zorder=5, label="WT")

    ax.set_xlabel("Time (ns)", fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=15, pad=12)
    ax.tick_params(labelsize=12)
    ax.set_xlim(0, t_max)
    lo = min(float((m - e).min()) for _g, m, e, _n, _t in traces.values())
    hi = max(float((m + e).max()) for _g, m, e, _n, _t in traces.values())
    pad = 0.08 * (hi - lo)
    data_floor = lo - pad
    # Extra room below the data so the inset sits in blank space instead of over
    # the traces. Ticks are not drawn into that band -- for RMSD it lies below
    # zero, where a label would imply a negative distance.
    ax.set_ylim(data_floor - args.inset_headroom * (hi - lo), hi + pad)
    ax.set_yticks([t for t in ax.get_yticks() if t >= data_floor])
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=11, frameon=False, loc="upper left")

    # Inset: the equilibration rise, which is compressed against the y-axis in the
    # main panel. Same series, same colours, no legend -- it is a magnification.
    inset = ax.inset_axes(list(inset_rect))
    for mutation in mutants:
        grid, mean, sem, _n, _end = traces[mutation]
        keep = grid <= args.inset_ns
        c = colour(folds[mutation]) if mutation in folds else "0.6"
        inset.fill_between(grid[keep], (mean - sem)[keep], (mean + sem)[keep], color=c, alpha=0.10, lw=0)
        inset.plot(grid[keep], mean[keep], color=c, lw=1.1)
    if "WT" in traces:
        grid, mean, sem, _n, _end = traces["WT"]
        keep = grid <= args.inset_ns
        inset.fill_between(grid[keep], (mean - sem)[keep], (mean + sem)[keep],
                           color="#B2182B", alpha=0.18, lw=0)
        inset.plot(grid[keep], mean[keep], color="#B2182B", lw=1.9)
    inset.set_xlim(0, args.inset_ns)
    inset.tick_params(labelsize=8, length=2.5, pad=1.5, colors="black")
    inset.set_xlabel("Time (ns)", fontsize=8.5, labelpad=1.0)
    inset.set_ylabel(spec["inset_ylabel"], fontsize=8.5, labelpad=1.0)
    inset.set_xticks([0, 1, 2, 3, 4, 5])
    inset.locator_params(axis="y", nbins=4)
    # Label inside the frame so it cannot collide with the main panel.
    inset.text(0.035, 0.93, f"first {args.inset_ns:g} ns", transform=inset.transAxes,
               fontsize=8.5, color="black", ha="left", va="top")
    inset.grid(alpha=0.2, lw=0.5)
    inset.set_axisbelow(True)
    # Opaque, fully framed: it is a separate panel, not an overlay.
    inset.patch.set_facecolor("white")
    inset.patch.set_alpha(1.0)
    for side in ("top", "right"):
        inset.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        inset.spines[side].set_visible(True)
        inset.spines[side].set_linewidth(0.9)
        inset.spines[side].set_color("black")

    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=ramp), ax=ax, pad=0.015, aspect=28)
    cbar.set_label("DOR fold reduction", fontsize=12)
    ticks = [1, 3, 10, 30, 100]
    cbar.set_ticks(ticks); cbar.set_ticklabels([str(t) for t in ticks])
    cbar.ax.tick_params(labelsize=11)
    # Match the truncated ramp used for the lines.
    cbar.solids.set_alpha(1.0)

    n_mut = len(mutants)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {output}  ({n_mut} mutants + WT)")
    trunc = {m: round(v[4], 1) for m, v in traces.items() if v[4] < 95}
    if trunc:
        print(f"truncated below 95 ns: {trunc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
