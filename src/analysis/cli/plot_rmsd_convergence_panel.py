#!/usr/bin/env python3
"""Single-panel Calpha RMSD convergence figure for the whole genotype panel.

One mean trace per genotype with a shaded across-replicate SEM band, coloured by
DOR fold reduction so the reader can see that structural drift carries no
resistance ordering. Fold values come from the authoritative susceptibility
spreadsheet, never from the manifest.

Replicates run to slightly different lengths, so each genotype's trace stops at
its shortest replicate: the mean is always over n = 3, never a changing
denominator.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


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
    p.add_argument("--output", type=Path,
                   default=root / "results/analysis/md_convergence/plots/rmsd_convergence.png")
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


def genotype_traces(rmsd: pd.DataFrame, grid_points: int):
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
            np.interp(grid, r["time_ns"].to_numpy(), r["ca_rmsd_angstrom"].to_numpy())
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
    rmsd = pd.read_csv(args.rmsd_csv)
    if args.exclude:
        rmsd = rmsd[~rmsd["mutation"].isin(args.exclude)]
        print(f"excluded: {', '.join(args.exclude)}")
    folds = load_fold_map(args.susceptibility_xlsx)
    traces = genotype_traces(rmsd, args.grid_points)

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
    ax.set_ylabel(r"C$\alpha$ RMSD ($\mathrm{\AA}$)", fontsize=14)
    ax.set_title("RT–DOR backbone convergence across the genotype panel", fontsize=15, pad=12)
    ax.tick_params(labelsize=12)
    ax.set_xlim(0, t_max)
    lo = min(float((m - e).min()) for _g, m, e, _n, _t in traces.values())
    hi = max(float((m + e).max()) for _g, m, e, _n, _t in traces.values())
    pad = 0.08 * (hi - lo)
    ax.set_ylim(max(0.0, lo - pad), hi + pad)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=11, frameon=False, loc="upper left")

    # Inset: the equilibration rise, which is compressed against the y-axis in the
    # main panel. Same series, same colours, no legend -- it is a magnification.
    inset = ax.inset_axes([0.575, 0.105, 0.40, 0.295])
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
    inset.set_ylabel(r"C$\alpha$ RMSD ($\mathrm{\AA}$)", fontsize=8.5, labelpad=1.0)
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.output}  ({n_mut} mutants + WT)")
    trunc = {m: round(v[4], 1) for m, v in traces.items() if v[4] < 95}
    if trunc:
        print(f"truncated below 95 ns: {trunc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
