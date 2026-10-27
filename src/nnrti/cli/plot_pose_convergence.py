#!/usr/bin/env python3
"""Does a perturbed doravirine pose fall into the wild-type basin?

Reviewer 1 (comment 2) asks whether the same DOR-bound state is recovered from
independent starting poses. A start placed far enough away is simply trapped,
which answers nothing; the informative experiment is a ladder of starts at
increasing distance from the crystallographic pose, read for where convergence
stops.

This draws that ladder on the two pocket-frame coordinates that separate the
two ways a rigid ligand can move: centroid displacement (translation) and
long-axis tilt (rotation). Behind them sits the region the wild-type
crystal-seeded replicates actually visit. Each ladder run is drawn as an arrow
from where it started to where it ended, so a pose that converges is an arrow
landing inside that region and one that does not is an arrow that stays out.

Starting coordinates are the pose **as generated**, before minimization. What
matters is that a pose at a given distance converged, not whether the
minimizer or the dynamics did the work.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_pose_convergence
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]

DENSITY = "#2c7fb8"
CRYSTAL = "#1a1a1a"
LANDED = "#2b8a3e"      # converged into the wild-type region
STUCK = "#c92a2a"       # did not
TAIL = "#9e9e9e"


def density_contours(ax, x, y, fractions=(0.50, 0.90, 0.99)):
    """Nested highest-density contours enclosing the given fractions of frames."""
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(np.vstack([x, y]))
    at_samples = kde(np.vstack([x, y]))
    levels = sorted({float(np.quantile(at_samples, 1.0 - f)) for f in fractions})

    pad_x = max(0.35, 0.18 * (x.max() - x.min()))
    pad_y = max(6.0, 0.18 * (y.max() - y.min()))
    gx, gy = np.mgrid[
        max(0.0, x.min() - pad_x) : x.max() + pad_x : 240j,
        max(0.0, y.min() - pad_y) : y.max() + pad_y : 240j,
    ]
    grid = kde(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    ax.contourf(gx, gy, grid, levels=levels + [grid.max()],
                colors=DENSITY, alpha=0.18, zorder=1)
    cs = ax.contour(gx, gy, grid, levels=levels, colors=DENSITY,
                    linewidths=1.2, zorder=2)
    ax.clabel(cs, fmt={lvl: f"{int(f*100)}%" for lvl, f in
                       zip(levels, sorted(fractions, reverse=True))},
              fontsize=8, inline=True)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--wt-frames", type=Path,
                    default=REPO_ROOT / "results/analysis/pose_coordinates/tables/wt/sampled_frames.csv")
    ap.add_argument("--ladder", type=Path,
                    default=REPO_ROOT / "results/analysis/pose_coordinates/tables/ladder_summary.csv",
                    help="One row per ladder run: pose, start_rmsd_angstrom, "
                         "start_displacement, start_tilt, end_displacement, end_tilt, end_rmsd.")
    ap.add_argument("--basin-rmsd", type=float, default=1.78,
                    help="Upper edge of the wild-type sampled pose RMSD.")
    ap.add_argument("--output", type=Path,
                    default=REPO_ROOT / "results/analysis/pose_coordinates/pose_convergence.png")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args(argv)

    wt = pd.read_csv(args.wt_frames)
    if "min_contact_angstrom" in wt.columns:
        wt = wt[wt["min_contact_angstrom"] >= 2.2]
    lad = pd.read_csv(args.ladder)

    fig, ax = plt.subplots(figsize=(8.8, 7.0))

    x = wt["centroid_displacement_angstrom"].to_numpy()
    y = wt["long_axis_tilt_deg"].to_numpy()
    density_contours(ax, x, y)
    ax.scatter(x, y, s=4, color=DENSITY, alpha=0.22, linewidths=0, zorder=3,
               label=f"Wild-type frames (3 × 100 ns, n={len(wt)})")

    # the crystallographic pose: the origin of these coordinates by construction
    ax.scatter([0.0], [0.0], marker="*", s=520, color=CRYSTAL,
               edgecolor="white", linewidth=1.5, zorder=9,
               label="Crystallographic pose (4NCG)")

    for _, r in lad.sort_values("start_rmsd_angstrom").iterrows():
        landed = bool(r["end_rmsd_angstrom"] <= args.basin_rmsd)
        colour = LANDED if landed else STUCK
        ax.annotate(
            "", xy=(r["end_displacement_angstrom"], r["end_tilt_deg"]),
            xytext=(r["start_displacement_angstrom"], r["start_tilt_deg"]),
            arrowprops=dict(arrowstyle="-|>", color=TAIL, linewidth=1.1,
                            shrinkA=4, shrinkB=5, alpha=0.9), zorder=4,
        )
        ax.scatter([r["start_displacement_angstrom"]], [r["start_tilt_deg"]],
                   marker="o", s=70, facecolor="white", edgecolor=colour,
                   linewidth=1.8, zorder=6)
        ax.scatter([r["end_displacement_angstrom"]], [r["end_tilt_deg"]],
                   marker="o", s=95, color=colour, edgecolor="white",
                   linewidth=1.2, zorder=7)
        ax.annotate(f"{r['start_rmsd_angstrom']:.1f}",
                    xy=(r["start_displacement_angstrom"], r["start_tilt_deg"]),
                    xytext=(5, 6), textcoords="offset points",
                    fontsize=8, color=colour, zorder=8)

    # legend proxies
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles += [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor=CRYSTAL, markeredgewidth=1.6, markersize=9,
               label="Start (label: pose RMSD as generated, Å)"),
        Line2D([], [], marker="o", linestyle="none", color=LANDED, markersize=9,
               label=f"End, inside the basin (RMSD ≤ {args.basin_rmsd:.2f} Å)"),
        Line2D([], [], marker="o", linestyle="none", color=STUCK, markersize=9,
               label="End, outside"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, framealpha=0.95)

    ax.set_xlabel("DOR centroid displacement from the crystallographic pose (Å)")
    ax.set_ylabel("DOR long-axis tilt from the crystallographic pose (°)")
    ax.set_title("Perturbed doravirine poses against the wild-type sampled basin")
    ax.margins(0.06)
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    fig.savefig(args.output.with_suffix(".pdf"))
    n_landed = int((lad["end_rmsd_angstrom"] <= args.basin_rmsd).sum())
    logging.info("%d of %d ladder runs ended inside the basin", n_landed, len(lad))
    logging.info("wrote %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
