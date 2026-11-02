#!/usr/bin/env python3
"""Where the perturbed starting poses sit relative to the sampled ensemble.

Reviewer 1 (comment 2) asks whether the same DOR-bound state is recovered from
independent starting poses.  Answering that needs two things on one pair of axes: the
region the crystal-seeded simulations actually visit, and the starting poses we intend
to run.  If a start falls inside the sampled region it is not independent of what we
already have, whatever its nominal RMSD.

Axes are the pocket-frame pose coordinates from compute_pose_coordinates: centroid
displacement separates translation, long-axis tilt separates rotation.  The published
coordinates cannot be used here - a rotation about the ligand centroid leaves both of
them unchanged, so every rotated start would plot on top of the crystal point.

The replicates are written at different DCD strides (500 frames against 180 for the
same 100 ns), so frames are subsampled to the shortest replicate before the density is
drawn.  Otherwise one replicate carries most of the weight and the density describes it
rather than the ensemble.

The right-hand panel carries pose RMSD, the single scalar, for the same comparison.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_pose_landscape
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Categorical hues follow the repo convention; validated for CVD separation together.
FAMILY_STYLE = {
    "rotation": {"color": "#d7301f", "marker": "o", "label": "Rotated start"},
    "alternative": {"color": "#756bb1", "marker": "D", "label": "Alternative placement"},
}
CRYSTAL_COLOR = "#1a1a1a"
DENSITY_COLOR = "#2c7fb8"
LABEL_INK = "#404040"


def _balance_replicates(frames: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """Equal frames per replicate, so the density is the ensemble's and not one run's."""
    smallest = frames.groupby("replicate").size().min()
    rng = np.random.default_rng(seed)
    return pd.concat(
        [
            group.iloc[np.sort(rng.choice(len(group), size=smallest, replace=False))]
            for _, group in frames.groupby("replicate")
        ]
    )


def _density_contours(ax, x, y, fractions=(0.50, 0.90, 0.99)):
    """Nested highest-density contours enclosing the given fractions of sampled frames."""
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(np.vstack([x, y]))
    at_samples = kde(np.vstack([x, y]))
    # The level enclosing fraction f is the (1-f) quantile of the density at the samples.
    levels = sorted({float(np.quantile(at_samples, 1.0 - f)) for f in fractions})

    pad_x = 0.12 * (x.max() - x.min() or 1.0)
    pad_y = 0.12 * (y.max() - y.min() or 1.0)
    gx, gy = np.mgrid[
        x.min() - pad_x : x.max() + pad_x : 220j,
        y.min() - pad_y : y.max() + pad_y : 220j,
    ]
    grid = kde(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)

    ax.contourf(gx, gy, grid, levels=levels + [grid.max()], colors=DENSITY_COLOR,
                alpha=0.16, zorder=1)
    contours = ax.contour(gx, gy, grid, levels=levels, colors=DENSITY_COLOR,
                          linewidths=1.3, zorder=2)
    ax.clabel(
        contours,
        fmt={lvl: f"{int(f * 100)}%" for lvl, f in zip(levels, sorted(fractions, reverse=True))},
        fontsize=8,
        inline=True,
    )


def _stacked_labels(ax, points, x_label, min_gap, fontsize=8.0, ha="left"):
    """Annotate points in one column, pushing labels apart so none overlaps.

    points: list of (x, y, text, color), drawn with a leader line back to the mark.
    """
    ordered = sorted(points, key=lambda p: p[1])
    placed: list[float] = []
    for _, y, _, _ in ordered:
        target = y if not placed else max(y, placed[-1] + min_gap)
        placed.append(target)
    # Centre the column on the original spread rather than letting it drift upward.
    shift = (np.mean([p[1] for p in ordered]) - np.mean(placed)) if placed else 0.0
    for (px, py, text, color), label_y in zip(ordered, placed):
        ax.annotate(
            text,
            xy=(px, py),
            xytext=(x_label, label_y + shift),
            textcoords="data",
            fontsize=fontsize,
            color=color,
            va="center",
            ha=ha,
            arrowprops=dict(arrowstyle="-", color="#b0b0b0", linewidth=0.6,
                            shrinkA=0.0, shrinkB=3.0),
            zorder=8,
        )


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    tables = REPO_ROOT / "results" / "analysis" / "pose_coordinates" / "tables"
    parser.add_argument("--frames", type=Path, default=tables / "sampled_frames.csv")
    parser.add_argument("--starts", type=Path, default=tables / "starting_poses.csv")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "analysis" / "pose_coordinates" / "pose_landscape.png",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--selected", nargs="*", default=None,
        help="Pose stems chosen for simulation; drawn filled, the rest faded.",
    )
    parser.add_argument(
        "--min-contact",
        type=float,
        default=2.2,
        help="Drop frames whose closest ligand-protein heavy-atom approach is below this "
             "(A). Such frames are not equilibrated structures; boundness_qc tests only "
             "the 6 A upper bound and passes them.",
    )
    args = parser.parse_args(argv)

    frames = pd.read_csv(args.frames)
    starts = pd.read_csv(args.starts)
    if frames.empty:
        raise SystemExit(f"No sampled frames in {args.frames}")

    if "min_contact_angstrom" in frames.columns and args.min_contact > 0:
        rejected = frames[frames["min_contact_angstrom"] < args.min_contact]
        if not rejected.empty:
            logging.warning(
                "Dropping %d/%d frames with contacts below %.1f A: %s",
                len(rejected), len(frames), args.min_contact,
                ", ".join(
                    f"rep {r:02d} n={n}" for r, n in rejected.groupby("replicate").size().items()
                ),
            )
        frames = frames[frames["min_contact_angstrom"] >= args.min_contact]

    label = frames["mutation"].iloc[0]
    n_reps = frames["replicate"].nunique()
    balanced = _balance_replicates(frames)
    logging.info(
        "%d frames pooled, balanced to %d (%d per replicate)",
        len(frames), len(balanced), len(balanced) // n_reps,
    )

    fig, (ax, ax_rmsd) = plt.subplots(
        1, 2, figsize=(13.4, 6.8), gridspec_kw={"width_ratios": [2.2, 1.0]}
    )

    # ---- main panel: the two pocket-frame coordinates ----
    x = balanced["centroid_displacement_angstrom"].to_numpy()
    y = balanced["long_axis_tilt_deg"].to_numpy()

    _density_contours(ax, x, y)
    ax.scatter(x, y, s=5, color=DENSITY_COLOR, alpha=0.30, linewidths=0, zorder=3,
               label=f"{label} frames ({n_reps} × 100 ns)")

    crystal = starts[starts["family"] == "reference"]
    ax.scatter(
        crystal["centroid_displacement_angstrom"], crystal["long_axis_tilt_deg"],
        marker="*", s=420, color=CRYSTAL_COLOR, edgecolor="white", linewidth=1.4,
        zorder=6, label="Crystallographic pose (4NCG)",
    )

    for family, style in FAMILY_STYLE.items():
        subset = starts[starts["family"] == family]
        if subset.empty:
            continue
        chosen = set(args.selected or [])
        if chosen:
            picked = subset[subset["pose"].isin(chosen)]
            rest = subset[~subset["pose"].isin(chosen)]
            ax.scatter(
                rest["centroid_displacement_angstrom"], rest["long_axis_tilt_deg"],
                marker=style["marker"], s=80, facecolor="none",
                edgecolor=style["color"], linewidth=1.3, alpha=0.55, zorder=4,
                label="Not selected",
            )
            ax.scatter(
                picked["centroid_displacement_angstrom"], picked["long_axis_tilt_deg"],
                marker=style["marker"], s=190, color=style["color"],
                edgecolor="#1a1a1a", linewidth=1.8, zorder=5,
                label="Selected for simulation",
            )
        else:
            ax.scatter(
                subset["centroid_displacement_angstrom"], subset["long_axis_tilt_deg"],
                marker=style["marker"], s=115, color=style["color"],
                edgecolor="white", linewidth=1.5, zorder=5, label=style["label"],
            )

    # Rotations all sit at displacement ~0 and would collide; label them as one column.
    rotations = starts[starts["family"] == "rotation"]
    x_right = max(starts["centroid_displacement_angstrom"].max(), x.max())
    _stacked_labels(
        ax,
        [(r["centroid_displacement_angstrom"], r["long_axis_tilt_deg"], r["pose"],
          FAMILY_STYLE["rotation"]["color"]) for _, r in rotations.iterrows()],
        x_label=0.10 * x_right,
        min_gap=9.0,
    )
    chosen = set(args.selected or [])
    for _, row in starts[starts["family"] == "alternative"].iterrows():
        picked = row["pose"] in chosen
        ax.annotate(
            row["pose"],
            (row["centroid_displacement_angstrom"], row["long_axis_tilt_deg"]),
            textcoords="offset points", xytext=(12, 3),
            fontsize=8.8 if picked else 7.6,
            fontweight="bold" if picked else "normal",
            color="#1a1a1a" if picked else "#9590b8", zorder=7,
        )

    ax.set_xlabel("DOR centroid displacement from crystal pose (Å)", fontsize=11)
    ax.set_ylabel("DOR long-axis tilt from crystal orientation (°)", fontsize=11)
    ax.set_title(
        f"Starting poses against the region sampled by {label}\n"
        "pocket-frame coordinates, NNIBP Cα superposition",
        fontsize=12, loc="left",
    )
    ax.grid(alpha=0.18, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5, loc="upper left")

    # ---- right panel: the single scalar ----
    sampled_rmsd = balanced["pose_rmsd_angstrom"].to_numpy()
    lo, hi = sampled_rmsd.min(), sampled_rmsd.max()
    ax_rmsd.hist(sampled_rmsd, bins=46, color=DENSITY_COLOR, alpha=0.70,
                 orientation="horizontal", zorder=2)
    ax_rmsd.axhspan(lo, hi, color=DENSITY_COLOR, alpha=0.10, zorder=0)

    count_max = np.histogram(sampled_rmsd, bins=46)[0].max()
    ax_rmsd.set_xlim(0, count_max * 2.05)
    for _, row in starts.iterrows():
        color = (CRYSTAL_COLOR if row["family"] == "reference"
                 else FAMILY_STYLE[row["family"]]["color"])
        ax_rmsd.plot([0, count_max * 1.04], [row["pose_rmsd_angstrom"]] * 2,
                     color=color, linewidth=1.3, alpha=0.85, zorder=3)
    _stacked_labels(
        ax_rmsd,
        [(count_max * 1.04, r["pose_rmsd_angstrom"], r["pose"],
          CRYSTAL_COLOR if r["family"] == "reference" else FAMILY_STYLE[r["family"]]["color"])
         for _, r in starts.iterrows()],
        x_label=count_max * 1.20,
        min_gap=0.30,
        fontsize=7.8,
    )

    ax_rmsd.set_ylabel("DOR pose RMSD from crystal pose (Å)", fontsize=11)
    ax_rmsd.set_xlabel(f"{label} frames", fontsize=11)
    ax_rmsd.set_title(f"Sampled range {lo:.2f}–{hi:.2f} Å", fontsize=11, loc="left")
    ax_rmsd.grid(alpha=0.18, linewidth=0.6, axis="y")
    ax_rmsd.set_axisbelow(True)
    for spine in ("top", "right"):
        ax_rmsd.spines[spine].set_visible(False)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    logging.info("Wrote %s", args.output)

    inside = starts[(starts["pose_rmsd_angstrom"] >= lo) & (starts["pose_rmsd_angstrom"] <= hi)]
    logging.info(
        "Starts inside the sampled pose-RMSD range: %d/%d (%s)",
        len(inside), len(starts), ", ".join(inside["pose"]) or "none",
    )
    tilt_max = balanced["long_axis_tilt_deg"].max()
    outside_tilt = starts[starts["long_axis_tilt_deg"] > tilt_max]
    logging.info(
        "Starts beyond the sampled tilt range (%.1f°): %d/%d",
        tilt_max, len(outside_tilt), len(starts),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
