#!/usr/bin/env python3
"""Torsion trajectories and distributions, DOR free in water vs bound to RT.

The headline is kinetic, not thermodynamic. Binding does not merely narrow the
torsional distributions -- it abolishes rotamer exchange. Across 60 bound
replicates x 4 torsions x 100 ns, DOR leaves its bound rotamer once. Free in
water it leaves 2-50 times per torsion per 100 ns.

Two traps this script is built around:

  * Frame spacing is not matched. The solvent runs write 2000 frames per 100 ns
    (0.05 ns); the bound runs write 139-540 (0.19-0.72 ns), recorded per
    replicate as ``analysis_dt_ps``. Coarse sampling can only undercount
    transitions, so free trajectories are subsampled to the bound resolution
    before any rate is compared. The figure states the native spacing of each.

  * A fixed 120 degree sector binning puts a well edge at -60, which is 2.3
    degrees from where the bound tau1 distribution is centred, so rattling
    about its own mean registered as rotamer hopping. Occupancy is therefore
    defined as a basin centred on the bound circular mean, which has no edge
    near the populated region.

  * An angle is circular, so a plot of it always has a branch cut somewhere and
    the only question is where. Plotting the deviation from the bound value put
    the cut 180 degrees from the bound rotamer -- which is exactly where tau3's
    free rotamer sits, giving 594 seam crossings per 100 ns and a solid block of
    ink. Measured over the free trajectories, crossings per 100 ns are:

        frame          tau1  tau2  tau3  tau4
        -180..180         2   342     5   860
        0..360           53   286     0    75
        deviation         6   220   594   258
        emptiest cut      2    34     5    75

    The trajectories sidestep the cut entirely by unrolling the angle, so a
    trace is never broken and the plateau-and-jump structure is visible. The
    price is drift: free tau4 makes ~14 net turns and spans 5600 degrees. That
    is why the two columns of a row SHARE a y-range. Letting each autoscale
    draws the bound trace's ~50 degrees of jitter at the same visual amplitude
    as the free trace's thousands, so the bound ligand looks as mobile as the
    free one -- the exact opposite of the result. ``--trace cylinder`` keeps
    continuity on a one-period axis instead; ``--trace broken`` wraps and cuts.
    The ``--angle-frame`` cut still sets the period shown by those two.

Time axes come from each run's recorded ``analysis_dt_ps``, cross-checked
against ``production_ns_completed``; never from the DCD header or frame index.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_dor_torsion_kinetics
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TAUS = ["tau1", "tau2", "tau3", "tau4"]
TAU_LABELS = {
    "tau1": r"$\tau_1$  pyridinone–O ether",
    "tau2": r"$\tau_2$  O–aryl",
    "tau3": r"$\tau_3$  pyridinone-N–CH$_2$",
    "tau4": r"$\tau_4$  CH$_2$–triazolinone",
}

#: Basin half-width about the bound circular mean, degrees. Wide enough that
#: ordinary libration never leaves it (bound circular SDs are 8-12 degrees) and
#: narrow enough that a genuine rotamer change always does.
BASIN_HALF_WIDTH = 60.0

#: Frames a state must persist before a change of state is accepted.
MIN_RUN_FRAMES = 3

FREE = "dor_solvent"
BOUND_REF = "wt"

C_FREE, C_BOUND = "#C44E52", "#4C72B0"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


def circ_mean_deg(a) -> float:
    r = np.deg2rad(np.asarray(a, float)); r = r[np.isfinite(r)]
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circ_dist_deg(a, centre) -> np.ndarray:
    return np.abs((np.asarray(a, float) - centre + 180.0) % 360.0 - 180.0)




def _branch_cut(frame: str, angles_for_scan: np.ndarray | None = None) -> float:
    """Lower edge of the plotted range, in degrees.

    The plotted value is ``(angle - cut) % 360 + cut``, so the seam falls at
    ``cut``. "auto" scans for the angle with the least density within +/-15
    degrees, pooling free and bound so the cut is empty for both.
    """
    if frame == "positive":
        return 0.0
    if frame == "signed":
        return -180.0
    cands = np.arange(-180.0, 180.0, 5.0)
    dens = [float((np.abs((angles_for_scan - c + 180.0) % 360.0 - 180.0) <= 15.0).mean())
            for c in cands]
    return float(cands[int(np.argmin(dens))])


def _wrap(a: np.ndarray, cut: float) -> np.ndarray:
    return (np.asarray(a, float) - cut) % 360.0 + cut




def _cylinder_copies(t: np.ndarray, ang_deg: np.ndarray, lo: float, hi: float
                     ) -> list[tuple[np.ndarray, np.ndarray]]:
    """Continuous copies of a trace, offset by whole turns, that cross [lo, hi].

    A dihedral is periodic, so a plot of one has a branch cut and a trace that
    crosses it must either break or drift. Unwrapping removes the break but not
    the drift: free tau4 makes ~14 net turns in 100 ns and spans 5600 degrees,
    which flattens every panel. Drawing the unwrapped trace once per turn, on an
    axis one period tall, keeps both -- the line leaves the top and re-enters at
    the bottom, which is the correct reading of a cylinder cut open, and no
    segment is ever artificially broken.
    """
    u = np.rad2deg(np.unwrap(np.deg2rad(np.asarray(ang_deg, float))))
    kmin = int(np.floor((lo - u.max()) / 360.0))
    kmax = int(np.ceil((hi - u.min()) / 360.0))
    out = []
    for k in range(kmin, kmax + 1):
        y = u + 360.0 * k
        if y.max() < lo or y.min() > hi:
            continue
        out.append((t, y))
    return out


def _basin_bands(centre: float, half: float, lo: float,
                 hi: float) -> list[tuple[float, float]]:
    """The basin as 1 or 2 spans, splitting it if it runs off the plotted range."""
    a, b = centre - half, centre + half
    if a < lo:
        return [(lo, b), (a + 360.0, hi)]
    if b > hi:
        return [(a, hi), (lo, b - 360.0)]
    return [(a, b)]


def _break_at_seam(t: np.ndarray, y: np.ndarray,
                   jump_deg: float = 180.0) -> tuple[np.ndarray, np.ndarray]:
    """Insert NaN where consecutive angles wrap, so no false vertical stroke."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    brk = np.flatnonzero(np.abs(np.diff(y)) > jump_deg) + 1
    if brk.size == 0:
        return t, y
    return np.insert(t, brk, np.nan), np.insert(y, brk, np.nan)


def count_escapes(angles, centre) -> tuple[int, float]:
    """Returns to the basin after a departure, and the fraction of time inside."""
    inb = circ_dist_deg(angles, centre) <= BASIN_HALF_WIDTH
    edges = np.concatenate(([0], np.flatnonzero(np.diff(inb.astype(int))) + 1,
                            [len(inb)]))
    segs = [(bool(inb[edges[i]]), edges[i + 1] - edges[i])
            for i in range(len(edges) - 1)]
    kept = [s for s, n in segs if n >= MIN_RUN_FRAMES]
    n = sum(1 for i in range(1, len(kept)) if not kept[i] and kept[i - 1])
    return n, float(inb.mean())


def frame_dt_ns(runs_root: Path, label: str, rep: int) -> float | None:
    """Frame spacing from the run JSON, cross-checked against its own length."""
    for cand in sorted((runs_root / label).glob(f"rep_{rep:02d}/*.json")):
        if cand.name.endswith("_md.state.json"):
            continue
        try:
            d = json.loads(cand.read_text())
        except Exception:
            continue
        dt_ps, n = d.get("analysis_dt_ps"), d.get("analysis_n_frames")
        total = d.get("production_ns_completed")
        if not (dt_ps and n and total):
            continue
        if abs(dt_ps * n / 1000.0 - total) > 0.05 * total:
            raise ValueError(
                f"{cand}: analysis_dt_ps x analysis_n_frames = "
                f"{dt_ps * n / 1000.0:.1f} ns but production_ns_completed = "
                f"{total:.1f} ns. Fix the frame schedule before plotting a time axis.")
        return dt_ps / 1000.0
    return None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--torsions-csv", type=Path,
                   default=Path("results/analysis/dor_torsions_solution/tables/"
                                "dor_torsions_per_frame.csv"))
    p.add_argument("--runs-root", type=Path, default=Path("results/md_runs"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/dor_torsions_solution/plots"))
    p.add_argument("--angle-frame", choices=("positive", "signed", "auto"),
                   default="positive",
                   help="Lower edge of the plotted period. positive: 0-360 "
                        "(default). signed: -180..180. auto: per-torsion cut at "
                        "the emptiest angle.")
    p.add_argument("--trace", choices=("broken", "cylinder"), default="broken",
                   help="broken (default): wrapped to one period, line cut at "
                        "the seam. cylinder: unwrapped but drawn once per turn "
                        "on the same one-period axis, so a trace leaves the top "
                        "and re-enters at the bottom instead of breaking.")
    args = p.parse_args(argv)

    plt = _style()
    d = pd.read_csv(args.torsions_csv)
    sol = d[d.safe_label == FREE].sort_values(["replicate", "frame"])
    bound_all = d[d.safe_label != FREE]
    wt = bound_all[bound_all.safe_label == BOUND_REF].sort_values(["replicate", "frame"])
    ref = {t: circ_mean_deg(wt[t]) for t in TAUS}

    dt_free = frame_dt_ns(args.runs_root, "dor_solvent", 2) or 0.05
    dt_wt = {r: frame_dt_ns(args.runs_root, "wt", r) for r in sorted(wt.replicate.unique())}
    print(f"frame spacing: free {dt_free * 1000:.0f} ps; "
          f"WT bound " + ", ".join(f"rep{r} {v * 1000:.0f} ps" for r, v in dt_wt.items()))

    # Branch cut per torsion, and the seam-crossing count it costs, so the
    # choice is recorded in the output rather than argued about.
    cuts, seam = {}, {}
    for t in TAUS:
        pool = np.concatenate([sol[t].to_numpy(), bound_all[t].to_numpy()])
        cuts[t] = _branch_cut(args.angle_frame, pool)
        seam[t] = int(np.mean([
            np.sum(np.abs(np.diff(_wrap(g.sort_values("frame")[t].to_numpy(),
                                        cuts[t]))) > 180.0)
            for _, g in sol.groupby("replicate")]))
    print(f"angle frame {args.angle_frame!r}: " + ", ".join(
        f"{t} cut {cuts[t]:.0f}° ({seam[t]} seam crossings/100 ns free)" for t in TAUS))

    # ---------------- Figure 1: trajectories, each with its projection -------
    # Four column pairs: a trace and, immediately right of it on a shared y
    # axis, the same data projected onto the angle axis. The marginal is the
    # trace's own histogram, so what the eye integrates vertically is exactly
    # what the bars show.
    fig, axes = plt.subplots(
        4, 4, figsize=(8.0, 7.6), sharex="col",
        gridspec_kw={"width_ratios": [3.0, 0.85, 3.0, 0.85], "wspace": 0.08})
    for i, t in enumerate(TAUS):
        c, cut = ref[t], cuts[t]
        lo, hi = cut, cut + 360.0
        c_plot = float(_wrap(np.array([c]), cut)[0])
        bands = _basin_bands(c_plot, BASIN_HALF_WIDTH, lo, hi)
        bins = np.arange(lo, hi + 0.1, 7.5)
        peak = {}
        for j, (which, g_all, colour, dts) in enumerate([
                ("DOR free in water", sol, C_FREE,
                 {r: dt_free for r in sol.replicate.unique()}),
                ("DOR bound to WT RT", wt, C_BOUND, dt_wt)]):
            ax, axm = axes[i, 2 * j], axes[i, 2 * j + 1]
            for a in (ax, axm):
                # The basin is a band about the bound value; where it runs off
                # an edge the remainder is drawn at the other, because the axis
                # is a circle cut open at `cut`.
                for b0, b1 in bands:
                    a.axhspan(b0, b1, color=colour, alpha=0.08, lw=0)
                a.axhline(c_plot, color="0.4", lw=0.5, ls="--", zorder=1)

            pooled = []
            for rep, g in g_all.groupby("replicate"):
                g = g.sort_values("frame")
                dt = dts.get(rep) or dt_free
                time_ns = (np.arange(len(g)) + 1) * dt
                ang = _wrap(g[t].to_numpy(), cut)
                pooled.append(ang)
                if args.trace == "cylinder":
                    segs = _cylinder_copies(time_ns, ang, lo, hi)
                else:
                    segs = [_break_at_seam(time_ns, ang)]
                for xs, ys in segs:
                    ax.plot(xs, ys, lw=0.35, alpha=0.75, color=colour,
                            rasterized=True)

            dens, _, _ = axm.hist(np.concatenate(pooled), bins=bins, density=True,
                                   orientation="horizontal", color=colour,
                                   alpha=0.75, lw=0)
            peak[j] = float(np.max(dens))

            for a in (ax, axm):
                a.set_ylim(lo, hi)
            ax.set_yticks(np.arange(lo, hi + 1, 90))
            ax.set_xlim(0, 100)
            axm.set_yticks(np.arange(lo, hi + 1, 90))
            axm.set_yticklabels([])
            axm.set_xticks([])
            axm.spines["bottom"].set_visible(False)
            if i == 0:
                ax.set_title(which, color=colour, loc="left")
            if j == 0:
                ax.set_ylabel(f"{TAU_LABELS[t]}\n(°)", fontsize=7)
            if i == 3:
                ax.set_xlabel("time (ns)")
                axm.set_xlabel("density", fontsize=6)
        # One density scale for the row. Letting each marginal autoscale draws
        # the bound peak and the much flatter free distribution at the same
        # height, which hides the concentration the figure is about.
        xmax = max(peak.values()) * 1.08
        for j in (0, 1):
            axes[i, 2 * j + 1].set_xlim(0, xmax)
    fig.suptitle("Binding abolishes rotamer exchange in doravirine", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.trace == "broken" else f"_{args.trace}"
    f1 = args.output_dir / f"dor_torsion_trajectories{suffix}.png"
    fig.savefig(f1, dpi=300)
    plt.close(fig)

    # ---------------- Figure 2: distributions ----------------
    fig, axes = plt.subplots(1, 4, figsize=(7.6, 2.1), sharey=True)
    for i, t in enumerate(TAUS):
        ax = axes[i]
        c, cut = ref[t], cuts[t]
        lo, hi = cut, cut + 360.0
        c_plot = float(_wrap(np.array([c]), cut)[0])
        bins = np.arange(lo, hi + 0.1, 7.5)
        for which, g, colour in [("free", sol, C_FREE), ("bound (all)", bound_all, C_BOUND)]:
            ax.hist(_wrap(g[t].to_numpy(), cut), bins=bins, density=True,
                    color=colour, alpha=0.55, label=which, lw=0)
        for b0, b1 in _basin_bands(c_plot, BASIN_HALF_WIDTH, lo, hi):
            ax.axvspan(b0, b1, color="0.5", alpha=0.07, lw=0)
        ax.axvline(c_plot, color="0.4", lw=0.5, ls="--")
        ax.set_xlim(lo, hi)
        ax.set_xticks(np.arange(lo, hi + 1, 90))
        ax.set_title(TAU_LABELS[t], fontsize=7)
        ax.set_xlabel("dihedral (°)", fontsize=7)
        if i == 0:
            ax.set_ylabel("density")
            ax.legend(frameon=False, loc="upper left", fontsize=6)
    fig.suptitle("Torsional distributions, free vs bound", y=1.06, fontsize=10)
    fig.tight_layout()
    f2 = args.output_dir / "dor_torsion_distributions.png"
    fig.savefig(f2, dpi=300)
    plt.close(fig)

    # ---------------- the table behind the figures ----------------
    rows = []
    for t in TAUS:
        c = ref[t]
        fr = []
        for rep, g in sol.groupby("replicate"):
            a = g.sort_values("frame")[t].to_numpy()
            idx = np.linspace(0, len(a) - 1, 180).astype(int)
            fr.append(count_escapes(a[idx], c))
        frn = [x[0] for x in fr]
        frfull = [count_escapes(g.sort_values("frame")[t].to_numpy(), c)
                  for _, g in sol.groupby("replicate")]
        bn = [count_escapes(g.sort_values("frame")[t].to_numpy(), c)
              for _, g in bound_all.groupby(["safe_label", "replicate"])]
        rows.append({
            "torsion": t, "bound_centre_deg": round(c, 1),
            "free_escapes_per_100ns_matched": float(np.mean(frn)),
            "free_escapes_per_100ns_full": float(np.mean([x[0] for x in frfull])),
            "bound_escapes_total_60reps": int(np.sum([x[0] for x in bn])),
            "free_time_in_basin": float(np.mean([x[1] for x in frfull])),
            "bound_time_in_basin": float(np.mean([x[1] for x in bn])),
        })
    tab = pd.DataFrame(rows)
    tdir = args.output_dir.parent / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(tdir / "dor_torsion_kinetics.csv", index=False)

    print("\n=== rotamer kinetics ===")
    print(tab.to_string(index=False))
    tot = tab.bound_escapes_total_60reps.sum()
    print(f"\nacross 60 bound replicates x 4 torsions x 100 ns = 24 us of torsion-time:"
          f" {tot} escape event{'s' if tot != 1 else ''}")
    print(f"\nwrote {f1}\n      {f2}\n      {tdir / 'dor_torsion_kinetics.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
