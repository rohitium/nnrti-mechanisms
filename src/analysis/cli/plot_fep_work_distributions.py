#!/usr/bin/env python3
"""Combined pmx non-equilibrium work-distribution figure for every leg.

One panel per alchemical leg and phase, showing the forward and reverse
switching-work histograms pooled over the three replicates, with each
replicate's BAR free energy marked. The forward/reverse overlap is what
justifies BAR: well-overlapped, similarly-shaped distributions mean the estimate
is limited by statistics, whereas widely separated ones mean dissipation is
large and BAR is biased as well as noisy.

pmx writes the reverse work already sign-flipped, so both directions share an
axis and dG falls between the two means; their separation is the hysteresis.
Work values are stored in kJ/mol and converted here to kcal/mol.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

KJ_PER_KCAL = 4.184
FWD_COLOUR = "#2166AC"
REV_COLOUR = "#B2182B"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--legs-dir", type=Path, default=root / "results/analysis/fep_pmx/legs")
    p.add_argument("--output", type=Path,
                   default=root / "results/analysis/fep_pmx/work_distributions_all_legs.png")
    p.add_argument("--ncols", type=int, default=4, help="Panels per row; keep even so holo/apo stay paired.")
    p.add_argument("--bins", type=int, default=26)
    p.add_argument("--dpi", type=int, default=250)
    return p.parse_args()


def collect(legs_dir: Path) -> dict[tuple[str, str], dict]:
    """(leg, phase) -> pooled forward/reverse work in kcal/mol plus per-replicate BAR."""
    out: dict[tuple[str, str], dict] = {}
    for analysis in sorted(legs_dir.glob("*/*/rep_*/neq/analysis")):
        meta = analysis / "analysis.json"
        fwd, rev = analysis / "integ_fwd.dat", analysis / "integ_rev.dat"
        if not (meta.exists() and fwd.exists() and rev.exists()):
            continue
        payload = json.loads(meta.read_text())
        key = (str(payload["leg_id"]), str(payload["phase"]))
        entry = out.setdefault(key, {"fwd": [], "rev": [], "bar": [], "reps": []})
        entry["fwd"].append(np.loadtxt(fwd, usecols=1) / KJ_PER_KCAL)
        entry["rev"].append(np.loadtxt(rev, usecols=1) / KJ_PER_KCAL)
        if payload.get("bar_dg") is not None:
            entry["bar"].append(float(payload["bar_dg"]))
            entry["reps"].append(int(payload["replicate"]))
    for entry in out.values():
        entry["fwd"] = np.concatenate(entry["fwd"])
        entry["rev"] = np.concatenate(entry["rev"])
    return out


def leg_order(keys) -> list[tuple[str, str]]:
    """Primary legs from WT first, then the secondary legs; holo before apo."""
    legs = sorted({leg for leg, _phase in keys})
    primary = [l for l in legs if l.startswith("wt_to_")]
    secondary = [l for l in legs if not l.startswith("wt_to_")]
    ordered = []
    for leg in [*primary, *secondary]:
        for phase in ("holo", "apo"):
            if (leg, phase) in keys:
                ordered.append((leg, phase))
    return ordered


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    args = parse_args()
    data = collect(args.legs_dir)
    if not data:
        raise SystemExit(f"no work data under {args.legs_dir}")
    panels = leg_order(set(data))

    ncols = max(2, args.ncols)
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.05 * ncols, 1.75 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, key in zip(axes, panels):
        leg, phase = key
        e = data[key]
        lo = min(e["fwd"].min(), e["rev"].min())
        hi = max(e["fwd"].max(), e["rev"].max())
        bins = np.linspace(lo, hi, args.bins)
        ax.hist(e["fwd"], bins=bins, color=FWD_COLOUR, alpha=0.55, density=True, lw=0)
        ax.hist(e["rev"], bins=bins, color=REV_COLOUR, alpha=0.55, density=True, lw=0)
        for dg in e["bar"]:
            ax.axvline(dg, color="0.25", lw=0.9, ls="--", zorder=4)
        hyst = float(e["fwd"].mean() - e["rev"].mean())
        # Standardised separation of the two work distributions. STATUS.md treats
        # ~4-5 sigma as the regime where BAR is biased rather than merely noisy,
        # so the annotation is coloured on that scale.
        pooled = float(np.sqrt((e["fwd"].var(ddof=1) + e["rev"].var(ddof=1)) / 2.0))
        sep = hyst / pooled if pooled > 0 else float("nan")
        colour = "0.35" if sep < 2 else ("#B35806" if sep < 3.5 else "#B2182B")
        label = leg.replace("wt_to_", "WT→").replace("_to_", "→").replace("_", "+")
        ax.set_title(f"{label}  ({phase})", fontsize=8.5, pad=2.5)
        ax.text(0.03, 0.93, f"{hyst:.1f} kcal · {sep:.1f}σ", transform=ax.transAxes,
                fontsize=7, va="top", color=colour,
                fontweight="bold" if sep >= 3.5 else "normal")
        ax.tick_params(labelsize=7, length=2.5, pad=1.5)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
    for ax in axes[len(panels):]:
        ax.set_visible(False)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=FWD_COLOUR, alpha=0.55),
        plt.Rectangle((0, 0), 1, 1, color=REV_COLOUR, alpha=0.55),
        Line2D([0], [0], color="0.25", lw=0.9, ls="--"),
    ]
    fig.legend(handles, ["forward (0→1)", "reverse (1→0)", "BAR ΔG, per replicate"],
               loc="upper center", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 1.0))
    fig.text(0.5, 0.9955,
             "Each panel is annotated with the forward−reverse hysteresis and its separation in "
             "pooled standard deviations; amber ≥ 2σ, red ≥ 3.5σ, where BAR becomes biased "
             "rather than merely noisy.",
             ha="center", fontsize=8.5, color="0.3")
    fig.supxlabel("Switching work (kcal/mol)", fontsize=12, y=0.004)
    fig.suptitle("Non-equilibrium switching work distributions, all alchemical legs",
                 fontsize=13, y=1.013)
    fig.tight_layout(rect=(0, 0.012, 1, 0.985))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.output}  ({len(panels)} panels, "
          f"{len({l for l, _p in panels})} legs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
