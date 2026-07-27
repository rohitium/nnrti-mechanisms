#!/usr/bin/env python3
"""QC for NEQ switch work distributions (PLAN.md §4.6 gates).

Per leg/phase/replicate, reads the per-switch work values dumped by
``analyze_neq.py`` and reports:
  - forward/reverse Crooks overlap coefficient (P(W_f) vs P(-W_r))
  - work outlier fraction (MAD-based) per direction
  - BAR vs Jarzynski agreement

Outputs (under ``results/analysis/fep_pmx/``):
  panel_qc.csv                one row per leg/phase/replicate, with flags
  panel_crooks_overlap.png    forward vs reverse work histograms, per leg/phase
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.analyze_neq import ensure_leg_analysis, read_work_values_kcal
from scripts.fep_pmx.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K, P0_LEGS

# PLAN.md §4.6 pass thresholds
OVERLAP_MIN = 0.3
OUTLIER_FRAC_MAX = 0.05
BAR_JARZ_MAX_KCAL = 1.0


def outlier_fraction(work: np.ndarray) -> float:
    """MAD-based outlier fraction: |w - median| > 3 · 1.4826 · MAD."""
    if len(work) < 3:
        return 0.0
    med = np.median(work)
    mad = np.median(np.abs(work - med))
    if mad == 0:
        return 0.0
    return float(np.mean(np.abs(work - med) > 3.0 * 1.4826 * mad))


def overlap_coefficient(wf: np.ndarray, neg_wr: np.ndarray, nbins: int = 30) -> float:
    """Histogram overlap coefficient between P(W_f) and P(-W_r), in [0, 1]."""
    if len(wf) < 2 or len(neg_wr) < 2:
        return 0.0
    lo = min(wf.min(), neg_wr.min())
    hi = max(wf.max(), neg_wr.max())
    if hi <= lo:
        return 0.0
    bins = np.linspace(lo, hi, nbins + 1)
    hf, _ = np.histogram(wf, bins=bins, density=True)
    hr, _ = np.histogram(neg_wr, bins=bins, density=True)
    width = bins[1] - bins[0]
    return float(np.sum(np.minimum(hf, hr)) * width)


def qc_unit(leg_id: str, phase: str, replicate: int, *, temperature_k: float, nboots: int, auto: bool) -> dict:
    meta = ensure_leg_analysis(
        leg_id, phase=phase, replicate=replicate,
        temperature_k=temperature_k, nboots=nboots, auto=auto,
    )
    wf = np.array(read_work_values_kcal(Path(meta["integ_fwd"])))
    wr = np.array(read_work_values_kcal(Path(meta["integ_rev"])))
    neg_wr = -wr

    overlap = overlap_coefficient(wf, neg_wr)
    of_fwd = outlier_fraction(wf)
    of_rev = outlier_fraction(wr)
    bar_dg = meta.get("bar_dg")
    jarz = meta.get("jarz_dg_mean")
    bar_minus_jarz = (bar_dg - jarz) if (bar_dg is not None and jarz is not None) else None

    flags = []
    if overlap < OVERLAP_MIN:
        flags.append("low_overlap")
    if max(of_fwd, of_rev) > OUTLIER_FRAC_MAX:
        flags.append("outliers")
    if bar_minus_jarz is not None and abs(bar_minus_jarz) > BAR_JARZ_MAX_KCAL:
        flags.append("bar_jarz_disagree")

    return {
        "leg_id": leg_id, "phase": phase, "replicate": replicate,
        "n_fwd": len(wf), "n_rev": len(wr),
        "wf_mean": float(wf.mean()) if len(wf) else None,
        "wr_mean": float(wr.mean()) if len(wr) else None,
        "overlap": overlap, "outlier_frac_fwd": of_fwd, "outlier_frac_rev": of_rev,
        "bar_dg": bar_dg, "jarz_dg_mean": jarz, "bar_minus_jarz": bar_minus_jarz,
        "flags": ";".join(flags) if flags else "ok",
        "_wf": wf, "_neg_wr": neg_wr,
    }


def _plot(units: list[dict], output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # group by (leg, phase): pool reps
    groups: dict[tuple[str, str], list[dict]] = {}
    for u in units:
        groups.setdefault((u["leg_id"], u["phase"]), []).append(u)
    if not groups:
        return

    keys = sorted(groups)
    ncols = 2
    nrows = (len(keys) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 3.2 * nrows), squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)

    for idx, key in enumerate(keys):
        ax = axes[idx // ncols][idx % ncols]
        ax.set_visible(True)
        us = groups[key]
        wf = np.concatenate([u["_wf"] for u in us]) if us else np.array([])
        neg_wr = np.concatenate([u["_neg_wr"] for u in us]) if us else np.array([])
        if len(wf):
            ax.hist(wf, bins=25, alpha=0.55, density=True, label="forward $W_f$", color="#2c6fbb")
        if len(neg_wr):
            ax.hist(neg_wr, bins=25, alpha=0.55, density=True, label="reverse $-W_r$", color="#d1642f")
        bar_vals = [u["bar_dg"] for u in us if u["bar_dg"] is not None]
        if bar_vals:
            ax.axvline(float(np.mean(bar_vals)), color="0.2", lw=1.2, ls="--",
                       label=f"BAR ΔG = {np.mean(bar_vals):.2f}")
        ov = np.mean([u["overlap"] for u in us])
        ax.set_title(f"{key[0]}  {key[1]}  (overlap {ov:.2f}, {len(us)} rep)", fontsize=9)
        ax.set_xlabel("work (kcal/mol)")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)

    fig.suptitle("NEQ Crooks work overlap (forward vs reverse)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="QC pmx NEQ work distributions.")
    parser.add_argument("--legs", nargs="+", default=list(P0_LEGS))
    parser.add_argument("--phases", nargs="+", default=["holo", "apo"])
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    parser.add_argument("--nboots", type=int, default=100)
    parser.add_argument("--no-auto", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=FEP_PMX_ROOT)
    args = parser.parse_args(argv)

    units: list[dict] = []
    for leg_id in args.legs:
        for phase in args.phases:
            for replicate in range(1, args.replicates + 1):
                try:
                    units.append(qc_unit(
                        leg_id, phase, replicate,
                        temperature_k=args.temperature_k, nboots=args.nboots, auto=not args.no_auto,
                    ))
                except (FileNotFoundError, RuntimeError) as exc:
                    print(f"skip {leg_id} {phase} rep{replicate}: {exc}", file=sys.stderr)

    if not units:
        print("No completed NEQ units to QC.", file=sys.stderr)
        return 1

    csv_path = args.output_dir / "panel_qc.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["leg_id", "phase", "replicate", "n_fwd", "n_rev", "wf_mean", "wr_mean",
              "overlap", "outlier_frac_fwd", "outlier_frac_rev",
              "bar_dg", "jarz_dg_mean", "bar_minus_jarz", "flags"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for u in units:
            writer.writerow(u)

    plot_path = args.output_dir / "panel_crooks_overlap.png"
    _plot(units, plot_path)

    print(f"{'leg':<16}{'phase':<6}{'rep':>4}{'overlap':>9}{'out_f':>7}{'out_r':>7}{'BAR-Jarz':>10}  flags")
    for u in units:
        bj = f"{u['bar_minus_jarz']:+.2f}" if u["bar_minus_jarz"] is not None else "   -"
        print(f"{u['leg_id']:<16}{u['phase']:<6}{u['replicate']:>4}{u['overlap']:>9.2f}"
              f"{u['outlier_frac_fwd']:>7.2f}{u['outlier_frac_rev']:>7.2f}{bj:>10}  {u['flags']}")
    n_flagged = sum(1 for u in units if u["flags"] != "ok")
    print(f"\n{len(units)} units, {n_flagged} flagged.\nWrote {csv_path}\nWrote {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
