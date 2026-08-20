#!/usr/bin/env python3
"""Panel ΔΔG_bind vs experimental fold, split by clinical resistance category.

Why a subset analysis
---------------------
The full panel mixes genotypes whose clinical phenotype is well established with
genotypes whose phenotypic data is sparse or ambiguous (the ``Uncertain`` class).
The question the panel is meant to answer -- does the computed ΔΔG_bind track
resistance *where resistance is actually known* -- is therefore best asked on the
``Susceptible`` + ``Resistant`` genotypes alone.

This is a **pre-specified** subset, not a post-hoc one: the three category sets are
hardcoded in ``plot_wt_referenced_occupancy_tick_lines`` on clinical-evidence
grounds and predate these FEP results (Table-2-energetics.csv still lists G190E's
ΔΔG as "not determined"). They are imported from there rather than redefined, so
the two analyses cannot drift apart.

Robustness note that belongs beside the result
----------------------------------------------
On the current panel the numerical difference between the full set and the
Susceptible+Resistant subset is dominated by a single genotype, V106M
(ΔΔG +6.10 ± 0.16 at fold 3.4 -- the tightest value in the panel, and the clearest
binding-vs-phenotype case). Dropping V106M alone from the full panel gives R² 0.254,
essentially the subset value; adding it back into the subset collapses R² to 0.087.
That is not a reason to discard either analysis, but it should be reported: the
effect size rests heavily on one point, and a reviewer running leave-one-out will
find it immediately.

Usage
-----
    PYTHONPATH=. python -m src.analysis.cli.plot_panel_by_resistance_category
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .plot_wt_referenced_occupancy_tick_lines import RESISTANT, SUSCEPTIBLE, UNCERTAIN

# Reuse combine_neq's force-repel labeller so this figure places the clustered
# high-fold V106A compounds the same way the main panel does, instead of the
# fixed-offset annotations that overlapped them.
_FEP_PMX_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts" / "fep_pmx"

FEP_ROOT = Path("results/analysis/fep_pmx")
COLORS = {"Susceptible": "#2E7D32", "Resistant": "#C62828", "Uncertain": "#9E9E9E"}
MARKERS = {"Susceptible": "o", "Resistant": "s", "Uncertain": "^"}

# Genotypes whose label the automatic repel places poorly. Offsets are in points
# relative to the marker; these are annotated by hand (right-aligned, with a leader)
# and withheld from _repel_labels so it does not fight the manual placement.
LABEL_OVERRIDES = {
    "K103N+P225H": (-14, -16),   # left and down, clear of G190S
    "V106A+F227L": (-14, 6),
    "Y188L": (-14, 10),
    "V106A+P225H": (-16, -4),
}
AXIS_LABEL_FONTSIZE = 20  # 2x the matplotlib default
TICK_LABEL_FONTSIZE = 15  # scaled to stay proportionate to the axis labels


def _load_repel():
    """Import ``_repel_labels`` from scripts/fep_pmx/combine_neq (not a package)."""
    import importlib.util
    import sys

    path = _FEP_PMX_SCRIPTS / "combine_neq.py"
    sys.path.insert(0, str(_FEP_PMX_SCRIPTS))
    spec = importlib.util.spec_from_file_location("_fep_combine_neq", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod._repel_labels


_repel_labels = _load_repel()


def category_for(genotype: str) -> str:
    if genotype in SUSCEPTIBLE:
        return "Susceptible"
    if genotype in RESISTANT:
        return "Resistant"
    if genotype in UNCERTAIN:
        return "Uncertain"
    raise ValueError(
        f"{genotype!r} is not in SUSCEPTIBLE/RESISTANT/UNCERTAIN. Add it to "
        "plot_wt_referenced_occupancy_tick_lines so both analyses stay consistent."
    )


def load_rows(panel_csv: Path) -> list[dict]:
    rows = []
    for r in csv.DictReader(panel_csv.open()):
        fold = (r.get("dor_fold_reduction") or "").strip()
        if not fold:
            continue  # F227C has no clinical fold
        rows.append(
            {
                "genotype": r["genotype"],
                "fold": float(fold),
                "ddg": float(r["ddg_bind_kcal"]),
                "sem": float(r["sem_kcal"]),
                "category": category_for(r["genotype"]),
            }
        )
    return rows


def fit(rows: list[dict]) -> dict:
    x = np.array([math.log10(r["fold"]) for r in rows])
    y = np.array([r["ddg"] for r in rows])
    lr = stats.linregress(x, y)
    rho, rho_p = stats.spearmanr(x, y)
    return {
        "n": len(rows),
        "r2": lr.rvalue ** 2,
        "p": lr.pvalue,
        "spearman_rho": rho,
        "spearman_p": rho_p,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-csv", type=Path, default=FEP_ROOT / "panel_ddg.csv")
    ap.add_argument("--out-png", type=Path,
                    default=FEP_ROOT / "panel_ddg_vs_experiment_by_category.png")
    ap.add_argument("--out-csv", type=Path,
                    default=FEP_ROOT / "panel_category_subset_stats.csv")
    args = ap.parse_args(argv)

    rows = load_rows(args.panel_csv)
    known = [r for r in rows if r["category"] != "Uncertain"]

    subsets = {
        "all": rows,
        "susceptible_plus_resistant": known,
        "resistant_only": [r for r in rows if r["category"] == "Resistant"],
        "susceptible_only": [r for r in rows if r["category"] == "Susceptible"],
    }
    stats_rows = []
    for name, sub in subsets.items():
        if len(sub) < 3:
            continue
        s = fit(sub)
        s["subset"] = name
        s["genotypes"] = " ".join(sorted(r["genotype"] for r in sub))
        stats_rows.append(s)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["subset", "n", "r2", "p", "spearman_rho", "spearman_p", "genotypes"]
        )
        w.writeheader()
        for s in stats_rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()})

    for s in stats_rows:
        print(
            f"{s['subset']:28s} n={s['n']:2d}  R2={s['r2']:.3f}  p={s['p']:.4f}   "
            f"Spearman={s['spearman_rho']:+.3f} (p={s['spearman_p']:.4f})"
        )

    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    for cat in ("Susceptible", "Resistant", "Uncertain"):
        sub = [r for r in rows if r["category"] == cat]
        if not sub:
            continue
        x = [math.log10(r["fold"]) for r in sub]
        y = [r["ddg"] for r in sub]
        e = [r["sem"] for r in sub]
        ax.errorbar(x, y, yerr=e, fmt=MARKERS[cat], color=COLORS[cat], ecolor=COLORS[cat],
                    elinewidth=1.2, capsize=2, markersize=7, linestyle="none",
                    alpha=0.45 if cat == "Uncertain" else 1.0,
                    label=f"{cat} (n={len(sub)})")
    ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)

    auto = [r for r in rows if r["genotype"] not in LABEL_OVERRIDES]
    _repel_labels(
        ax, fig,
        [math.log10(r["fold"]) for r in auto],
        [r["ddg"] for r in auto],
        [r["genotype"] for r in auto],
        fontsize=7,
    )
    for r in rows:
        off = LABEL_OVERRIDES.get(r["genotype"])
        if off is None:
            continue
        ax.annotate(
            r["genotype"],
            xy=(math.log10(r["fold"]), r["ddg"]),
            xytext=off, textcoords="offset points",
            fontsize=7, ha="right", va="center",
            arrowprops=dict(arrowstyle="-", linewidth=0.5, color="0.5",
                            shrinkA=0.0, shrinkB=3.0),
        )

    # No title: the statistics live in panel_category_subset_stats.csv and belong in
    # the figure caption, not baked into the image.
    ax.set_xlabel(r"$\log_{10}$(experimental DOR fold reduction)",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(r"Computed $\Delta\Delta G_{\mathrm{bind}}$ (kcal/mol)",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=200)
    print(f"\nWrote {args.out_png}\nWrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
