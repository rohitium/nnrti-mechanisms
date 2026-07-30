#!/usr/bin/env python3
"""Combine NEQ per-replicate BAR results into ΔΔG_bind per genotype.

ΔΔG_bind(G) = ΔG_mut^holo − ΔG_mut^apo, summed over the additive legs of a
compound genotype (see ``mutations.py::MANUSCRIPT_PLANS``). Uncertainty is the
SEM across replicates of the per-replicate target ΔΔG — not the pooled within-
replicate BAR error, which underestimates it (``pmx-neq-fep-plan.md`` §3.2).

Outputs (under ``results/analysis/fep_pmx/``):
  targets/{genotype}/summary.json   per-target ΔΔG_bind ± SEM + per-rep table
  panel_ddg.csv                     one row per genotype, with experimental fold
  panel_ddg_vs_experiment.png       ΔΔG_bind vs experimental fold (Spearman ρ)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_jorgensen.mutations import MANUSCRIPT_PLANS
from scripts.fep_pmx.analyze_neq import ensure_leg_analysis
from scripts.fep_pmx.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K

DEFAULT_TARGETS = ("V106A", "Y188L")
EXPERIMENTAL_CSV = Path(
    "results/analysis/dor_susceptibility_bar_chart/tables/dor_susceptibility_values.csv"
)
P0_SIGN_GATE = ("V106A", "Y188L")


def load_experimental(csv_path: Path) -> dict[str, float]:
    fold: dict[str, float] = {}
    if not csv_path.is_file():
        return fold
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                fold[row["mutation"]] = float(row["dor_fold_reduction"])
            except (KeyError, ValueError):
                continue
    return fold


def _rank(values: np.ndarray) -> np.ndarray:
    """Average ranks (ties shared), matching scipy.stats.rankdata('average')."""
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(1, len(values) + 1)
    # average tied ranks
    _, inv, counts = np.unique(values, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 3:
        return None
    rx, ry = _rank(x), _rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def target_ddg(
    genotype: str,
    *,
    replicates: range,
    temperature_k: float,
    nboots: int,
    auto: bool,
    force: bool = False,
) -> dict:
    """Per-replicate and combined ΔΔG_bind for one genotype (sum of its legs)."""
    if genotype not in MANUSCRIPT_PLANS:
        raise ValueError(f"Unknown genotype {genotype}; not in MANUSCRIPT_PLANS")
    legs = MANUSCRIPT_PLANS[genotype].legs

    per_rep: list[dict] = []
    for replicate in replicates:
        legs_detail = []
        rep_ddg = 0.0
        complete = True
        for leg in legs:
            holo = ensure_leg_analysis(
                leg.leg_id, phase="holo", replicate=replicate,
                temperature_k=temperature_k, nboots=nboots, auto=auto, force=force,
            )
            apo = ensure_leg_analysis(
                leg.leg_id, phase="apo", replicate=replicate,
                temperature_k=temperature_k, nboots=nboots, auto=auto, force=force,
            )
            if holo.get("bar_dg") is None or apo.get("bar_dg") is None:
                complete = False
                break
            leg_ddg = holo["bar_dg"] - apo["bar_dg"]
            rep_ddg += leg_ddg
            legs_detail.append(
                {
                    "leg_id": leg.leg_id,
                    "holo_dg": holo["bar_dg"],
                    "apo_dg": apo["bar_dg"],
                    "holo_err": holo.get("bar_err_analytical"),
                    "apo_err": apo.get("bar_err_analytical"),
                    "leg_ddg": leg_ddg,
                }
            )
        if complete:
            per_rep.append({"replicate": replicate, "ddg_bind": rep_ddg, "legs": legs_detail})

    if not per_rep:
        return {"genotype": genotype, "n_reps": 0, "ddg_bind": None, "sem": None, "per_rep": []}

    ddgs = np.array([r["ddg_bind"] for r in per_rep])
    n = len(ddgs)
    mean = float(ddgs.mean())
    sem = float(ddgs.std(ddof=1) / math.sqrt(n)) if n > 1 else None
    return {
        "genotype": genotype,
        "n_reps": n,
        "ddg_bind": mean,
        "sem": sem,
        "units": "kcal/mol",
        "per_rep": per_rep,
    }


def _plot(rows: list[dict], rho: float | None, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    computed = [r for r in rows if r["ddg_bind"] is not None]
    if not computed:
        return
    with_fold = [r for r in computed if r.get("fold") is not None]

    fig, ax = plt.subplots(figsize=(6, 5))
    if with_fold:
        # ΔΔG_bind vs experimental fold reduction (the ranking gate)
        x = [r["ddg_bind"] for r in with_fold]
        y = [r["fold"] for r in with_fold]
        xerr = [r["sem"] if r["sem"] is not None else 0.0 for r in with_fold]
        ax.errorbar(x, y, xerr=xerr, fmt="o", capsize=3, color="#2c6fbb", ecolor="#9bbce0")
        for r in with_fold:
            ax.annotate(r["genotype"], (r["ddg_bind"], r["fold"]),
                        textcoords="offset points", xytext=(5, 4), fontsize=8)
        ax.axvline(0.0, color="0.6", lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel(r"Computed $\Delta\Delta G_{\mathrm{bind}}$ (kcal/mol)")
        ax.set_ylabel("Experimental DOR fold reduction")
        title = "NEQ ΔΔG_bind vs experiment"
        if rho is not None:
            title += f"  (Spearman ρ = {rho:.2f}, n = {len(with_fold)})"
    else:
        # No experimental values yet (e.g. P0): show ΔΔG_bind ± SEM per genotype
        labels = [r["genotype"] for r in computed]
        vals = [r["ddg_bind"] for r in computed]
        errs = [r["sem"] if r.get("sem") is not None else 0.0 for r in computed]
        ax.bar(labels, vals, yerr=errs, capsize=4, color="#2c6fbb")
        ax.axhline(0.0, color="0.4", lw=0.8)
        ax.set_ylabel(r"Computed $\Delta\Delta G_{\mathrm{bind}}$ (kcal/mol)")
        ax.set_xlabel("genotype")
        title = "NEQ ΔΔG_bind (no experimental values loaded)"
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Combine NEQ BAR results into ΔΔG_bind per genotype.")
    parser.add_argument("--targets", nargs="+", default=list(DEFAULT_TARGETS))
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    parser.add_argument("--nboots", type=int, default=100)
    parser.add_argument("--experimental-csv", type=Path, default=EXPERIMENTAL_CSV)
    parser.add_argument("--no-auto", action="store_true", help="Do not run analyze_neq for missing legs")
    parser.add_argument("--force", action="store_true", help="Re-run pmx analyse even if analysis.json exists (use after re-running switches)")
    parser.add_argument("--strict", action="store_true", help="Exit non-zero if a P0 sign gate fails")
    parser.add_argument("--output-dir", type=Path, default=FEP_PMX_ROOT)
    args = parser.parse_args(argv)

    replicates = range(1, args.replicates + 1)
    experimental = load_experimental(args.experimental_csv)

    rows: list[dict] = []
    for genotype in args.targets:
        try:
            res = target_ddg(
                genotype, replicates=replicates,
                temperature_k=args.temperature_k, nboots=args.nboots, auto=not args.no_auto,
                force=args.force,
            )
        except (FileNotFoundError, ValueError) as exc:
            print(f"skip {genotype}: {exc}", file=sys.stderr)
            continue
        res["fold"] = experimental.get(genotype)
        rows.append(res)
        if res["n_reps"]:
            tgt_dir = args.output_dir / "targets" / genotype
            tgt_dir.mkdir(parents=True, exist_ok=True)
            (tgt_dir / "summary.json").write_text(json.dumps(res, indent=2) + "\n")

    computed = [r for r in rows if r["ddg_bind"] is not None]
    rho = None
    if computed:
        x = np.array([r["ddg_bind"] for r in computed if r.get("fold") is not None])
        y = np.array([r["fold"] for r in computed if r.get("fold") is not None])
        if len(x) >= 3:
            rho = spearman_rho(x, y)

    # Panel CSV
    panel_csv = args.output_dir / "panel_ddg.csv"
    panel_csv.parent.mkdir(parents=True, exist_ok=True)
    with panel_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genotype", "ddg_bind_kcal", "sem_kcal", "n_reps", "dor_fold_reduction"])
        for r in rows:
            writer.writerow([
                r["genotype"],
                f"{r['ddg_bind']:.3f}" if r["ddg_bind"] is not None else "",
                f"{r['sem']:.3f}" if r.get("sem") is not None else "",
                r["n_reps"],
                r.get("fold", ""),
            ])

    plot_path = args.output_dir / "panel_ddg_vs_experiment.png"
    _plot(rows, rho, plot_path)

    # Report
    print(f"{'genotype':<14} {'ddg_bind':>10} {'sem':>7} {'n':>3}  {'fold':>7}")
    for r in rows:
        ddg = f"{r['ddg_bind']:+.2f}" if r["ddg_bind"] is not None else "  n/a"
        sem = f"{r['sem']:.2f}" if r.get("sem") is not None else "  - "
        fold = f"{r['fold']:.1f}" if r.get("fold") is not None else "   -"
        print(f"{r['genotype']:<14} {ddg:>10} {sem:>7} {r['n_reps']:>3}  {fold:>7}")
    if rho is not None:
        print(f"\nSpearman ρ (ΔΔG_bind vs fold) = {rho:.3f}  (n = {len(x)})")
    print(f"\nWrote {panel_csv}\nWrote {plot_path}")

    # P0 sign gate
    gate_fail = False
    for r in rows:
        if r["genotype"] in P0_SIGN_GATE and r["ddg_bind"] is not None:
            ok = r["ddg_bind"] > 0
            gate_fail = gate_fail or not ok
            print(f"sign gate {r['genotype']}: ΔΔG_bind = {r['ddg_bind']:+.2f} "
                  f"→ {'PASS' if ok else 'FAIL (expected positive / resistance)'}")

    if args.strict and gate_fail:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
