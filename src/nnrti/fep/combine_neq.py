#!/usr/bin/env python3
"""Combine NEQ per-replicate BAR results into ΔΔG_bind per genotype.

ΔΔG_bind(G) = ΔG_mut^holo − ΔG_mut^apo, summed over the additive legs of a
compound genotype (see ``mutations.py::MANUSCRIPT_PLANS``). Uncertainty is the
SEM across replicates of the per-replicate target ΔΔG — not the pooled within-
replicate BAR error, which underestimates it.

Outputs (under ``results/analysis/fep_pmx/``):
    targets/{genotype}/summary.json   per-target ΔΔG_bind ± SEM + per-rep table
    panel_ddg.csv                     one row per genotype, with experimental fold
    panel_ddg_vs_experiment.png       ΔΔG_bind vs experimental fold (no fitted line)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.fep.mutations import MANUSCRIPT_PLANS
from nnrti.fep.analyze_neq import ensure_leg_analysis
from nnrti.fep.charge_correction import charge_leg_correction
from nnrti.fep.config import CHARGE_LEG_DELTA_Q, FEP_PMX_ROOT, NEQ_TEMPERATURE_K


# Genotypes in the manuscript panel that carry an experimental DOR fold reduction.
# 19 legs total; F227C has no clinical fold, so 18 points appear on the scatter.
MANUSCRIPT_PANEL_N_WITH_FOLD = 18


def _leg_charge_correction(leg_id: str, replicate: int) -> float:
    """Analytical net-charge (finite-size) correction to add to a charge leg's ΔΔG.

    Zero for neutral legs. For charge legs (CHARGE_LEG_DELTA_Q) it is the leading
    Rocklin/Hunenberger periodicity self-energy term, ΔG_holo_self - ΔG_apo_self,
    read from the leg's built holo/apo boxes. For our ~12 nm boxes it is ~1e-4
    kcal/mol (the per-phase terms cancel); see charge_correction.py. Falls back to
    0.0 if the built boxes are not present (e.g. on a machine without gromacs_build).
    """
    delta_q = CHARGE_LEG_DELTA_Q.get(leg_id)
    if delta_q is None:
        return 0.0
    leg_dir = FEP_PMX_ROOT / "legs" / leg_id
    holo_gro = leg_dir / "holo" / f"rep_{replicate:02d}" / "gromacs_build" / "system.gro"
    apo_gro = leg_dir / "apo" / f"rep_{replicate:02d}" / "gromacs_build" / "system.gro"
    if not (holo_gro.is_file() and apo_gro.is_file()):
        return 0.0
    return charge_leg_correction(holo_gro, apo_gro, delta_q=delta_q)["ddg_correction_kcal"]

DEFAULT_TARGETS = ("V106A", "Y188L")
EXPERIMENTAL_CSV = Path(
    "results/analysis/susceptibility/tables/dor_susceptibility_values.csv"
)
P0_SIGN_GATE = ("V106A", "Y188L")
# Main-text discussion gate (plan §2B). Confidence = SEM + estimator agreement,
# not experimental match. V106M is included (tight SEM) but framed as a
# binding-vs-phenotype finding, not a pipeline failure.
SEM_MAIN_TEXT_MAX_KCAL = 0.6
# The analytical Rocklin/Hunenberger correction makes the charge legs usable:
# holo and apo box volumes differ by under 1%, so the correction cancels to below
# 1e-3 kcal/mol in the double difference.
#
# The SEM threshold alone now gives the right answer for both -- K103N (0.23)
# main_text, G190E (1.63) show -- so the override is redundant as well as stale.
# Kept as an empty set rather than deleted so a genotype can be forced out of
# the main text again if a specific reason arises; state the reason here.
# (The former OMIT_MAIN_TEXT / discussion-tier machinery has been removed. It
# classified genotypes as main_text / show / omit_main from their SEM, which was
# a triage aid during the campaign and never belonged in the published output:
# every table already prints the SEM beside the value, so a reader can judge
# precision directly, and the tier strings went stale whenever a leg was
# re-run -- G190E was still marked omit_main at a SEM of 0.38.)


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

    # Independent-leg propagation: each leg is a separate simulation, so we
    # aggregate per leg (mean ± SEM over that leg's own replicates) and combine
    # legs as sum of means with SEM added in quadrature. This (a) is the correct
    # estimator for physically independent legs — rep-index pairing injects a
    # spurious per-rep correlation that inflates the SEM at small n — and (b)
    # lets a noisy leg carry more replicates than an already-tight partner leg.
    legs_agg: list[dict] = []
    for leg in legs:
        rep_vals: list[dict] = []
        for replicate in replicates:
            holo = ensure_leg_analysis(
                leg.leg_id, phase="holo", replicate=replicate,
                temperature_k=temperature_k, nboots=nboots, auto=auto, force=force,
            )
            apo = ensure_leg_analysis(
                leg.leg_id, phase="apo", replicate=replicate,
                temperature_k=temperature_k, nboots=nboots, auto=auto, force=force,
            )
            if holo.get("bar_dg") is None or apo.get("bar_dg") is None:
                continue
            charge_corr = _leg_charge_correction(leg.leg_id, replicate)
            leg_ddg = holo["bar_dg"] - apo["bar_dg"] + charge_corr
            rep_vals.append(
                {
                    "replicate": replicate,
                    "holo_dg": holo["bar_dg"],
                    "apo_dg": apo["bar_dg"],
                    "charge_correction": charge_corr,
                    "leg_ddg": leg_ddg,
                }
            )
        if not rep_vals:  # a leg with no complete replicate => genotype incomplete
            return {"genotype": genotype, "n_reps": 0, "ddg_bind": None, "sem": None,
                    "legs": legs_agg, "incomplete_leg": leg.leg_id}
        v = np.array([r["leg_ddg"] for r in rep_vals])
        nrep = len(v)
        legs_agg.append(
            {
                "leg_id": leg.leg_id,
                "n_reps": nrep,
                "leg_ddg_mean": float(v.mean()),
                "leg_ddg_sd": float(v.std(ddof=1)) if nrep > 1 else None,
                "leg_ddg_sem": float(v.std(ddof=1) / math.sqrt(nrep)) if nrep > 1 else None,
                "per_rep": rep_vals,
            }
        )

    mean = float(sum(la["leg_ddg_mean"] for la in legs_agg))
    leg_sems = [la["leg_ddg_sem"] for la in legs_agg]
    sem = float(math.sqrt(sum(s * s for s in leg_sems))) if all(s is not None for s in leg_sems) else None
    return {
        "genotype": genotype,
        "n_reps": min(la["n_reps"] for la in legs_agg),
        "n_reps_per_leg": {la["leg_id"]: la["n_reps"] for la in legs_agg},
        "ddg_bind": mean,
        "sem": sem,
        "units": "kcal/mol",
        "legs": legs_agg,
    }


def _repel_labels(ax, fig, xs, ys, labels, *, fontsize=8, n_iter=400) -> None:
    """Place scatter-point labels with a dependency-free force repel + leader lines.

    Clustered genotypes (e.g. the high-fold V106A compounds) collide with a fixed
    offset annotation; here labels are pushed apart from each other and off the
    markers in pixel space, then a thin leader connects each label to its point.
    """
    import numpy as np

    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    pts = ax.transData.transform(np.column_stack([xs, ys]))  # marker pixel coords
    lab = pts + np.array([7.0, 9.0])                          # initial label anchors (px)

    # measure each label's pixel box (left-center anchored)
    texts, sizes = [], []
    for s, (lx, ly) in zip(labels, lab):
        t = ax.text(*inv.transform((lx, ly)), s, fontsize=fontsize, ha="left",
                    va="center", zorder=6)
        bb = t.get_window_extent(renderer)
        texts.append(t)
        sizes.append((bb.width, bb.height))
    sizes = np.array(sizes)
    half = sizes / 2.0

    for _ in range(n_iter):
        moved = False
        cen = lab + half * np.array([1.0, 0.0])  # label box centers (anchor is left-center)
        for i in range(len(texts)):
            fx = fy = 0.0
            for j in range(len(texts)):
                if i == j:
                    continue
                dx = cen[i, 0] - cen[j, 0]
                dy = cen[i, 1] - cen[j, 1]
                ox = half[i, 0] + half[j, 0] - abs(dx) + 2.0
                oy = half[i, 1] + half[j, 1] - abs(dy) + 2.0
                if ox > 0 and oy > 0:  # boxes overlap -> push along cheaper axis
                    if oy <= ox:
                        fy += (np.sign(dy) or 1.0) * oy
                    else:
                        fx += (np.sign(dx) or 1.0) * ox
            for p in pts:  # keep labels off the markers
                dx = cen[i, 0] - p[0]
                dy = cen[i, 1] - p[1]
                d2 = dx * dx + dy * dy
                if 1e-6 < d2 < 20.0 ** 2:
                    d = np.sqrt(d2)
                    fx += dx / d * (20.0 - d) * 0.6
                    fy += dy / d * (20.0 - d) * 0.6
            if abs(fx) > 0.4 or abs(fy) > 0.4:
                lab[i] += np.array([fx, fy]) * 0.5
                moved = True
        if not moved:
            break

    for t, (lx, ly), px, py in zip(texts, lab, xs, ys):
        dx, dy = inv.transform((lx, ly))
        t.set_position((dx, dy))
        px_disp, py_disp = ax.transData.transform((px, py))
        if np.hypot(lx - px_disp, ly - py_disp) > 14.0:  # draw a leader when offset
            ax.annotate("", xy=(px, py), xytext=(dx, dy), zorder=1,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="0.6"))


def _plot(rows: list[dict], rho: float | None, output: Path) -> dict:
    """Scatter ΔΔG vs log10(fold). No fitted line — weak correlation is the finding."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    computed = [r for r in rows if r["ddg_bind"] is not None]
    stats: dict = {"n_with_fold": 0, "pearson_r2": None, "pearson_p": None, "spearman_rho": rho}
    if not computed:
        return stats
    with_fold = [r for r in computed if r.get("fold") is not None]

    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    if with_fold:
        x = np.array([math.log10(r["fold"]) for r in with_fold])
        y = np.array([r["ddg_bind"] for r in with_fold])
        yerr = np.array([r["sem"] if r["sem"] is not None else 0.0 for r in with_fold])
        ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, color="#2c6fbb", ecolor="#9bbce0", zorder=3)
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
        _repel_labels(ax, fig, x, y, [r["genotype"] for r in with_fold], fontsize=8)
        corr_label = None
        if len(with_fold) >= 3 and x.std() > 0:
            try:
                from scipy.stats import pearsonr
                r_, p_ = pearsonr(x, y)
                stats["pearson_r2"] = float(r_ ** 2)
                stats["pearson_p"] = float(p_)
                corr_label = f"Pearson R² = {r_**2:.2f}, p = {p_:.2g}  (no fit line)"
            except Exception:
                r_ = float(np.corrcoef(x, y)[0, 1])
                stats["pearson_r2"] = float(r_ ** 2)
                corr_label = f"Pearson R² = {r_**2:.2f}  (no fit line)"
        if rho is not None:
            spearman_bit = f"Spearman ρ = {rho:.2f}"
            corr_label = f"{corr_label}  ·  {spearman_bit}" if corr_label else spearman_bit
        ax.set_xlabel(r"$\log_{10}$(experimental DOR fold reduction)")
        ax.set_ylabel(r"Computed $\Delta\Delta G_{\mathrm{bind}}$ (kcal/mol)")
        # MANUSCRIPT_PANEL_N_WITH_FOLD is how many of the 19 legs carry an
        # experimental fold (F227C has none).
        n_done = len(with_fold)
        completeness = "" if n_done >= MANUSCRIPT_PANEL_N_WITH_FOLD else "; incomplete panel"
        title = f"NEQ ΔΔG_bind vs experiment (n = {n_done}{completeness})"
        if corr_label:
            title += f"\n{corr_label}"
        stats["n_with_fold"] = len(with_fold)
    else:
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
    return stats


def rows_from_panel_csv(panel_csv: Path) -> list[dict]:
    """Rebuild plot rows from an existing panel_ddg.csv (no re-analysis)."""
    rows = []
    with panel_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            ddg = row.get("ddg_bind_kcal", "").strip()
            sem = row.get("sem_kcal", "").strip()
            fold = row.get("dor_fold_reduction", "").strip()
            n = row.get("n_reps", "").strip()
            rows.append(
                {
                    "genotype": row["genotype"],
                    "ddg_bind": float(ddg) if ddg else None,
                    "sem": float(sem) if sem else None,
                    "n_reps": int(n) if n else 0,
                    "fold": float(fold) if fold else None,
                }
            )
    return rows


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
    parser.add_argument(
        "--replot-only",
        action="store_true",
        help="Rebuild panel_ddg_vs_experiment.png from existing panel_ddg.csv",
    )
    args = parser.parse_args(argv)

    if args.replot_only:
        panel_csv = args.output_dir / "panel_ddg.csv"
        rows = rows_from_panel_csv(panel_csv)
        with_fold = [r for r in rows if r.get("ddg_bind") is not None and r.get("fold") is not None]
        rho = None
        if len(with_fold) >= 3:
            rho = spearman_rho(
                np.array([r["ddg_bind"] for r in with_fold]),
                np.array([r["fold"] for r in with_fold]),
            )
        plot_path = args.output_dir / "panel_ddg_vs_experiment.png"
        stats = _plot(rows, rho, plot_path)
        print(f"replot n={stats.get('n_with_fold')}  R²={stats.get('pearson_r2')}  "
              f"p={stats.get('pearson_p')}  Spearman ρ={rho}")
        print(f"Wrote {plot_path}")
        return 0

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
    x = y = np.array([])
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
                r.get("fold", "") if r.get("fold") is not None else "",
            ])

    plot_path = args.output_dir / "panel_ddg_vs_experiment.png"
    stats = _plot(rows, rho, plot_path)

    # Report
    print(f"{'genotype':<14} {'ddg_bind':>10} {'sem':>7} {'n':>3}  {'fold':>7}")
    for r in rows:
        ddg = f"{r['ddg_bind']:+.2f}" if r["ddg_bind"] is not None else "  n/a"
        sem = f"{r['sem']:.2f}" if r.get("sem") is not None else "  - "
        fold = f"{r['fold']:.1f}" if r.get("fold") is not None else "   -"
        print(f"{r['genotype']:<14} {ddg:>10} {sem:>7} {r['n_reps']:>3}  {fold:>7}")
    if rho is not None:
        print(f"\nSpearman ρ (ΔΔG_bind vs fold) = {rho:.3f}  (n = {len(x)})")
    if stats.get("pearson_r2") is not None:
        print(f"Pearson R² = {stats['pearson_r2']:.3f}  p = {stats.get('pearson_p')}")
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
