#!/usr/bin/env python3
"""Classification performance of FEP and MM/GBSA against established phenotype.

Turns the paper's counts -- "3 of 4 susceptible, 7 of 9 resistant" -- into
standard classification metrics. Two framings are computed, answering different
questions:

**Threshold classification.** Applies the manuscript's own rule -- a genotype is
called "strong impact" when its estimate lies more than one standard error above
0.5 kcal/mol -- and scores it against the established Susceptible/Resistant
labels. Reports the 2x2 table, sensitivity, specificity, accuracy, Matthews
correlation coefficient and a Fisher exact test. MCC is the right headline for a
small, imbalanced panel: accuracy alone is flattered by the 9:4 class ratio,
and MCC is only high when all four cells behave.

**Threshold-free ranking.** ROC AUC treats the raw energy as a continuous score,
so it measures whether the method *ranks* resistant above susceptible genotypes
regardless of where the cutoff sits. This is the fairer question for a physical
observable that was never calibrated to a clinical threshold, and it does not
reward a threshold that happens to suit this panel. Confidence interval by
bootstrap; p-value by exact permutation, which is feasible at n = 13.

Uncertain-phenotype genotypes are excluded throughout: they have no ground truth
to score against.

Usage
-----
    python -m nnrti.cli.classification_performance
    python -m nnrti.cli.classification_performance --include-uncertain
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN

#: thermal energy at 300 K; the manuscript's "strong impact" threshold
THRESHOLD = 0.5


def confusion(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    n = tp + tn + fp + fn
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "accuracy": (tp + tn) / n if n else float("nan"),
        "sensitivity": tp / (tp + fn) if (tp + fn) else float("nan"),
        "specificity": tn / (tn + fp) if (tn + fp) else float("nan"),
        "ppv": tp / (tp + fp) if (tp + fp) else float("nan"),
        "npv": tn / (tn + fn) if (tn + fn) else float("nan"),
        "mcc": (tp * tn - fp * fn) / denom if denom else float("nan"),
        "fisher_p": float(stats.fisher_exact([[tp, fp], [fn, tn]], alternative="greater")[1]),
    }


def auc(y_true: np.ndarray, score: np.ndarray) -> float:
    """Mann-Whitney U form of the ROC AUC; ties count as half."""
    pos, neg = score[y_true == 1], score[y_true == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (pos.size * neg.size))


def auc_permutation_p(y_true: np.ndarray, score: np.ndarray, max_exact: int = 20000) -> tuple[float, str]:
    """Exact permutation p for AUC when the panel is small enough to enumerate."""
    n, k = len(y_true), int(y_true.sum())
    observed = auc(y_true, score)
    n_comb = math.comb(n, k)
    if n_comb <= max_exact:
        count = 0
        for idx in itertools.combinations(range(n), k):
            lab = np.zeros(n, dtype=int)
            lab[list(idx)] = 1
            if auc(lab, score) >= observed:
                count += 1
        return count / n_comb, f"exact ({n_comb} labelings)"
    rng = np.random.default_rng(0)
    draws = [auc(rng.permutation(y_true), score) for _ in range(max_exact)]
    return float((np.asarray(draws) >= observed).mean()), f"monte carlo ({max_exact})"


def auc_ci(y_true: np.ndarray, score: np.ndarray, n_boot: int = 10000) -> tuple[float, float]:
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        vals.append(auc(y_true[idx], score[idx]))
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def evaluate(name: str, est: pd.Series, sem: pd.Series, include_uncertain: bool) -> dict:
    # sorted so the genotype order in the output does not depend on set iteration
    labels = {g: 1 for g in sorted(RESISTANT)} | {g: 0 for g in sorted(SUSCEPTIBLE)}
    genos = [g for g in labels if g in est.index and pd.notna(est[g])]
    y = np.array([labels[g] for g in genos])
    score = est[genos].to_numpy(dtype=float)
    # manuscript rule: strong when the estimate lies more than one SEM above 0.5
    pred = ((est[genos] - sem[genos]) > THRESHOLD).to_numpy().astype(int)

    res = {"method": name, "n": len(genos),
           "n_resistant": int(y.sum()), "n_susceptible": int((y == 0).sum())}
    res.update(confusion(y, pred))
    res["auc"] = auc(y, score)
    lo, hi = auc_ci(y, score)
    res["auc_ci_low"], res["auc_ci_high"] = lo, hi
    p, how = auc_permutation_p(y, score)
    res["auc_p"], res["auc_p_method"] = p, how
    res["misclassified"] = ", ".join(
        f"{g} ({'FN' if labels[g] == 1 else 'FP'})"
        for g, pr in zip(genos, pred) if pr != labels[g])
    return res


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ddg-csv", type=Path,
                    default=root / "results/analysis/binding_energy/tables/ddg_full.csv")
    ap.add_argument("--fep-csv", type=Path,
                    default=root / "results/analysis/fep_pmx/panel_ddg.csv")
    ap.add_argument("--include-uncertain", action="store_true")
    ap.add_argument("--output-dir", type=Path,
                    default=root / "results/analysis/classification_performance")
    args = ap.parse_args()

    fep = pd.read_csv(args.fep_csv).set_index("genotype")
    mm = pd.read_csv(args.ddg_csv)
    mm = mm[mm.mutation != "WT"].groupby("mutation")
    mm_est = mm["ddg"].mean()
    mm_sem = mm["ddg"].apply(lambda x: x.std(ddof=1) / np.sqrt(x.size))

    rows = [
        evaluate("FEP ddG_bind", fep["ddg_bind_kcal"], fep["sem_kcal"], args.include_uncertain),
        evaluate("MM/GBSA ddE_Total", mm_est, mm_sem, args.include_uncertain),
    ]

    for r in rows:
        print(f"\n=== {r['method']}  (n={r['n']}: {r['n_resistant']} resistant, "
              f"{r['n_susceptible']} susceptible) ===")
        print(f"  confusion      TP {r['TP']}  FP {r['FP']}  FN {r['FN']}  TN {r['TN']}")
        print(f"  sensitivity    {r['sensitivity']:.2f}   ({r['TP']}/{r['TP']+r['FN']} resistant detected)")
        print(f"  specificity    {r['specificity']:.2f}   ({r['TN']}/{r['TN']+r['FP']} susceptible correct)")
        print(f"  accuracy       {r['accuracy']:.2f}")
        print(f"  PPV / NPV      {r['ppv']:.2f} / {r['npv']:.2f}")
        print(f"  MCC            {r['mcc']:.2f}")
        print(f"  Fisher exact   p = {r['fisher_p']:.3f}")
        print(f"  ROC AUC        {r['auc']:.2f}  (95% CI {r['auc_ci_low']:.2f}-{r['auc_ci_high']:.2f}), "
              f"p = {r['auc_p']:.3f} [{r['auc_p_method']}]")
        if r["misclassified"]:
            print(f"  misclassified  {r['misclassified']}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "classification_metrics.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    (args.output_dir / "classification_metrics.json").write_text(json.dumps(rows, indent=2))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
