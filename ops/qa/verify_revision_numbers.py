#!/usr/bin/env python3
"""Check that numbers quoted in the revision documents match the committed data.

Every figure written into the manuscript, the response or the Supporting
Information came from a CSV or JSON in results/. Those files get regenerated;
the prose does not. This re-derives each quoted number from its source and
fails loudly on a mismatch, so a stale number cannot survive a rerun unnoticed.

Run it before any submission.

Usage
-----
    PYTHONPATH=src python ops/qa/verify_revision_numbers.py
"""

from __future__ import annotations

import json
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def doc_text(path: Path) -> str:
    x = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    return re.sub(r"\s+", " ", "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", x, re.S)))


def main() -> int:
    ms = doc_text(ROOT / "revisions/Manuscript-Revision-1.docx")
    checks: list[tuple[str, bool]] = []

    j = json.loads((ROOT / "results/analysis/classification_performance/"
                    "classification_metrics.json").read_text())[0]
    checks += [("classification sensitivity 0.78", f"{j['sensitivity']:.2f}" == "0.78"),
               ("classification specificity 0.75", f"{j['specificity']:.2f}" == "0.75"),
               ("classification MCC 0.50", f"{j['mcc']:.2f}" == "0.50"),
               ("AUC 0.78", f"{j['auc']:.2f}" == "0.78"),
               ("AUC permutation p 0.074", f"{j['auc_p']:.3f}" == "0.074")]

    es = pd.read_csv(ROOT / "results/analysis/effect_sizes/tables/effect_sizes.csv")
    r = es[(es.comparison == "V106A -> V106I") & (es.metric == "dor_to_S105_mindist")].iloc[0]
    checks.append(("V106A->V106I Ser105 +1.85 [1.25, 2.45]",
                   abs(r.difference - 1.85) < 0.01 and abs(r.welch_ci_low - 1.25) < 0.01
                   and abs(r.welch_ci_high - 2.45) < 0.01))
    r = es[(es.comparison == "wt -> V106A") & (es.metric == "dor_to_S105_mindist")].iloc[0]
    checks.append(("wt->V106A Ser105 -1.35 [-1.73, -0.97]",
                   abs(r.difference + 1.35) < 0.01 and abs(r.welch_ci_low + 1.73) < 0.01))
    checks.append(("7 of 42 intervals exclude zero",
                   int((~es.crosses_zero).sum()) == 7 and len(es) == 42))

    fs = json.loads((ROOT / "results/analysis/fep_pmx/charge_correction_residual.json").read_text())
    checks.append(("finite-size residual <= 0.003 kcal/mol",
                   round(max(v["residual_bound_kcal"] for v in fs["legs"].values()), 3) == 0.003))
    checks.append(("box mismatch 0.13% (K103N) and 0.31% (G190E)",
                   abs(fs["legs"]["wt_to_K103N"]["box_mismatch_fraction"] - 0.001295) < 1e-5
                   and abs(fs["legs"]["wt_to_G190E"]["box_mismatch_fraction"] - 0.003146) < 1e-5))

    t5 = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                     "relaxation_sensitivity_ddE_total.csv", index_col=0)
    checks.append(("Table 5 unrestrained SEMs span 0.13-0.54",
                   abs(t5.unrestrained_sem.min() - 0.13) < 0.01
                   and abs(t5.unrestrained_sem.max() - 0.54) < 0.01))
    checks.append(("Table 5 G190E h_relax 28.3 +/- 15.1",
                   abs(t5.loc["G190E", "h_relax_mean"] - 28.32) < 0.05
                   and abs(t5.loc["G190E", "h_relax_sem"] - 15.12) < 0.05))

    apo = ROOT / "results/analysis/apo_ensembles/tables/apo_matched_window_summary.csv"
    if apo.exists():
        a = pd.read_csv(apo, index_col=0)
        if "dVolume_vs_WT" in a:
            for g, v in (("g190e", -49.0), ("a98g_f227c", -28.4), ("v106i_f227c", -36.6)):
                if g in a.index:
                    checks.append((f"apo dVolume {g} {v}",
                                   abs(float(a.loc[g, "dVolume_vs_WT"]) - v) < 0.6))

    for s in ("0.13 to 0.54", "0.003 kcal/mol", "0.019%", "1.85", "−1.35"):
        checks.append((f"manuscript contains {s!r}", s in ms))

    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"{name:{width}s}  {'PASS' if ok else 'FAIL'}")
    failed = [n for n, ok in checks if not ok]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} consistent")
    if failed:
        print("\nFAILED:")
        for n in failed:
            print(f"  - {n}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
