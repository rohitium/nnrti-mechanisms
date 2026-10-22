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
    # <w:t[^>]*> also matches <w:tc>, <w:tbl>, <w:tcPr> ... so a document
    # containing tables came back with markup in its "text". The space is what
    # separates the text element from its attributes.
    return re.sub(r"\s+", " ",
                  "".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>", x, re.S)))


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

    # Charge-model comparison. The headline claim is a bound on how far the
    # WT-referenced ddE moves between charge sets, so the bound is checked
    # against the table rather than against a remembered value.
    cm = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                     "charge_model_comparison.csv")
    tot = cm[cm.component == "ddE_Total"]
    checks += [
        ("charge model: all legs have 3 replicates",
         bool((cm.n_gast == 3).all() and (cm.n_am1bcc == 3).all())),
        ("charge model: ddE_Total shift at most 0.81",
         abs(tot.difference.abs().max() - 0.807) < 0.005),
        ("charge model: sign preserved both genotypes",
         bool((tot.mean_gast * tot.mean_am1bcc > 0).all())),
        ("charge model: Y188L 2.78 -> 2.34",
         abs(float(tot[tot.genotype == "Y188L"].mean_am1bcc.iloc[0]) - 2.337) < 0.005),
        ("charge model: V106A 1.73 -> 2.54",
         abs(float(tot[tot.genotype == "V106A"].mean_am1bcc.iloc[0]) - 2.537) < 0.005),
    ]

    supp = doc_text(ROOT / "revisions/Supplementary-Text-Revision-1.docx")
    resp = doc_text(ROOT / "revisions/Response-to-Reviewers.docx")
    for s in ("2.78 \u2192 2.34", "up to 1.17", "0.23\u20130.29"):
        checks.append((f"supplementary contains {s!r}", s in supp))
    for s in ("2.78 to 2.34", "up to 1.17", "0.23 to 0.29"):
        checks.append((f"response contains {s!r}", s in resp))
    checks.append(("response no longer says a replicate is pending",
                   "still computing" not in resp))

    for s in ("1.85", "−1.35"):
        checks.append((f"manuscript contains {s!r}", s in ms))
    # The finite-size derivation lives in Supplementary Note 1; Methods keeps
    # only the bound and a pointer to it, so these are checked there. That
    # way shortening Methods does not fail the check.
    for s in ("0.003 kcal/mol", "0.019%", "0.13%", "0.31%"):
        checks.append((f"supplementary contains {s!r}", s in supp))
    checks.append(("manuscript keeps the finite-size bound or a pointer",
                   "0.003 kcal/mol" in ms or "Supplementary Note 1" in ms))
    # The relaxation-protocol numbers moved out of Methods and into the SI in
    # revision, so they are checked where they now live, not where they were.
    rel = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                      "relaxation_protocol_agreement.csv").set_index("component")
    ph = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                     "relaxation_phenotype_correlation.csv")
    tot = ph[ph.component == "ddE_Total"].set_index("protocol").spearman
    checks += [
        ("relaxation: protocols agree on elec and GB",
         rel.loc["ddE_elec", "protocol_spearman"] >= 0.6
         and rel.loc["ddE_GB", "protocol_spearman"] >= 0.6),
        ("relaxation: protocols disagree on Total",
         abs(rel.loc["ddE_Total", "protocol_spearman"]) < 0.3),
        ("relaxation: phenotype correlations have opposite sign",
         float(tot["full_min"]) > 0 > float(tot["h_relax"])),
        ("supplementary quotes the protocol correlations",
         "0.69" in supp and "0.81" in supp),
        ("supplementary quotes the phenotype correlations",
         "+0.45" in supp and "-0.32" in supp),
    ]
    checks.append(("methods no longer reports the sensitivity result",
                   "materially sensitive" not in ms))
    checks.append(("supplementary has the apo note and table",
                   "Supplementary Note 3" in supp and "Supplementary Table 7" in supp))

    # Literal citation numbers in added (red) text. A digit typed after a word
    # is not a Zotero field: it points at whatever reference held that number
    # when it was typed, and it does not move when the bibliography renumbers.
    # Three got in this way; "pmx,46" and "form48" pointed at MM-PBSA papers
    # and "Rocklin scheme,53" pointed at Crooks.
    import html as _html
    for doc_name in ("Manuscript-Revision-1.docx", "Supplementary-Text-Revision-1.docx"):
        raw = zipfile.ZipFile(ROOT / "revisions" / doc_name).read(
            "word/document.xml").decode("utf-8")
        stray = []
        for m in re.finditer(
                r"<w:r>(<w:rPr>(?:(?!</w:rPr>).)*<w:color w:val=\"FF0000\"/>.*?</w:rPr>)"
                r"<w:t(?: [^>]*)?>(.*?)</w:t></w:r>", raw, re.S):
            t = _html.unescape(m.group(2))
            for c in re.finditer(r"(\S{0,12}?)(\d{1,3})(?=[ .,;:]|$)", t):
                word, num = c.group(1), c.group(2)
                # residue names (Ser105, Tyr188) and mutation codes (V106,
                # G190E, K103N) end in digits without being citations
                if re.search(r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
                             r"Phe|Pro|Ser|Thr|Trp|Tyr|Val)$", word):
                    continue
                if re.search(r"log$", word, re.I):      # log10, log2
                    continue
                if re.fullmatch(r"[A-Z]?", word) or re.search(r"[A-Z]$", word):
                    continue
                if not re.search(r"[a-z],?$", word):
                    continue
                ctx = t[max(0, c.start() - 45):c.end() + 10]
                if re.search(r"\d\s*(ns|\u00c5|kcal|nm|%|\u00d7|fs|bar|K|genotypes|"
                             r"replicates|frames|histidines|of|\u00b0)", ctx):
                    continue
                stray.append(ctx)
        checks.append((f"{doc_name.split('-')[0]} has no literal citation numbers "
                       f"in added text", not stray))

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
