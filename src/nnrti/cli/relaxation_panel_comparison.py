#!/usr/bin/env python3
"""h-only relaxation against full minimization, across the whole panel.

Supplementary Table 5 originally compared three snapshot-relaxation protocols
on six genotypes. On that subset the two relaxed protocols looked like noisy
versions of each other. They are not: across the panel they disagree on
ddE_Total, and the disagreement sits entirely in the van der Waals term, which
is the one residual steric contacts enter. The electrostatic and polar
solvation terms, which they do not, agree well.

Reports three correlations: between the two protocols, and between each
protocol and the measured phenotype.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.relaxation_panel_comparison
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
KJ_PER_KCAL = 4.184
COMPONENTS = {"binding_dg": "ddE_Total", "binding_dg_vdw": "ddE_vdW",
              "binding_dg_electrostatic": "ddE_elec", "binding_dg_gb": "ddE_GB"}


def wt_referenced(d: pd.DataFrame, col: str) -> tuple[pd.Series, pd.Series]:
    wt = d[d.safe_label == "wt"][col].mean()
    g = d[d.safe_label != "wt"].groupby("safe_label")[col]
    return (g.mean() - wt) / KJ_PER_KCAL, g.std(ddof=1) / np.sqrt(g.size()) / KJ_PER_KCAL


def main(argv: list[str] | None = None) -> int:
    from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--h-relax", type=Path,
                    default=ROOT / "results/.checkpoints/.checkpoint_mmgbsa_relaxsens_h_relax.csv")
    ap.add_argument("--full", type=Path,
                    default=ROOT / "results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv")
    ap.add_argument("--ddg-csv", type=Path,
                    default=ROOT / "results/analysis/binding_energy/tables/ddg_full.csv")
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "results/analysis/relaxation_sensitivity/tables")
    args = ap.parse_args(argv)

    h = pd.read_csv(args.h_relax)
    f = pd.read_csv(args.full)

    # fold reduction, keyed on the safe label used by the checkpoints
    fold = (pd.read_csv(args.ddg_csv)[["mutation", "fold_reduction"]]
            .dropna().drop_duplicates("mutation"))
    fold["safe_label"] = fold.mutation.str.replace("+", "_", regex=False)
    fold = fold.set_index("safe_label").fold_reduction

    rows, summary = [], []
    for col, name in COMPONENTS.items():
        hm, hs = wt_referenced(h, col)
        fm, fs = wt_referenced(f, col)
        common = sorted(set(hm.index) & set(fm.index))
        for g in common:
            rows.append({"genotype": g.replace("_", "+"), "component": name,
                         "h_relax_mean": round(float(hm[g]), 2),
                         "h_relax_sem": round(float(hs[g]), 2),
                         "full_min_mean": round(float(fm[g]), 2),
                         "full_min_sem": round(float(fs[g]), 2)})
        a, b = hm[common].to_numpy(float), fm[common].to_numpy(float)
        sp = stats.spearmanr(a, b)
        summary.append({"component": name, "n": len(common),
                        "protocol_spearman": round(float(sp.statistic), 2),
                        "protocol_p": round(float(sp.pvalue), 3),
                        "sign_agreement": f"{int((np.sign(a) == np.sign(b)).sum())}/{len(common)}"})

    # against phenotype, on the genotypes with an established label
    established = {g.replace("+", "_") for g in RESISTANT | SUSCEPTIBLE}
    pheno = []
    for col, name in COMPONENTS.items():
        hm, _ = wt_referenced(h, col)
        fm, _ = wt_referenced(f, col)
        for label, series in (("h_relax", hm), ("full_min", fm)):
            k = sorted(set(series.index) & set(fold.index) & established)
            x = np.log10(fold[k].to_numpy(float))
            y = series[k].to_numpy(float)
            lr = stats.linregress(x, y)
            sp = stats.spearmanr(x, y)
            pheno.append({"component": name, "protocol": label, "n": len(k),
                          "R2_vs_log10_fold": round(float(lr.rvalue ** 2), 3),
                          "p": round(float(lr.pvalue), 3),
                          "spearman": round(float(sp.statistic), 2)})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "relaxation_panel_long.csv", index=False)
    s = pd.DataFrame(summary); p = pd.DataFrame(pheno)
    s.to_csv(args.output_dir / "relaxation_protocol_agreement.csv", index=False)
    p.to_csv(args.output_dir / "relaxation_phenotype_correlation.csv", index=False)

    print("agreement between the two protocols:")
    print(s.to_string(index=False))
    print("\ncorrelation with phenotype (established genotypes only):")
    print(p.to_string(index=False))
    print(f"\nwrote 3 tables to {args.output_dir.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
