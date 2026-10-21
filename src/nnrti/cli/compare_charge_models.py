#!/usr/bin/env python3
"""Compare MM/GBSA energies under Gasteiger and AM1-BCC ligand charges.

Reviewer 2 (comment 5) asked us either to compare against a more accurate
charge model for a subset of systems, or to discuss how the choice affects the
computed energies. This does the comparison.

The same trajectories are re-scored under both charge sets with everything else
held fixed -- same snapshots, same sampling, same relaxation protocol -- so the
only difference between the two columns is the ligand charges. The Gasteiger
column is the published computation, because that is the number in the main
text rather than a re-run of it.

What to read: the absolute energies are expected to move, since AM1-BCC is
markedly more polarised (mean |q| 0.228 e against 0.138 e) and the
electrostatic and polar-solvation terms scale with the ligand charges. The
question that matters for the paper is whether the WT-referenced ddE, which is
what Table 2 reports, moves with them.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_charge_models
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

KJ_PER_KCAL = 4.184
COMPONENTS = {"binding_dg": "ddE_Total", "binding_dg_vdw": "ddE_vdW",
              "binding_dg_electrostatic": "ddE_elec", "binding_dg_gb": "ddE_GB"}


def wt_referenced(d: pd.DataFrame) -> pd.DataFrame:
    wt = d[d.safe_label == "wt"]
    rows = []
    for col, name in COMPONENTS.items():
        if col not in d:
            continue
        ref = wt[col].mean()
        for lbl, g in d[d.safe_label != "wt"].groupby("safe_label"):
            dd = (g[col] - ref) / KJ_PER_KCAL
            rows.append({"genotype": lbl, "component": name,
                         "mean": dd.mean(), "sem": dd.sem(), "n": len(dd)})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--am1bcc", type=Path,
                   default=Path("results/.checkpoints/.checkpoint_mmgbsa_am1bcc.csv"))
    p.add_argument("--gasteiger", type=Path,
                   default=Path("results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/relaxation_sensitivity/tables"))
    args = p.parse_args(argv)

    a = pd.read_csv(args.am1bcc)
    subset = sorted(a.safe_label.unique())
    g = pd.read_csv(args.gasteiger)
    g = g[g.safe_label.isin(subset)]

    print(f"subset: {subset}")
    print(f"replicates: AM1-BCC {len(a)}, Gasteiger {len(g)}\n")

    print("absolute binding energy, kcal/mol (mean over replicates):")
    abs_cmp = pd.DataFrame({
        "gasteiger": g.groupby("safe_label").binding_dg.mean() / KJ_PER_KCAL,
        "am1bcc": a.groupby("safe_label").binding_dg.mean() / KJ_PER_KCAL,
    })
    abs_cmp["difference"] = abs_cmp.am1bcc - abs_cmp.gasteiger
    print(abs_cmp.round(1).to_string())

    ga, aa = wt_referenced(g), wt_referenced(a)
    m = ga.merge(aa, on=["genotype", "component"], suffixes=("_gast", "_am1bcc"))
    m["difference"] = m["mean_am1bcc"] - m["mean_gast"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    m.round(3).to_csv(args.output_dir / "charge_model_comparison.csv", index=False)

    print("\nWT-referenced ddE, kcal/mol:")
    show = m[["genotype", "component", "mean_gast", "sem_gast",
              "mean_am1bcc", "sem_am1bcc", "difference"]].copy()
    for c in show.columns[2:]:
        show[c] = show[c].round(2)
    print(show.to_string(index=False))

    tot = m[m.component == "ddE_Total"]
    if len(tot):
        worst = tot.difference.abs().max()
        print(f"\nddE_Total shifts by at most {worst:.2f} kcal/mol between charge models.")
        same_sign = (np.sign(tot.mean_gast) == np.sign(tot.mean_am1bcc)).all()
        print(f"Sign preserved for every genotype: {same_sign}")
    print(f"\nwrote {args.output_dir}/charge_model_comparison.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
