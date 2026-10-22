#!/usr/bin/env python3
"""Build the Supporting Information table workbooks.

Supplementary Tables 1-4 are XLSX, so Tables 5-7 are too rather than being
embedded in the Supplementary Text. Each workbook carries the summary the
text refers to, the per-genotype values behind it, and a Notes sheet stating
how the numbers were produced.

Usage
-----
    PYTHONPATH=src python ops/docx/build_si_workbooks.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revisions"
REL = ROOT / "results/analysis/relaxation_sensitivity/tables"


def _write(path: Path, sheets: dict[str, pd.DataFrame], notes: list[str]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name[:31], index=False)
        pd.DataFrame({"Notes": notes}).to_excel(
            xl, sheet_name="Notes", index=False, header=False)
    print(f"wrote {path.relative_to(ROOT)}  sheets: "
          f"{list(sheets) + ['Notes']}")


def table5() -> None:
    panel = pd.read_csv(REL / "relaxation_panel_long.csv")
    agree = pd.read_csv(REL / "relaxation_protocol_agreement.csv")
    pheno = pd.read_csv(REL / "relaxation_phenotype_correlation.csv")
    per_rep = pd.read_csv(
        ROOT / "results/.checkpoints/.checkpoint_mmgbsa_relaxsens_h_relax.csv")
    keep = ["structure", "mutation", "safe_label", "replicate", "fold_reduction",
            "binding_dg", "binding_dg_sem", "binding_dg_vdw", "binding_dg_vdw_sem",
            "binding_dg_electrostatic", "binding_dg_electrostatic_sem",
            "binding_dg_gb", "binding_dg_gb_sem", "mmgbsa_snapshot_relaxation"]
    _write(OUT / "Supplementary-Table-5.xlsx",
           {"Table 5": panel, "Table 5b agreement": agree,
            "Table 5b phenotype": pheno,
            "Per replicate h_relax": per_rep[[c for c in keep if c in per_rep]]},
           [
               "Supplementary Table 5. Sensitivity of MM/GBSA energies to the snapshot "
               "relaxation protocol.",
               "",
               "Sheet 'Table 5': WT-referenced ddE in kcal/mol for all 19 mutant "
               "genotypes, mean +/- SEM over three replicates, under two protocols.",
               "  h_relax    hydrogens relaxed, heavy atoms restrained (k = 10000 kJ/mol/nm^2)",
               "  full_min   all atoms free; the protocol used for Table 2 of the main text",
               "",
               "Sheet 'Table 5b agreement': Spearman correlation between the two protocols "
               "across the 19 genotypes, per energy component.",
               "Sheet 'Table 5b phenotype': each protocol against log10 fold reduction in "
               "susceptibility, over the 13 genotypes with an established phenotype.",
               "",
               "Sampling is identical across protocols and matches the main-text run: 100 "
               "snapshots spread evenly over the post-equilibration region of each 100 ns "
               "replicate, 0.25 discard fraction, 2000 iteration cap, CPU platform, double "
               "precision.",
               "",
               "Single-point evaluation, with no relaxation, was tested on six genotypes and "
               "is described in the Supplementary Text. It is not tabulated here: its mean is "
               "set by a small number of snapshots carrying steric contacts rather than by the "
               "protocol, so the values are not comparable with the two columns above.",
           ])


def main() -> int:
    table5()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
