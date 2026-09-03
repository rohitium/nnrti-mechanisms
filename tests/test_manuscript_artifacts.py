"""Consistency checks on the manuscript artifacts.

The manuscript's numbers are spread across several files that are generated
separately and pasted into the paper by hand. These tests assert the invariants
that tie them together:

- Table 2's energy components sum to its total.
- Table 2 and Supplementary Table 3 agree with the free energy panel they
  summarise.
- Table 3 and Supplementary Table 4 report the same values.
- Figure 2D takes its phenotype categories from the analysis library rather than
  from a plotting script.

Run with:  PYTHONPATH=src pytest tests -q
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
TABLE2 = REPO / "paper/tables/Table-2-energetics.csv"
TABLE3 = REPO / "paper/tables/Table-3-structural.csv"
SUPP3 = REPO / "paper/submission/Supplementary-Table-3.xlsx"
SUPP4 = REPO / "paper/submission/Supplementary-Table-4.xlsx"
PANEL_DDG = REPO / "results/analysis/fep_pmx/panel_ddg.csv"


def _val(cell: str) -> float:
    """Parse a ``mean ± sem`` cell into its mean."""
    return float(str(cell).split("±")[0].strip().replace("−", "-"))


def _sem(cell: str) -> float:
    return float(str(cell).split("±")[1].strip().replace("<", "").replace("−", "-"))


# --------------------------------------------------------------------------- #
# The package is importable and self-consistent
# --------------------------------------------------------------------------- #

def test_panel_categories_are_disjoint_and_complete():
    from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN

    assert not (SUSCEPTIBLE & RESISTANT)
    assert not (SUSCEPTIBLE & UNCERTAIN)
    assert not (RESISTANT & UNCERTAIN)
    assert len(SUSCEPTIBLE) + len(RESISTANT) + len(UNCERTAIN) == 18


def test_panel_categories_have_one_definition():
    """No CLI may keep its own copy of the phenotype categories."""
    for path in sorted((REPO / "src/nnrti/cli").glob("*.py")):
        src = path.read_text()
        for name in ("SUSCEPTIBLE", "RESISTANT", "UNCERTAIN"):
            assert f"{name} = [" not in src and f"{name} = {{" not in src, (
                f"{path.name} redefines {name}; import it from nnrti.analysis.panel"
            )


# --------------------------------------------------------------------------- #
# Table 2: the components must sum to the total
# --------------------------------------------------------------------------- #

def test_table2_components_sum_to_total():
    rows = list(csv.DictReader(TABLE2.open(encoding="utf-8-sig")))
    assert len(rows) == 18
    cols = list(rows[0])
    total, vdw, elec, gb, sa = cols[2], cols[3], cols[4], cols[5], cols[6]
    for r in rows:
        parts = _val(r[vdw]) + _val(r[elec]) + _val(r[gb]) + _val(r[sa])
        assert parts == pytest.approx(_val(r[total]), abs=0.021), (
            f"{r[cols[1]]}: vdW + elec + GB + SA = {parts:.3f}, "
            f"but the total column says {_val(r[total]):.3f}"
        )


# --------------------------------------------------------------------------- #
# Supplementary Table 3 must reproduce the FEP panel it is summarising
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not PANEL_DDG.exists(), reason="FEP panel not present")
def test_supp_table3_fep_sheet_matches_panel():
    import numpy as np

    fep = pd.read_excel(SUPP3, sheet_name="FEP")
    panel = pd.read_csv(PANEL_DDG).set_index("genotype")
    col = "Leg ddG bind, BAR (kcal/mol)"
    for genotype, d in fep.groupby("Mutation"):
        if genotype not in panel.index:
            continue
        total = sum(g[col].mean() for _, g in d.groupby("Alchemical leg"))
        sem = np.sqrt(sum(g[col].sem() ** 2 for _, g in d.groupby("Alchemical leg")))
        assert total == pytest.approx(panel.loc[genotype, "ddg_bind_kcal"], abs=0.005), genotype
        assert sem == pytest.approx(panel.loc[genotype, "sem_kcal"], abs=0.005), genotype


def test_table2_ddg_matches_the_fep_panel():
    """Table 2's last column must match the free energy panel."""
    panel = pd.read_csv(PANEL_DDG).set_index("genotype")
    rows = list(csv.DictReader(TABLE2.open(encoding="utf-8-sig")))
    ddg_col = list(rows[0])[-1]
    for r in rows:
        g = r[list(rows[0])[1]].replace(" ", "")
        # Table 2 spells compound genotypes in either order; try both.
        key = g if g in panel.index else "+".join(reversed(g.split("+")))
        assert key in panel.index, f"{g} missing from panel_ddg.csv"
        assert _val(r[ddg_col]) == pytest.approx(panel.loc[key, "ddg_bind_kcal"], abs=0.01), g


# --------------------------------------------------------------------------- #
# Table 3 and Supplementary Table 4 are emitted together and must agree
# --------------------------------------------------------------------------- #

def test_table3_matches_supp_table4_summary():
    t3 = pd.read_csv(TABLE3)
    s4 = pd.read_excel(SUPP4, sheet_name="Summary")
    assert "Replicates" not in t3.columns
    assert "F227C" not in set(t3["Genotype"]), "F227C is an alchemical intermediate, not a panel genotype"
    s4 = s4[s4["Genotype"].isin(t3["Genotype"])].set_index("Genotype")
    t3 = t3.set_index("Genotype")
    assert len(t3) == 19
    # Table 3 renames the columns for print, so compare positionally over the
    # nine metric columns, which are last in both tables and in the same order.
    t3_metrics = list(t3.columns)[1:]
    s4_metrics = [c for c in s4.columns if c not in ("Category", "Replicates")]
    assert len(t3_metrics) == len(s4_metrics) == 9
    for gt in t3.index:
        for a_col, b_col in zip(t3_metrics, s4_metrics):
            a, b = str(t3.loc[gt, a_col]), str(s4.loc[gt, b_col])
            assert a == b, f"{gt} / {a_col}: Table 3 says {a!r}, Supp Table 4 says {b!r}"


def test_supp_table4_has_three_replicates_per_genotype():
    detail = pd.read_excel(SUPP4, sheet_name="Per-replicate")
    counts = detail.groupby("Genotype")["Replicate"].nunique()
    assert set(counts) == {3}, f"not 3 replicates everywhere: {counts[counts != 3].to_dict()}"
    assert len(detail) == 60


def test_y188l_has_no_interplanar_angle():
    """Y188L has no Tyr188 ring: the interplanar angle is undefined, not zero."""
    s4 = pd.read_excel(SUPP4, sheet_name="Summary").set_index("Genotype")
    angle = [c for c in s4.columns if "interplanar" in c][0]
    assert pd.isna(s4.loc["Y188L", angle]) or s4.loc["Y188L", angle] in ("", "—")
    assert not pd.isna(s4.loc["WT", angle])


# --------------------------------------------------------------------------- #
# Fold-change values are never hardcoded
# --------------------------------------------------------------------------- #

def test_fold_changes_come_from_the_spreadsheet():
    xlsx = REPO / "data/DRM-susceptibilities.csv.xlsx"
    assert xlsx.exists(), "the authoritative susceptibility source is missing"
    values = REPO / "results/analysis/susceptibility/tables/dor_susceptibility_values.csv"
    if not values.exists():
        pytest.skip("run nnrti.cli.plot_dor_susceptibility_bars first")
    df = pd.read_csv(values)
    assert len(df) >= 18
