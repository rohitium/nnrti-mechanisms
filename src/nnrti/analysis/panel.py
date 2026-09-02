"""The genotype panel: phenotype categories and genotype-name normalisation.

These definitions are used by the manuscript figures and tables, so they live in
the analysis library rather than inside a plotting CLI. They were previously
declared in ``plot_wt_referenced_occupancy_tick_lines`` and
``plot_triplet_contact_story``, which meant that Figure 2D imported two
exploratory scripts in order to learn which genotypes are Resistant -- a
dependency that would have broken the moment those scripts were archived.

The three categories are the panel definition given in Results ("DOR phenotypic
susceptibility landscape") and must stay in step with Table 1. Fold-change
values are never hardcoded anywhere; they are read from
``data/DRM-susceptibilities.csv.xlsx`` (see :mod:`nnrti.analysis.susceptibility`).
"""

from __future__ import annotations

import re

#: Well-established low-impact mutations: little if any reduction in DOR susceptibility.
SUSCEPTIBLE = {"V106I", "K103N", "Y181C", "G190A"}

#: Canonical DOR-associated resistance mutations and patterns.
RESISTANT = {
    "V106A",
    "Y188L",
    "Y318F",
    "A98G+F227C",
    "V106A+F227L",
    "V106A+L234I",
    "V106A+P225H",
    "V106I+F227C",
    "K103N+M230L",
}

#: Limited or inconsistent reductions in susceptibility.
UNCERTAIN = {"L100I+K103N", "K103N+P225H", "V106M", "G190E", "G190S"}

#: Plot/table ordering for the three categories.
CATEGORY_ORDER = {"Susceptible": 0, "Uncertain": 1, "Resistant": 2}


def normalize_mutation_token(text: str) -> str:
    """Canonicalise a genotype label: upper case, no spaces, ``+`` separators.

    Accepts the comma- and space-separated spellings that appear in the
    susceptibility spreadsheet and the manifests.
    """
    t = str(text).strip().upper()
    if not t or t == "NAN":
        return ""
    t = t.replace(" ", "").replace(",", "+")
    return re.sub(r"\++", "+", t)


def category_of(mutation: str) -> str:
    """Return ``"Susceptible"``, ``"Resistant"``, ``"Uncertain"`` or ``"WT"``.

    Raises ``KeyError`` for anything outside the simulated panel, so a typo in a
    genotype name fails loudly instead of silently dropping out of a figure.
    """
    m = normalize_mutation_token(mutation)
    if m in ("WT", "WILDTYPE", "WILD-TYPE"):
        return "WT"
    if m in SUSCEPTIBLE:
        return "Susceptible"
    if m in RESISTANT:
        return "Resistant"
    if m in UNCERTAIN:
        return "Uncertain"
    raise KeyError(f"{mutation!r} is not part of the simulated genotype panel")
