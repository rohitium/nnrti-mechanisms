# Manuscript changes — G190E resolved, tier artifact removed (2026-09-01)

## 1. G190E: ∆∆G_bind = **−1.02 ± 0.38** kcal/mol

Was `+2.00 ± 1.77` in the draft, then `+0.99 ± 1.63` after the 20 ns
equilibration campaign. The Lever C repair — re-running holo replicate 2's λ=1
endpoint ensemble alone — resolves it.

| | 5 ns equil | 20 ns equil | **after Lever C** |
|---|---:|---:|---:|
| ∆∆G_bind | +2.00 | +0.99 | **−1.02** |
| SEM | 1.77 | 1.63 | **0.38** |
| σ_DDG | 3.06 | 2.82 | **0.66** |
| units failing BAR–Jarzynski | 1 | 1 | **0** |

**Two things change in the text, not one.**

The **magnitude** is now the most precise in the panel — SEM 0.38, comfortably
inside the 0.5 kcal/mol standard. The old sentence explaining G190E's weak
classification "because of large replicate variance" is obsolete and must go.

The **sign** flips. G190E slightly *favours* DOR binding against an 18-fold
measured loss of susceptibility, so binding energetics do not explain its
resistance at all. This strengthens the paper's central limitation argument and
now rests on the panel's best-determined value rather than its worst.

Text to replace, Results:

> "…in case of G190E because of large replicate variance ∆∆G<sub>bind</sub> =
> 2.00 ± 1.77 kcal/mol…"

becomes

> …and in the case of G190E the computed change is slightly favourable
> (∆∆G<sub>bind</sub> = −1.02 ± 0.38 kcal/mol) despite an 18-fold reduction in
> susceptibility, indicating that the resistance of this genotype is not
> mediated by a loss of DOR binding affinity.

Table 2's ∆∆G_bind cell for G190E: **−1.02 ± 0.38**.

## 2. Updated statistics — recomputed, all of them

| quantity | draft | **current** |
|---|---|---|
| ∆∆G<sub>bind</sub> vs ∆∆E<sub>Total</sub> (18) | R² = 0.25, p = 0.03 | **R² = 0.18, p = 0.077** |
| ∆∆G<sub>bind</sub> vs log₁₀(fold), all 18 | R² = 0.09, p = 0.21 | **R² = 0.079, p = 0.258** |
| ∆∆E<sub>Total</sub> vs log₁₀(fold), all 18 | R² = 0.02, p = 0.55–0.60 | **R² = 0.064, p = 0.310** |
| ∆∆G<sub>bind</sub> vs log₁₀(fold), 13 established | R² = 0.26, p = 0.07 | **R² = 0.261, p = 0.075** (unchanged) |
| ∆∆E<sub>Total</sub> vs log₁₀(fold), 13 established | R² = 0.05, p = 0.45 | **R² = 0.126, p = 0.234** |
| Spearman ρ, ∆∆G vs fold (18) | 0.351 | **0.295, p = 0.234** |

**Note the correlation between the two computed quantities weakened**
(R² 0.22 → 0.18, and p crossed 0.05 to 0.077). G190E moving from +0.99 to −1.02
while its ∆∆E_Total stays at +1.94 pulls that fit apart. The Discussion sentence
"∆∆E<sub>Total</sub> and ∆∆G<sub>bind</sub> showed a modest correlation
(R² = 0.25, p = 0.03)" must now read **R² = 0.18, p = 0.08** and should no
longer be described as statistically significant.

Classification performance is **unchanged** (FEP MCC 0.50, AUC 0.78; MM/GBSA MCC
0.03, AUC 0.69) — G190E is an Uncertain-phenotype genotype and is excluded from
the 13 scored genotypes.

Mean SEM: FEP **0.656** (was 0.722), MM/GBSA 0.518. Eleven of 19 FEP legs are
still above 0.5 kcal/mol.

## 3. The reporting "tier" has been deleted

`combine_neq.py` classified each genotype `main_text` / `show` / `omit_main`
from its SEM, and wrote `panel_discussion_tiers.csv`. That was campaign triage
and never belonged in the output: every table already prints the SEM beside the
value, and the labels went stale whenever a leg was re-run — G190E was still
marked `omit_main` at a SEM of 0.38, and K103N at 0.23.

Removed: `OMIT_MAIN_TEXT`, `discussion_tier()`, `write_discussion_tiers()`, the
tier column in the console summary, and `panel_discussion_tiers.csv` itself.
`build_table_2.py` and `build_supplementary_table_3.py` now read
`panel_ddg.csv`, which carries the same values plus `n_reps`.

## 4. Regenerated

| artifact | file |
|---|---|
| Table 2 | `manuscript/Table-2-energetics.csv` |
| Supplementary Table 3 | `manuscript/Supplementary-Table-3.xlsx` |
| FEP panel vs experiment | `results/analysis/fep_pmx/panel_ddg_vs_experiment.png` |
| by resistance category | `results/analysis/fep_pmx/panel_ddg_vs_experiment_by_category.png` |
| category subset stats | `results/analysis/fep_pmx/panel_category_subset_stats.csv` |
| Crooks overlap QC | `results/analysis/fep_pmx/panel_crooks_overlap.png` |
| per-leg protocol figures | `results/analysis/fep_pmx/protocol/` (18 legs) |
| classification metrics | `results/analysis/classification_performance/` |

MM/GBSA is untouched — no MM/GBSA quantity depends on the FEP result.

## 5. Still to do

- Apply the G190E sentence and Table 2 cell to the draft (Word, by hand).
- The Discussion correlation sentence must drop "statistically significant".
- `manuscript/abstract-jcim-unstructured.md` Option B quotes R² = 0.26/p = 0.07
  for the 13-genotype fit — that one is unchanged, so the abstract stands.
