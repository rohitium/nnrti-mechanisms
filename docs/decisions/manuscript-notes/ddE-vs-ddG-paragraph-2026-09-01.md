# Discussion paragraphs: why ∆∆E_Total and ∆∆G_bind disagree

Numbers verified 2026-09-01 against `panel_ddg.csv`, `ddg_full.csv` and the
per-leg `analysis.json` files.

| statistic | value |
|---|---|
| ∆∆G_bind vs ∆∆E_Total, 18 genotypes | R² = 0.183, p = 0.077 |
| sign agreement | 12 of 18 |
| opposite sign with **both** values resolved | A98G+F227C, G190E, K103N, K103N+M230L |
| opposite sign, one value unresolved | Y181C, L100I+K103N |

FEP leg decomposition (∆∆G_bind = ∆G_holo − ∆G_apo):

| leg | ∆G holo | ∆G apo | ∆∆G |
|---|---:|---:|---:|
| wt → G190E | −53.84 ± 1.22 | −52.82 ± 1.51 | **−1.02** |
| wt → Y188L | −0.40 ± 0.35 | −4.92 ± 0.32 | **+4.52** |
| K103N → K103N+M230L | −14.09 ± 1.94 | −14.14 ± 2.21 | **+0.05** |
| F227C → A98G+F227C | −8.97 ± 0.27 | −7.75 ± 0.19 | **−1.22** |

---

## Suggested text (two paragraphs)

> Across the 18 genotypes the two computed quantities are only weakly related
> (R² = 0.18, p = 0.08) and agree in sign for 12 of them. In four
> genotypes — K103N, K103N+M230L, G190E and A98G+F227C — the two are opposite in
> sign with both estimates resolved beyond their standard errors, so the
> disagreement cannot be attributed to replicate noise. This is expected, because
> the two quantities are not estimates of the same thing. ∆∆E_Total is an
> interaction energy evaluated on the bound complex alone, whereas ∆∆G_bind is a
> difference between the bound and unbound states. A substitution that
> destabilises the protein equally whether or not the drug is present contributes
> nothing to ∆∆G_bind however much it perturbs the interface, and conversely a
> substitution that reshapes the unliganded pocket can change ∆∆G_bind without
> altering the bound-state interaction energy at all.
>
> G190E illustrates this most clearly. Decomposing the alchemical calculation
> into its two legs, introducing the glutamate charge costs 53.84 ± 1.22 kcal/mol
> in the holo state and 52.82 ± 1.51 kcal/mol in the apo state; the binding term
> is the 1.02 kcal/mol residual, under 2% of either leg. Almost the entire
> energetic consequence of the mutation is an intrinsic charging and desolvation
> cost that the enzyme pays in both states and which therefore cancels in the
> double difference. MM/GBSA, which never samples the apo state, reports only the
> bound-state consequence — ∆∆E_Total = 1.94 ± 0.40 kcal/mol, dominated by a
> polar solvation penalty of 2.50 ± 0.32 kcal/mol as the new charge is buried at
> the interface. Y188L is the contrasting case in which the legs genuinely differ:
> removing the tyrosine costs 4.52 kcal/mol more in the bound state than in the
> unbound one, because the residue is expendable in an empty pocket but not when
> it is stacking against the drug, and here the two methods agree. A further
> consideration applies to the charged and polar genotypes, where ∆∆E_Total is a
> small residual of two large and opposing terms: for K103N+M230L a Coulombic
> gain of −2.74 ± 0.56 kcal/mol survives a generalised-Born desolvation penalty of
> +2.33 ± 0.57 kcal/mol, so the sign of the total is set by an implicit-solvent
> approximation in precisely the regime where it is least reliable. A98G+F227C is
> the one discrepancy not of this kind: its electrostatic and solvation terms are
> both negligible and its positive ∆∆E_Total is driven by van der Waals packing
> (1.81 ± 0.27 kcal/mol), yet the alchemical calculation finds the substitution
> more favourable in the bound state than in the unbound one, indicating that the
> mutation's effect on the protein's own conformational preferences — which the
> end-point calculation does not capture — outweighs the interfacial packing loss.

---

## Why this is worth saying plainly

The temptation is to present the disagreement as a limitation. It is more useful
as a result: the genotypes where the two methods agree (Y188L, the V106A series)
are exactly those where the mutation acts locally at the bound interface, and
the ones where they diverge are those where it does not. The pattern is
diagnostic rather than noisy, and it identifies which resistance mechanisms are
interfacial and which are not — which is the paper's central question.

Note also that this weakens a claim the draft currently makes: at R² = 0.18,
p = 0.077 the correlation between the two computed quantities is **no longer
statistically significant** and should not be described as such. See
`paper/CHANGES_2026-09-01_g190e_resolved.md`.
