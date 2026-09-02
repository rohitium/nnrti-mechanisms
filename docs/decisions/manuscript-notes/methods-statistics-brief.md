# Methods — Statistical analysis (brief)

Replaces the longer draft in `methods-statistical-analysis.md`, with the
classification metrics removed.

---

## Statistical analysis

> Uncertainties on ∆∆G<sub>bind</sub> and ∆∆E<sub>Total</sub> are standard
> errors of the mean across three independent replicate simulations. For
> genotypes constructed as sequential single-residue legs, the per-leg errors
> were propagated in quadrature. For the MM/GBSA components, each mutant
> replicate was referenced to the mean of the three wild-type replicates; the
> uncertainty of that reference is a common offset shared by every genotype and
> is therefore excluded from the per-genotype standard error, which is what
> governs comparisons between genotypes. Associations between the computed
> energies and phenotype were assessed by ordinary least-squares regression
> against log<sub>10</sub>-transformed fold-change in susceptibility, reported as
> the coefficient of determination R² with the two-sided p-value for a non-zero
> slope (`scipy.stats.linregress`); fold-change was log-transformed because
> susceptibility is measured as a ratio spanning two orders of magnitude, and
> because ∆∆G = −RT ln(fold) under a pure binding model makes the logarithm the
> scale on which the relationship is expected to be linear. Rank association was
> additionally checked with Spearman's ρ, which does not assume linearity, and
> the correlation between the two computed quantities was assessed the same way.
> Given the size of the panel, these statistics are reported as descriptive
> effect sizes rather than as hypothesis tests, and no correction for multiple
> comparisons was applied.

---

## Notes

- ~200 words, one paragraph.
- Covers only: replicate SEMs, quadrature propagation, the WT-referencing
  convention for MM/GBSA, the log-transform and its justification, OLS R²/p,
  Spearman, and the descriptive-not-inferential framing.
- Omits everything about sensitivity, specificity, MCC, ROC AUC, bootstrap
  confidence intervals and permutation tests, per request. If the classification
  result is later restored to the Results, the corresponding Methods text is
  preserved in `methods-statistical-analysis.md`.
- The 0.5 kcal/mol strong/weak reporting threshold is deliberately left where it
  currently sits in the Results, since it is a physical constant (k<sub>B</sub>T
  at 300 K) used as a reporting convention rather than a fitted statistic.
