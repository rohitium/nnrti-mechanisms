# JCIM abstract — unstructured rewrite

The 09-02 draft's abstract is **289 words under four subheadings**
(Background / Methods / Results / Conclusion). JCIM Articles take a single
unstructured paragraph, so this has to be restructured regardless of how G190E
resolves — which is why it is worth doing now.

It also carries two stale statistics: `R² = 0.02, p = 0.60` for ∆∆E<sub>Total</sub>
against susceptibility (now **0.06, p = 0.31** over all 18; **0.13, p = 0.23**
over the 13 established), and `R² = 0.26, p = 0.03` for the two computed
quantities against each other (now **0.22, p = 0.048**). Both are corrected below.

---

## Option A — tight (108 words)

Closest to the "3–4 sentences" reading of the guidelines. Drops the panel
composition and the methodological detail.

> Doravirine (DOR) retains activity against several common NNRTI resistance
> mutations, but public phenotypic susceptibility data remain sparse for many
> DOR-associated genotypes. We performed all-atom explicit-solvent molecular
> dynamics simulations of HIV-1 reverse transcriptase bound to DOR for the
> wild-type enzyme and 18 clinically observed genotypes, and estimated
> mutation-induced changes in binding using nonequilibrium alchemical free-energy
> perturbation (∆∆G<sub>bind</sub>) and wild-type-referenced MM/GBSA interaction
> energies (∆∆E<sub>Total</sub>). ∆∆G<sub>bind</sub> separated resistant from
> susceptible genotypes with a Matthews correlation coefficient of 0.50 and
> correlated only weakly with measured fold-change in susceptibility
> (R² = 0.26, p = 0.07), while MM/GBSA misclassified three of four susceptible
> mutations. Binding energetics alone therefore recover broad resistance trends
> but do not quantitatively reproduce phenotypic drug susceptibility.

---

## Option B — full (196 words), recommended

Keeps the mechanistic result, which is the part reviewers will want and which
Option A loses entirely. Still one paragraph, comfortably within typical JCIM
Article abstract length.

> Doravirine (DOR) is a second-generation non-nucleoside reverse transcriptase
> inhibitor that retains activity against several common NNRTI resistance
> mutations, but public phenotypic susceptibility data remain limited for many
> DOR-associated genotypes. To assess whether physics-based simulation could help
> fill this gap, we performed all-atom explicit-solvent molecular dynamics
> simulations of HIV-1 reverse transcriptase (RT) bound to DOR for the wild-type
> enzyme and a curated panel of 18 clinically observed genotypes spanning
> susceptible, resistant and uncertain phenotypes, and estimated mutation-induced
> changes in binding by nonequilibrium alchemical free-energy perturbation
> (∆∆G<sub>bind</sub>) and by wild-type-referenced MM/GBSA interaction energies
> (∆∆E<sub>Total</sub>). ∆∆G<sub>bind</sub> identified 7 of 9 resistant and 3 of 4
> susceptible genotypes correctly (Matthews correlation coefficient 0.50; ROC AUC
> 0.78), whereas MM/GBSA detected resistance equally often but misclassified three
> of four susceptible mutations; neither quantity correlated strongly with
> measured fold-change in susceptibility (R² = 0.26, p = 0.07 and R² = 0.13,
> p = 0.23 across the 13 genotypes with established phenotypes). The equilibrium
> trajectories nonetheless resolved distinct structural mechanisms, including loss
> of aromatic packing in Y188L and displacement of the drug's anchored pyridinone
> ring in V106A-containing genotypes. Relative binding free-energy calculations
> thus recover broad resistance-associated trends, but binding energetics alone
> are insufficient to reproduce phenotypic drug susceptibility quantitatively.

---

## Notes

- Both versions state the classification result as an **effect size** (MCC, AUC)
  rather than a significance claim, consistent with the Methods position that at
  n = 13 nothing reaches conventional significance.
- Neither mentions G190E, so neither needs revisiting when Lever C lands.
- The Conclusion sentence about "more comprehensive computational models
  incorporating conformational dynamics, enzyme function, and viral fitness" is
  dropped from both for length. It is a good sentence and it survives in the
  Discussion; if you want it in the abstract, Option B has the most room.
- **Unverified:** the "3–4 sentences" figure comes from a reading of the JCIM
  author guidelines earlier in this project. Worth re-checking against the live
  guidelines before submission, since it drives the choice between A and B.
