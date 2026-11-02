# JCIM ci-2026-030865 — Revision 1 status

"Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase".
Satija, Tao, Shafer. **Due 10-Nov-2026.** Last updated 02-Nov-2026.

This is the only live note in `revisions/`. The eleven earlier memos are in
`_archive/` with nothing removed; §9 says which two still hold text worth
re-reading.

---

## 1. Documents and house rules

| file | state |
|---|---|
| `Manuscript-Revision-1.docx` | working copy, additions in **red**. Discussion rework pending (§4) |
| `Supplementary-Text-Revision-1.docx` | apo numbers corrected 02-Nov |
| `Response to Reviewers.docx` | apo numbers corrected 02-Nov |
| `Supplementary-Table-*.xlsx` | tables 5–9. Table 6 needs the four-genotype AM1-BCC rebuild (§8) |
| `Supporting-Information-Cover-Page.docx`, `TOC-Graphic.{pdf,tif,png}` | done |

Clean copies sit in `_stale-clean-copies/` and get regenerated at the end by
stripping the red colour, so they are identical by construction. Do not hand-edit them.

- **The author writes the prose.** Supply numbers and technical content.
- **Edit in place, never rebuild a `.docx`.** Backups: `_backups/20261102-0246/`.
- **Citations are live Zotero fields.** Never renumber by hand. The Discussion
  alone holds 26 of them, so paragraph moves go through Word, not a script.
- **Response format:** "To address this comment, we have added/modified this
  paragraph in the \<section\> section(s):" then the manuscript text verbatim.
- Spread each edit across the sections it belongs in. Never introduce a result
  for the first time in the Discussion.

## 2. Reviewer coverage

Reviewers 2 and 3 are complete. Reviewer 1 is complete except for the Discussion
rework. The 7Z2G crystal-representativeness section and the perturbed-pose
ladder are both folded into the documents (5–6 mentions in each).

The one scientific request not acted on is **R1.2a, independent starting
poses**: replicates differ only by 0.1 Å jitter, and the response says so
plainly. The pose ladder tests pose sensitivity, not basin exploration, so it
does not cover this. If there is a second round, this is where it comes from.

## 3. Decisions waiting on you

**(a) The cancellation sentence in Discussion paragraph 595.** It says G190E and
A98G+F227C distort the interface "but comparable distortion is present in the
apo simulations, so the penalty largely cancels in the thermodynamic cycle."
The 100 ns apo data does not support that. Counting features whose 95% interval
excludes zero:

| | holo | apo |
|---|---|---|
| **G190E** | **8** resolved: pocket +55.8 ± 12.1 Å³, Val179 +2.84 ± 0.44 Å, RMSF in pocket lining, primer grip, β6, fingers, palm | **0** |
| **A98G+F227C** | 2: Tyr181 χ1 −75.8 ± 12.0°, primer-grip RMSF +0.09 ± 0.03 Å | 2, but **different features**: palm helix content, connection RMSF |

For G190E the distortion is bound-state only, so there is nothing to cancel
against. Figure 4 agrees from the other side: G190E's induced fit differs from
wild type by **+71 ± 26 Å³**, a discordant change, which is the case that should
*not* cancel. So G190E distorts the bound complex more than anything else in the
panel, it does not cancel, and ΔΔGbind is still favourable at −1.02 ± 0.38.
Three ways out, no recommendation from me:

1. Drop the cancellation clause and state the puzzle plainly.
2. Keep cancellation as a hypothesis, flagged as unsupported by the apo data at
   three replicates.
3. Point at the G190E confounds already in paragraph 593 (replication capacity,
   single-isolate phenotype) and treat ΔΔGbind as untestable for this genotype.

The same sentence is quoted in the R1.3 and R2.2 responses, so the choice
propagates to both.

**(b) The dropped-metric sentence in the Figure 4 caption.** The four metrics
removed are the subdomain and catalytic-motif fluctuations, which is exactly the
slow loop remodelling R1 asks about for the P225 and M230 genotypes. Saying so
in the caption answers that part of the comment honestly, but it is also a
concession and may read better in Results with a pointer in the caption.

**(c) A Methods half-sentence on the 95% criterion.** Figure 3 needs 2.78 SEM,
Figure 4 needs 2.31, because of the degrees of freedom (two groups of three
against four). A reader applying a mental factor of 2 will think more cells are
resolved than are.

## 4. The Discussion rework

Your structure: three opening paragraphs as they stand, then one paragraph per
mechanism, then one paragraph per limitation. Real paragraph indices in
`Manuscript-Revision-1.docx` (Discussion is 588–597, 26 Zotero fields):

**Opening, unchanged.** 588 context · 589 approach · 590 headline result.

**Mechanisms, limitations stripped out.**

| ¶ | topic | change |
|---|---|---|
| 591 | K103N and the backbone anchor | none |
| 592 | Tyr181/Tyr188 aromatic packing | none |
| 593 | the G190 series | delete the last two sentences (replication capacity; "might affect the primer grip, which is not modeled here") — they move to the limitations. "pocket enlarged by ~25%" reads as a new result; cite Table 3 instead |
| 594 | the V106 series | delete the last 90 words, "We quantified the uncertainty…" to "…which they do not", and replace with the sentence below |
| 595 | ΔΔETotal against ΔΔGbind | delete the final implicit-solvent sentence (moves to limitations) and settle §3(a) |
| 596 | K103N+M230L | delete the last two sentences on the binary complex (move to limitations) |

Replacement for the end of 594:

> The displacement is well resolved (−1.35 Å [−1.73, −0.97] for V106A relative
> to wild type, Hedges' g = −9.5); the change in pyridinone contact count is
> not, and the argument rests on the displacement.

**Limitations, each stated once.**

1. **Sampling.** Three replicates from one seeded pose sample one basin;
   convergence within a basin is not evidence the basin is the only relevant
   one; the unliganded equilibrium is slower than 100 ns. Add: of 1,287
   genotype-by-feature comparisons, 66 reach 95%, which is the null
   expectation, so most structural differences in the panel are unresolved at
   three replicates.
2. **A simplified model system.** The binary RT–DOR complex omits
   template/primer and dNTP, so primer-grip and assembly mechanisms are
   structurally absent (merge from 593 and 596, where it is currently stated
   twice); MM/GBSA uses implicit solvent, which may overestimate desolvation
   for a buried glutamate (from 595); the ligand charge model shifts ΔΔETotal
   by up to 1.33 kcal/mol without changing its sign.
3. **Mechanisms we cannot reach.** 597 as written, plus the
   replication-capacity sentence moved out of 593.

The 90-word statistical passage leaving 594 belongs in Results or the SI.

All 17 energetic and structural numbers quoted in the Discussion were checked
against Tables 2 and 3 and are correct. Only the apo-derived claims needed work.

## 5. Figure 3 and Figure 4 captions

Written against the figures as they stand: 18 genotypes, Figure 3 with 7 metrics,
Figure 4 with 5. Both drop the Tyr188 interplanar angle, which no genotype
resolves; reinstating it changes the dropped-metric sentence.

**Figure 3. Changes in RT–DOR contacts across the genotype panel, referenced to
wild type.** Each cell gives the shift of a structural metric from its wild-type
value, in that metric's own units, above the standard error of the shift. Errors
combine the genotype's and wild type's standard errors over three independent
100 ns replicates in quadrature, since wild type is itself a three-replicate
mean. Colour is scaled within each column, so the darkest cell in a column is
the largest shift for that metric; it carries no meaning across columns, which
are in ångström, degrees and counts. Genotypes are grouped by the resistance
category of Table 1 and ordered by fold reduction in doravirine susceptibility
within each group. Metrics for which no genotype's 95% confidence interval
excluded zero are not shown; on this criterion the Tyr188 / chlorocyanophenyl
interplanar angle was removed. Of the 126 genotype–metric comparisons plotted,
21 have a 95% interval excluding zero. Per-comparison effect sizes, confidence
intervals and p values for every metric, including those not plotted, are in
Supplementary Table X.

**Figure 4. Mutation-induced changes in RT's own conformation, as the difference
between the bound and unliganded states, referenced to wild type.** Each cell
gives (holo − apo) for a genotype minus (holo − apo) for wild type: the extent
to which binding doravirine reshapes the enzyme, and how far that reshaping
differs from wild type's. A value near zero means the mutation moves the bound
and unliganded states by the same amount, concordant remodelling, which cancels
in a thermodynamic cycle and is therefore invisible to ΔΔGbind by construction.
Errors propagate the standard errors of all four three-replicate means in
quadrature. Colour, grouping and ordering are as in Figure 3. Metrics for which
no genotype's 95% confidence interval excluded zero are not shown; on this
criterion the Cα RMSF of the YMDD motif, the p66 thumb and the p66 connection
subdomains, and the β9–β10 strand content were all removed, and we report these
as unresolved at three replicates rather than as unchanged. Of the 90
genotype–metric comparisons plotted, 11 have a 95% interval excluding zero.

## 6. The numbers that matter now

**Effect sizes, one method for every feature.**
`nnrti.cli.compute_structural_effect_sizes` →
`results/analysis/structural_metrics/tables/structural_effect_sizes.csv`:
shift, sem, Welch t/df/p, ci95, Hedges' g, for 1,287 comparisons in both states.
At n = 3, 2× SEM is p = 0.12; the 95% interval needs 2.78× SEM for a
mutant-vs-WT difference and 2.31× for Figure 4's difference of differences.

- V106A toward Ser105: **−1.35 Å [−1.73, −0.97]**, p = 0.004, g = −9.5.
- K103N side chain to DOR: −3.31 Å [−3.89, −2.73], p = 0.001, g = −14.2.
- G190E clears 95% on **eight** holo features.
- V106M clears it on **none**, so "within error on every metric we examined" stands.
- 66 of 1,287 reach 95%, the null expectation for that many tests.

**Apo pocket and opening cost** (3 × 100 ns, all 20 genotypes;
`pocket_opening_cost_100ns.csv`). The holo arm reproduces the published value to
0.6 Å³ (227.4 against 228), so the recompute is trustworthy. Wild-type apo is
**76 ± 9 Å³** and opening costs **147.5 Å³**, not 137. With errors propagated,
only **G190E (+71 ± 26 Å³)** differs from wild type; A98G+F227C is +22 ± 19 and
K103N+M230L is +2 ± 28, so the claim that all three mismatched genotypes pay
more does not survive, and a large nominal cost is not peculiar to them.

Three genotypes differ from wild type in apo volume, all V106A-containing and
all **more open**: V106A +52 ± 11 Å³ (p = 0.02), V106A+L234I, V106A+P225H. V106A
pre-organises the unliganded pocket toward the bound geometry, which is a
measured instance of what R1.2 asks about. Y181C is a milder second case
(apo 103.5 Å³, χ1 at +40° against WT's −68°).

One replicate is still excluded, G190E apo rep_02, which loses 158 of 200 frames
to the pocket-Cα Rg > 10 Å filter (local unfolding, not a pocket). That is now
the **only** exclusion; Methods still describes two. The apo cavity is not
re-solvated before measurement, so apo volumes are biased low and openings are
upper bounds. This is the one Methods caveat worth keeping. Delete the other,
"apo replicates differ in length", which is no longer true.

**DOR in water: binding stops the clock rather than narrowing the
distribution.** Across 60 bound replicates × 4 torsions × 100 ns, 24 µs of
cumulative torsion-time, DOR leaves its bound rotamer **three times**. Free in
water each torsion leaves several times per 100 ns at matched frame spacing, and
3–50 times at full resolution. Circular SD is 8–12° bound against 68–116° free.
The largest mutational effect on ligand rigidity is a quarter of the
bound-to-free gap, so **every genotype holds DOR in essentially the same
internal conformation** and resistance does not work by letting the ligand relax.
ΔG_preorg is 3–4 kcal/mol; do not quote a decimal, because τ3 crosses its
barrier only 3–6 times per 100 ns and the joint population multiplies through
it. ΔG_preorg is a free-ligand property, identical for every genotype, so it
cancels exactly in a mutant-minus-WT difference and cannot explain why ΔΔGbind
misses G190E.

**Charge model.** AM1-BCC now covers WT, V106A, Y188L, G190E and K103N+M230L at
3 replicates each (`charge_model_comparison.csv`, finished 02-Nov 02:59). The
largest ΔΔETotal shift is **1.33 kcal/mol (G190E)**, and all four signs are
preserved. The old "≤0.81 kcal/mol" in the coverage audit is superseded; the
documents already say 1.33.

## 7. Simulations

All 60 apo replicates are at 100 ns except **v106a rep_02**, job 50786462,
running on sh02-13n13 and at 44,952,000 steps (89.9 ns) as of 09:40, about 5 h
from 100 ns. It resumed correctly from the canonical tree
(`$SCRATCH/nnrti-mechanisms-git`) at step 43,080,000. When it lands, pull it and
regenerate the apo chain so V106A reaches n = 3; the V106A-class conclusion does
not depend on it, since three combination genotypes carry it at n = 3 each.

## 8. Open work

- Supplementary Figure 5 caption still says "at most 0.25 Å … for Y181C and
  G190S". After the p66-anchored re-imaging fix the maximum is 0.193 Å, with
  G190S top and Y181C down at 0.116.
- Supplementary Table 6 workbook needs rebuilding with the four-genotype AM1-BCC
  data (`ops/docx/build_si_workbooks.py`).
- Supplementary Table 4 workbook: the structural metrics table exists as CSV but
  `build_supplementary_table_4.py` is still the holo-only builder.
- The apo under-solvation caveat still needs a home in the Supplementary Text.
- `mutation_component_summary.csv` is stale and does not reconcile with the
  manuscript. Use `binding_energy/tables/ddg_full.csv`.
- Five `src/nnrti/cli/*.py` files are untracked.
- **Author actions:** ORCID iDs, the funding field (NIAID R24AI136618) on the
  submission form, and an English-language pass over the pre-existing text. The
  ~1,100 added words were audited; the original text has not been.

## 9. Where the detail went

`_archive/` holds the eleven superseded memos. Two are still worth opening:

- `R1-crystal-representativeness.md` — the full 4NCG vs 7Z2G evidence and the
  manuscript-ready SI text. The text is already in the documents; the file is
  the record of how the numbers were obtained.
- `deleted-table-notes.md` — verbatim text deleted from the supplementary
  workbook Notes tabs, including the method and limitation lines that have
  nowhere to live except the Supplementary Text.

Technical traps that cost time live in `CLAUDE.md`, not here.
