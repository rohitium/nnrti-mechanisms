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
| `Supplementary-Table-1..10.xlsx` | the complete set, all rebuilt or verified 02-Nov. See the table below |
| `Supporting-Information-Cover-Page.docx`, `TOC-Graphic.{pdf,tif,png}` | done |

### The Supplementary Table set

All ten now live in `revisions/`, so there is one place to submit from.
Previously Tables 1–4 existed only as the submitted copies and 5–9 only here.

| # | contents | state on 02-Nov |
|---|---|---|
| 1 | isolate-level DOR susceptibility | unchanged. All 18 fold values re-verified against `data/DRM-susceptibilities.csv.xlsx`; `manifests/md_manifest.csv` agrees too |
| 2 | RT p66 sequences | unchanged |
| 3 | per-replicate MM/GBSA and per-leg FEP | unchanged, and the builder reproduces it to 0.0 on all 138 rows. This is where the §3(a) cancellation numbers live |
| 4 | per-replicate structural observables | **rebuilt in both states.** The bound half reproduces the submitted workbook on all 180 summary cells; the unliganded half and the χ1, RMSF and secondary-structure metrics are new |
| 5 | MM/GBSA relaxation sensitivity | rebuilt, unchanged |
| 6 | ligand charge model | **rebuilt with four genotypes plus WT.** Was two. Largest ΔΔE_Total shift 1.33 kcal/mol (G190E), all signs preserved |
| 7 | apo and holo NNIBP volume | **rebuilt from the 100 ns apo ensembles.** The previous build predated them and disagreed with the documents |
| 8 | 4NCG vs 7Z2G | current against its sources (spot-checked: DOR 0.619 Å, pocket Cα 0.490 Å) |
| 9 | Tyr181/Tyr188 across four NNRTIs | current against its sources |
| 10 | **new.** Effect size and 95% interval behind every Figure 3 and Figure 4 cell | built. This is the "Supplementary Table X" the figure captions referred to |

Builders: `ops/docx/build_si_workbooks.py` (5–7),
`nnrti.cli.build_si_structural_tables` (4 and 10),
`nnrti.cli.build_supplementary_table_3` (3). **Tables 8 and 9 have no builder** —
they were written ad hoc, so they cannot be regenerated from a command. Worth
fixing before any of their sources change.

Previous copies: `_backups/20261102-xlsx/`.

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

## 3. Decisions taken, and what is still yours

All three of the decisions this file previously left open have been taken, with
the reasoning recorded here so they can be reversed.

**(a) The cancellation sentence is kept and corrected, not dropped.** It was
right; it was only ambiguous. ΔΔGbind is the difference of the bound and
unliganded alchemical legs, so anything the mutation does to both cancels by
construction, and Supplementary Table 3 already carried the numbers: G190E
costs −53.8 kcal/mol in the bound leg and −52.8 in the unliganded leg, leaving
−1.02. Across the panel the residual is a median 9% of the leg magnitude, and
for the three mismatched genotypes 0.3–2%. The phrase "the apo simulations"
could only mean the alchemical leg when it was written, but the revision adds
3 × 100 ns equilibrium apo ensembles and a Figure 4 built on holo − apo, so the
sentence now names the leg and quantifies it.

**(b) Main-text Table 3 stays.** Figures 3 and 4 were added as Reviewer 1
comment 4 asked, not substituted for the table, and the manuscript already
referenced both. Keeping it leaves every existing `(Table 3)` pointer correct.

**(c) The dropped-metric sentence stays in the Figure 4 caption,** and Methods
now states the 95% criterion (2.78 standard errors for Figure 3, 2.31 for
Figure 4), so a reader applying a mental factor of 2 is not misled.

**Still yours:** whether to say more about the primer-grip result in §6, and
whether the G190E replication-capacity specifics belong in the paper — see §8.

## 4. The Discussion rework — applied

The Discussion now separates mechanism from limitation, so no caveat is stated
more than once. Opening unchanged (context, approach, headline result), then
one paragraph per mechanism (K103N, the tyrosines, the G190 series, the V106
series, ΔΔETotal against ΔΔGbind, K103N+M230L), then three limitation
paragraphs, then the closing paragraph on mechanisms out of reach.

The three limitation paragraphs are **sampling** (one basin, 100 ns against a
slower apo equilibrium, and 66 of 1,287 comparisons resolved), **the model
system** (no template/primer or dNTP, with the cryo-EM evidence for why a
binary model is still defensible at the drug site), and **the energy functions**
— the new one, written for this revision. Nothing in the paper had interpreted
the force field or the charge model as a limitation even though the
charge-model sensitivity had been measured. It states that the two energetic
methods do not share a ligand charge model, that re-scoring under AM1-BCC moves
ΔΔETotal by up to 1.33 kcal/mol while changing no sign, and that 1.33 exceeds
every gap between adjacent genotypes in the ΔΔETotal ranking (median gap 0.16,
largest 1.29) — so the sign of an interface penalty survives the charge model
but the rank of a genotype does not. It also carries the implicit-solvent point
about G190E's buried glutamate, which used to sit in the mechanism paragraph.

Mechanism paragraphs lost their caveats; the 90-word statistics passage in the
V106 paragraph is now one sentence and the method moved to Methods. The
citations travelled with the sentences as field elements, so each is still
bound to its Zotero item.

**Zotero must be refreshed in Word.** Moved citations are out of numerical
order until then. The bibliography is a `ZOTERO_BIBL` field, so a refresh
renumbers the in-text marks and rebuilds the reference list; nothing was
renumbered by hand.

Applied by `ops/docx/apply_discussion_revision.py`, built on
`ops/docx/docx_surgery.py`, which moves citation fields rather than merely
avoiding them. Both scripts verify paragraph counts, field balance and that no
field was split before writing.

## 5. Figure 3 and Figure 4 captions — in the manuscript

Both captions are written into `Manuscript-Revision-1.docx`, replacing the
one-line stubs, and both reference Supplementary Table 10 for the per-comparison
statistics. Figure 3 is 18 genotypes by 7 metrics with 21 of 126 cells resolved;
Figure 4 is 18 by 5 with 11 of 88 non-blank cells resolved. Both drop the Tyr188
interplanar angle, which no genotype resolves; reinstating it means changing the
dropped-metric sentence with it.

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
`results/analysis/apo_ensembles/tables/pocket_opening_cost.csv`). The holo arm
reproduces the published value to 0.4 Å³ (228.4 against 228), so the recompute
is trustworthy. Wild-type apo is **76.4 ± 9.4 Å³** and opening costs
**152.1 ± 15.0 Å³**, not 137. With errors propagated,
only **G190E (+70.7 ± 26.2 Å³)**, G190S (+65.1 ± 27.0) and V106A+P225H
(−47.9 ± 18.3) differ from wild type. A98G+F227C is +22.5 ± 18.9 and
K103N+M230L is +1.8 ± 28.4, so the claim that all three mismatched genotypes pay
more does not survive, and a large nominal cost is not peculiar to them: G190S
pays the second most and ΔΔGbind tracks its phenotype.

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

**The primer grip loosens, consistently.** Of the 171 genotype-by-region
fluctuation comparisons in the bound state, 17 are resolved at 95%, and they are
not scattered: the primer grip (225–235) becomes more mobile in **nine of the
nineteen genotypes**, by 0.08 to 0.16 Å, always in the same direction. G190E is
more mobile than wild type in four regions at once — pocket lining, primer grip,
β6 and the p66 fingers, the last by 0.56 Å [0.12, 0.99]. Secondary structure
moves by at most 0.11 of strand character and only 4 of 342 comparisons are
resolved. This is the most direct thing the data says about Reviewer 1's
question on the P225 and M230 genotypes, and a binary complex cannot connect it
to a mechanism. It is currently only in the Supplementary Figure 5 caption.

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

- **Refresh Zotero in Word** and check the reference list, as §4 says.
- The G190E replication-capacity specifics the Response used to assert — 20–60%
  of reference replication, compensatory L74V or V75 — are **not** in the
  manuscript. I could not reach Huang 2003 to verify them, so the manuscript
  states only what reference 76's own title supports plus the methodological
  point that the clinical fold change is measured on viruses carrying other
  substitutions. If you can confirm the numbers against the paper, they are
  worth adding.
- The primer-grip finding in §6 is not yet in the main text. It is the most
  direct answer the data gives to Reviewer 1's question about the P225 and M230
  genotypes, and it currently appears only in the Supplementary Figure 5
  caption.
- `mutation_component_summary.csv` is stale and does not reconcile with the
  manuscript. Use `binding_energy/tables/ddg_full.csv`.
- Six `src/nnrti/cli/*.py` and `ops/docx/*.py` files are new this session and
  are committed; nothing is untracked.
- **Author actions:** ORCID iDs, the funding field (NIAID R24AI136618) on the
  submission form, and an English-language pass over the pre-existing text. The
  added words were audited; the original text has not been.

### Fixed this session

- The Response quoted 19 sentences that appeared in neither the manuscript nor
  the Supplementary Text, in four blocks drafted as answers to reviewers and
  never inserted. Methods gained the rest of the unliganded-simulation
  protocol, Results gained the subsection on the cost of opening the pocket,
  and the Supplementary Text gained the sampling discussion. Every quoted
  sentence now appears in the paper: **0 of 90**, checked by re-auditing.
- The figure renumbering was half-done: the mechanism figure was Figure 5 in
  Results but captioned Figure 3, which the new heatmaps also claimed, and the
  Discussion still pointed at Figure 3A/3B. Both heatmap captions were stubs.
- The Supporting Information list described the September submission: one
  Supplementary Text with two figures and four tables, against nine sections,
  six figures and ten tables.
- `apo_crystal_reference.csv` still carried the pre-100 ns wild-type apo volume
  of 91.0 Å³, which fed Supplementary Table 7. It is 76.4. The script that
  writes it also printed a hardcoded "matched 10 ns" label; that is gone.
- Supplementary Figure 5's caption claimed regional RMSF changes reach at most
  0.25 Å for Y181C and G190S. Both the figure and my own earlier note were
  wrong; see §6.
- Supplementary Text section 9 named V106I+F227C as one of the three genotypes
  whose ΔΔGbind misses the phenotype. It is K103N+M230L.

## 9. Where the detail went

`_archive/` holds the eleven superseded memos. Two are still worth opening:

- `R1-crystal-representativeness.md` — the full 4NCG vs 7Z2G evidence and the
  manuscript-ready SI text. The text is already in the documents; the file is
  the record of how the numbers were obtained.
- `deleted-table-notes.md` — verbatim text deleted from the supplementary
  workbook Notes tabs, including the method and limitation lines that have
  nowhere to live except the Supplementary Text.

Technical traps that cost time live in `CLAUDE.md`, not here.
