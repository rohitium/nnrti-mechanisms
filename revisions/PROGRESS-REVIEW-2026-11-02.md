# Where the revision stands, and what is stale

Written 02-Nov-2026, ~03:00, from the three `.docx` files as they stood at
02:40. Backups of all four documents are in
`revisions/_backups/20261102-0246/`. I did **not** edit any of them — see the
last section for why, and for what is ready to paste.

## 1. Progress

The Response document is in good shape structurally. Every reviewer comment
has a Response block, most of it marked red as required. Reviewer 1 is
answered in full, including the two new SI sections (7Z2G crystal
representativeness, the perturbed-pose ladder) that you wanted to lead with.
Reviewer 3 is essentially complete. Figures 1, 3 and 4 are referenced and the
Table 3 → Figure 3/4 replacement is recorded.

What is not done: the Discussion reorganisation, and the stale numbers below.

## 2. Stale numbers — the important part

**Every apo pocket-volume number quoted in the Response predates the 100 ns
apo runs.** They come from the old matched 10 ns window. The 100 ns values:

| genotype | quoted (10 ns) | current (100 ns) | |
|---|---|---|---|
| G190E | −49.0 | **−33.1** | still collapsed |
| V106I+F227C | −36.6 | **−17.7** | still collapsed, half as much |
| A98G+F227C | −28.4 | **−16.3** | still collapsed, half as much |
| V106A | +12.5 | **+54.8** | **4x larger — changes the story** |
| K103N+M230L | +2.0 | **+9.7** | no longer "unchanged" |

Two consequences for the quoted Discussion paragraph under Reviewer 2,
comment 2:

- The sentence "two of the three genotypes whose ΔΔGbind misses the phenotype
  therefore have markedly more collapsed apo pockets than WT" survives for
  G190E and A98G+F227C, but K103N+M230L is now **+9.7**, i.e. more open than
  wild type, not unchanged.
- The whole caveat beginning "Among the genotypes whose apo trajectories run
  well beyond 10 ns…" and ending "The matched 10 ns window is therefore a
  limitation of this comparison and not a controlled one" is **moot**. Every
  apo replicate is now 100 ns. That paragraph should go, not be updated.

The pocket-opening costs quoted under Reviewer 1, comment 3 are stale the same
way: WT 137 → **147.5**, G190E +105 → **+99.2**, A98G+F227C +38 → **+20.9**,
K103N+M230L +20 → **+8.5**. The last one matters: at +8.5 it is identical to
F227C's +8.5 and below ten other genotypes, so pocket-opening cost no longer
distinguishes K103N+M230L from genotypes ΔΔGbind predicts correctly. Dropping
it from that argument is the honest move and it agrees with Reviewer 2's
point that the primer-grip route is foreclosed anyway.

One reference looks wrong: the Reviewer 2 comment 2 text cites the 4NCG/7Z2G
agreement as "Supplementary Figure 5", but the Reviewer 1 response introduces
it as Supplementary Figure 1.

## 3. Effect sizes: one method for every feature

You asked for something simple and general built on the triplicate SEMs.
`nnrti.cli.compute_structural_effect_sizes` now does this for every genotype
and every structural feature in both states — 1,287 comparisons —
`results/analysis/structural_metrics/tables/structural_effect_sizes.csv`:

    shift      difference of replicate means, in the feature's own units
    sem        standard error of that difference
    t, df, p   Welch, unequal variance
    ci95       shift +/- t(0.975, df) x sem
    hedges_g   standardised effect size, small-sample corrected

**A threshold correction you need to know about.** Twice the standard error is
not a 95% criterion at n = 3. For a difference of two three-replicate means
Welch gives df ≈ 4, so 2x SEM is p = 0.12; the 95% interval needs **2.78x
SEM**. The figures used 2x until tonight and now use the right value per
panel — 2.776 for Figure 3 (mutant vs WT) and 2.306 for Figure 4, which is a
difference of differences over four groups and so has df ≈ 8.

**What this does to the claims.** Across all 1,287 comparisons, 66 reach 95% —
which is the null expectation for that many tests. That is worth stating
plainly about the panel as a whole: most structural features do not differ
between genotypes at three replicates. The pre-specified comparisons are a
different matter and are strong:

- V106A toward Ser105: **−1.35 Å [−1.73, −0.97]**, p = 0.004, g = −9.5.
  This reproduces the interval already in your R2.3 response exactly.
- K103N side-chain to DOR: −3.31 Å [−3.89, −2.73], p = 0.001, g = −14.2.
- G190E clears 95% on **eight** holo features, including Val179 distance
  (+2.84 Å [1.01, 4.68]) and NNIBP volume.
- V106M: **zero** features reach 95%. The existing sentence "the V106A versus
  V106M contrast is within error on every metric we examined" is correct.

## 4. Reviewer 2, comments 2 and 3

**R2.2** — the response text is fine; the Discussion paragraph it quotes is
not, for the reasons in §2. It also quotes the same paragraph that R1.3 quotes,
so the two responses currently duplicate each other. Suggest: R1.3 keeps the
mechanism argument, R2.2 keeps only the foreclosed-mechanisms point plus the
subdomain numbers in [86], which are unaffected.

**R2.3** — stronger than you think, and only needs two edits. The numbers
reproduce. But the opening sentence says "bootstrap confidence intervals"
while the quoted text says "Welch intervals on the difference of replicate
means"; they are not the same method and the quoted text is the one that
matches the computation. And the response answers only the two comparisons
the reviewer named, where you can now point to the full table.

Draft replacement for the R2.3 opening paragraph:

> Response: We accept that several mechanistic statements rested on means with
> replicate standard errors. Rather than treat the comparisons the reviewer
> names as special cases, we now report an effect size and a 95% confidence
> interval for every genotype and every structural feature, in both the bound
> and unliganded states (Supplementary Table X). The replicate, not the frame,
> is the unit of replication. The intervals changed what we claim.

Then keep [91]/[92] as they stand, and add after them:

> Applied across the panel, this criterion is also what decides which features
> appear in Figures 3 and 4: a feature whose shift from wild type does not
> exceed its 95% interval for any genotype is not plotted. Five such features
> were removed, among them the subdomain and YMDD-motif fluctuations, which we
> report as unresolved at three replicates rather than as unchanged.

## 5. Discussion

**Not started.** You said twice that my prose aims at the wrong things, so I
have deliberately not written any of it, including tonight. The reorganisation
you described — three opening paragraphs, then one paragraph per mechanism,
then one paragraph per limitation (sampling, simplified model system, missed
mechanisms) — is yours to write.

What I can do without writing it: the stale numbers in §2 all live in
Discussion text quoted by the responses, so that text has to change regardless
of how you reorganise. If you tell me which paragraphs survive the
reorganisation, I will check every number in them against the current tables.

## 6. Why I did not edit the .docx files

You said to go ahead, and the Response document is safe to edit — I checked,
and none of the paragraphs in question contain Zotero fields. The manuscript is
a different matter: it holds **97 live Zotero fields**, and a programmatic edit
that lands near one can break the citation.

The larger reason is that every change above is a judgement about what the
paper claims, not a typo fix, and you have twice told me that my placement
decisions make more work for you. The text in §4 is ready to paste; nothing in
it depends on my choosing where it goes.
