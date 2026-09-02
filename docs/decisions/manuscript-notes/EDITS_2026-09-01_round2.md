# Round-2 edits — paste-ready

**The .docx is open in Word right now** (`~$rDRM-MD-09-02-26.docx`, lock 18:03; file
saved 18:39). I did not touch it — anything I wrote would be clobbered on your next save.
Everything below is find/replace text. The two supplementary workbooks I *did* regenerate,
because those are mine.

---

# Done by me

| item | status |
|---|---|
| **A1** Table 2 ∆∆E_SA | **verified fixed.** All 18 values now match `Table-2-energetics.csv`, and additivity (vdW+elec+GB+SA = Total) holds on every row, max residual 0.011. |
| **A2** Supp Table 3 | **regenerated.** G190E FEP rows are now −0.57 / −1.78 / −0.71 → **−1.02 ± 0.38**. All 18 genotypes reproduce Table 2 exactly. Written to `paper/` and `paper/submission/`. |
| **A3** Figure 2D | **slide 9 is correct** — `image25.png` is byte-identical to the regenerated `panel_ddg_vs_experiment_by_category.png`. One caution below. |
| **A5 / A7** Supp Table 4 | **built.** New workbook, `paper/submission/Supplementary-Table-4.xlsx` (+ submission copy). |
| **B1** missing negation | **verified fixed** — "which are not captured in our calculations due to cost considerations and limitations of classical MD." The comma splice is gone too. |

### A3 — one caution
Slide 9 carries **both** `image25.png` (correct) and `image22.png` (the old panel, also on
slide 8). Figure 2 in the .docx is a pasted **EMF**, so I cannot read its content to confirm
which one you copied. Quick visual check: in the correct panel **G190E sits below zero**.

### A5 — Supplementary Table 4
`paper/submission/Supplementary-Table-4.xlsx`, three sheets:

- **Summary** — mean ± SEM per genotype, 20 rows (19 + WT + the F227C alchemical intermediate).
- **Per-replicate** — 60 rows, one per genotype × replicate, plus frames analysed.
- **Definitions** — what each column is, the aggregation rule, PBC handling, the ring-vs-moiety
  distinction, and full provenance.

Nine columns, exactly the list you gave: Res103 C=O→triazolinone N; Res103 side-chain polar
atoms→DOR; chlorocyanophenyl ring burial; Tyr188/chlorocyanophenyl interplanar angle; Val179→DOR;
Ser105→DOR; residue 106→DOR; pyridinone moiety contacts; NNIBP volume.

Built by a new committed script, `src/analysis/cli/build_supplementary_table_4.py`, from three
canonical sources in one pass. **Every structural number in the Discussion now traces to one file.**

Two consequences worth knowing:

1. **The moiety analysis had only ever been run on 4 of 20 systems.** I reran it across the whole
   panel, so the pyridinone column is now internally consistent. Five numbers in the Discussion
   shift slightly — see edit **#7**.
2. **G190E pocket volume is 284.3 ± 3.2 Å³, not 286 ± 3** (A7). See edit **#8**.

### Housekeeping
`paper/Supplementary Table 4.xlsx` (with spaces, dated June) is an orphan from an abandoned
linear-regression analysis — unrelated to the new table and cited nowhere. Delete or rename it
before it gets confused with `Supplementary-Table-4.xlsx`.

---

# Edits for you to paste

## 1. A4 — Methods, MM/GBSA minimisation

**Find:**
> Each snapshot was energy-minimized for 100 iterations with all atoms free before energy evaluation, to relieve steric overlaps that would otherwise dominate the van der Waals term. Replicate means are reported from 100 evenly spaced snapshots across the final 75ns of each production trajectory.

**Replace:**
> Each snapshot was energy-minimized with all atoms free before energy evaluation, to relieve steric overlaps that would otherwise dominate the van der Waals term; minimization was run to convergence over 2,000 iterations. Replicate means are reported from 100 snapshots spaced evenly across the post-equilibration portion (final 75 ns) of each production trajectory. No frames were excluded: snapshots containing close interatomic contacts were relieved by this minimization rather than discarded, since excluding configurations on the basis of their energy biases an ensemble average. All energies were evaluated in double precision.

**Why:** the reported data were produced with `relaxation_iterations=2000`, `frame_sampling=even`,
`contact_screened=False`, double precision (verified in
`.checkpoint_mmgbsa_cpu2000_2026-08-29.csv`). The 100-iteration protocol gives materially
different numbers — K103N −0.04 instead of −0.75, G190E +3.66 instead of +1.94.

> ⚠️ **`CHANGES_2026-08-30_final_mmgbsa.md` §7 suggests adding "verified to be sufficient for the
> total binding energy to be stable to within 0.1 kcal/mol out to 5000 iterations." Do not paste
> that clause — no 5,000-iteration run exists in `results/.checkpoints/`.** The only archived runs
> are at 100 and 2,000 iterations. I have left the claim out. Say the word and I will run the
> 5,000-iteration comparison so you can make it truthfully; it is cheap.

---

## 2. B8 — three genotypes have net negative interface energy, not two

**Find:**
> Relative to WT, K103N (∆∆ETotal = -0.75 ± 0.13 kcal/mol) and K103N+M230L (∆∆ETotal = -1.42 ± 0.67 kcal/mol) had net negative interface energy, both driven by significant electrostatic stabilization, ∆∆Eelec < -2 kcal/mol.

**Replace:**
> Relative to WT, the three K103N-containing genotypes K103N (∆∆ETotal = -0.75 ± 0.13 kcal/mol), K103N+M230L (∆∆ETotal = -1.42 ± 0.67 kcal/mol) and K103N+L100I (∆∆ETotal = -0.59 ± 1.01 kcal/mol) had net negative interface energy, all three driven by electrostatic stabilization, ∆∆Eelec < -1.8 kcal/mol.

**Why:** K103N+L100I is −0.59 ± 1.01 with ∆∆Eelec = −1.82. Threshold loosened to −1.8 so it
covers all three. Note this also makes the point cleaner — the effect is specific to K103N.

---

## 3. B9 — Abstract, Conclusion

**Find:**
> Our findings demonstrate that in silico relative binding free-energy calculations can correlate with experimentally observed phenotypic data for canonical DOR susceptible and resistant genotypes.

**Replace:**
> Our findings suggest that in silico relative binding free-energy calculations partially track experimentally observed phenotypic data for canonical DOR susceptible and resistant genotypes, but do not, at this panel size, constitute a robust proxy for in vitro susceptibility testing.

**Why:** removes the contradiction with the Introduction's "neither alchemical FEP nor MM/GBSA
calculations demonstrate their utility as a robust proxy," and stops "demonstrate" from carrying
a p = 0.07 result.

---

## 4. B10 + B11 — Abstract, Results

**Find:**
> Among genotypes with established susceptible or resistant phenotypes, a weak association between relative change in binding free energy, ΔΔGbind, and with log10-transformed susceptibility fold-change was observed (R² = 0.26, p = 0.07). Including the 5 genotypes with uncertain phenotype, however, this correlation was negligible (R² = 0.08, p = 0.26). In contrast, while MM/GBSA also correctly identified unfavorable energetic changes for 7 of 9 resistant genotypes, only 1 out of 4 susceptible mutations was identified as low-impact and there was no correlation observed between ΔΔETotal and experimental susceptibility (R² = 0.13, p = 0.23). ΔΔGbind and ΔΔETotal were themselves modestly correlated (R² = 0.18, p = 0.08), with van der Waals interactions frequently dominating MM/GBSA energetic shifts.

**Replace:**
> Across the 13 genotypes with established susceptible or resistant phenotypes, ΔΔGbind accounted for 26% of the variance in log10-transformed susceptibility fold-change (R² = 0.26, p = 0.07); across the full 18-genotype panel this fell to R² = 0.08, p = 0.26. MM/GBSA likewise identified unfavorable energetic changes for 7 of 9 resistant genotypes, but only 1 of 4 susceptible mutations as low-impact, and ΔΔETotal accounted for less of the variance on both the same 13-genotype subset (R² = 0.13, p = 0.23) and the full panel (R² = 0.06, p = 0.31). ΔΔGbind and ΔΔETotal were themselves only loosely related (R² = 0.18, p = 0.08), with van der Waals interactions frequently dominating MM/GBSA energetic shifts. At this panel size none of these associations is statistically resolved.

**This fixes three things at once:** the "between … and with" grammar; the verbal grading
(*weak / negligible / no / modest* → variance explained, plus one honest closing sentence);
and B11 — the ΔΔETotal R² = 0.13 is now explicitly labelled as the 13-genotype subset, with
the full-panel R² = 0.06 given alongside.

---

## 5. B10 — Results, correlation paragraph

**Find:**
> While ∆∆Gbind and ∆∆ETotal exhibit a weak correlation with each other (R2 = 0.18, p = 0.08), both of these metrics show no correlation with DOR fold-change values observed experimentally - ∆∆Gbind vs log10(Fold-change) (R2 = 0.08, p = 0.26) and ∆∆ETotal vs log10(Fold-change) (R2 = 0.06, p = 0.31). However, amongst the set of the 13 established DOR Susceptible and DOR Resistant phenotype mutations (i.e. dropping the Uncertain phenotype set), ∆∆Gbind exhibits weak correlation with log10(Fold-change) (R2 = 0.26, p = 0.07), while ∆∆ETotal still shows no correlation with log10(Fold-change) (R2 = 0.13, p = 0.23).

**Replace:**
> ∆∆Gbind and ∆∆ETotal track each other only loosely (R2 = 0.18, p = 0.08). Across the full panel, neither tracks the experimentally observed DOR fold-change: ∆∆Gbind vs log10(Fold-change) gives R2 = 0.08, p = 0.26, and ∆∆ETotal vs log10(Fold-change) gives R2 = 0.06, p = 0.31. Restricting to the 13 genotypes with established DOR Susceptible or DOR Resistant phenotypes (i.e. dropping the Uncertain phenotype set) raises the ∆∆Gbind association to R2 = 0.26, p = 0.07 and the ∆∆ETotal association to R2 = 0.13, p = 0.23. At this panel size none of these associations is statistically resolved, and consistent with the Statistical analysis section we report them as descriptive effect sizes rather than as evidence for or against a correlation.

---

## 6. B13 — remove the DNA force field

**Find:**
> Protein and nucleic-acid parameters were assigned from the Amber14 OpenMM force-field files, using ff14SB for protein,37 bsc1 for DNA,38 and TIP3P water.39

**Replace:**
> Protein parameters were assigned from the Amber14 OpenMM force-field files, using ff14SB37 with TIP3P water.39

**Why:** 4NCG is the RT–DOR binary complex; there is no nucleic acid in any system simulated here.

> ⚠️ **This orphans reference 38 (Ivani et al., parmbsc1).** If your bibliography is a
> Zotero/Mendeley field, deleting the citation renumbers 39–79 automatically. If the reference
> list is hand-numbered, every subsequent callout shifts by one — check before you delete.
>
> Do not simply leave ref 38 in the list uncited; that is worse than the original problem.

**And the absence of the template/primer is not just a tidy-up — it does real work below (B7).**

---

## 7. A5 — updated pyridinone contacts (V106A series)

Now that the moiety analysis has been run across all 20 systems from one source, five numbers move
slightly. All are in Supplementary Table 4.

**Find:**
> we observe a drop in the number of pyridinone moiety contacts within 4.0 Å from 14.7 ± 1.4 in WT to 11.7 ± 1.8 in V106A, 10.3 ± 0.7 in V106A+F227L, 10.5 ± 0.8 in V106A+L234I, and 11.8 ± 0.1 in V106A+P225H.

**Replace:**
> we observe a drop in the number of pyridinone moiety contacts within 4.0 Å from 14.8 ± 1.4 in WT to 11.6 ± 1.7 in V106A, 10.3 ± 0.6 in V106A+F227L, 10.3 ± 0.8 in V106A+L234I, and 11.8 ± 0.1 in V106A+P225H (Supplementary Table 4).

---

## 8. A7 — G190E pocket volume

**Find:**
> and the NNIBP volume expands from 230 ± 12 Å³ to 286 ± 3 Å³.

**Replace:**
> and the NNIBP volume expands from 228 ± 12 Å³ to 284 ± 3 Å³ (Supplementary Table 4).

---

## 9. B16 — V106I/V106M paragraph

Two problems: the "rise" is inside the error bars, and "anchor strain" is asserted with the
evidence missing. Both fixed, and the anchor numbers restored.

**Find:**
> Consistent with this crowding, pyridinone contacts rise to 15.9 ± 0.6 in V106I and 16.7 ± 2.3 in V106M against 14.7 ± 1.4 in WT. Both genotypes nonetheless remain phenotypically susceptible or near-susceptible (1.1- and 3.4-fold), indicating that this strain may be tolerated in vivo, and they account for the two largest overestimates of ∆∆Gbind in our panel (+2.27 ± 0.74 and +6.10 ± 0.16 kcal/mol).

**Replace:**
> Pyridinone contacts are correspondingly unchanged or slightly higher, 16.0 ± 0.8 in V106I and 16.8 ± 2.4 in V106M against 14.8 ± 1.4 in WT, a shift that is not resolved against replicate variation. The clearer signal is at the Lys103 anchor: the residue 103 main-chain carbonyl to triazolinone nitrogen distance, 2.97 ± 0.01 Å in WT and within 0.25 Å of that value in every other genotype we simulated, lengthens to 3.62 ± 0.30 Å in V106I and 3.56 ± 0.20 Å in V106M — the two longest in the panel (Supplementary Table 4). Both genotypes nonetheless remain phenotypically susceptible or near-susceptible (1.1- and 3.4-fold in vitro), indicating that this steric penalty is tolerated by the virus, and they account for the two largest overestimates of ∆∆Gbind in our panel (+2.27 ± 0.74 and +6.10 ± 0.16 kcal/mol).

**Also fixes B19** ("in vivo" → "in vitro"; "this strain" → "this steric penalty"). The closing
sentence's "steric crowding and anchor strain when it is enlarged" is now earned.

*Verified:* 3.62 ± 0.30 and 3.56 ± 0.20 are the two largest values in the panel; every other
genotype lies between 2.80 and 3.22 Å.

---

# B14 — Gasteiger charges: what the literature actually says

You asked me to check whether this is justifiable. **Honest answer: not on charge-quality grounds.**
The one systematic benchmark that tested exactly this — Gasteiger against RESP, ESP and AM1-BCC
for MM/PBSA and MM/GBSA, 46 ligands across 5 receptors — ranks Gasteiger last:

> Xu, L.; Sun, H.; Li, Y.; Wang, J.; Hou, T. *Assessing the Performance of MM/PBSA and MM/GBSA
> Methods. 3. The Impact of Force Fields and Ligand Charge Models.* **J. Phys. Chem. B 2013**,
> *117* (28), 8408–8421. DOI: 10.1021/jp404160y

RESP performs best; AM1-BCC and ESP are "fairly satisfactory"; Gasteiger is not among the
recommended models. Citing that paper would argue *against* you, so don't.

**What is defensible is a cancellation argument, and it is genuinely strong for most of this panel** —
but state it as a limitation, not as a justification:

- The ligand is chemically identical in all 19 systems. Every DOR self-term — intramolecular
  energy, ligand internal GB, ligand SA — cancels *exactly* in the WT-referenced ∆∆E and in the
  ∆∆G thermodynamic cycle. This is arithmetic, not an approximation.
- The residual sensitivity is confined to the ligand's **electrostatic coupling to the mutated
  residue**. For most of the panel that coupling is small, because the substitutions are
  hydrophobic↔hydrophobic (V106A/I/M, Y188L, L234I, F227L, L100I, A98G).
- It is **not** small for the polar and charge-changing legs — K103N, G190E, G190S, Y318F, Y181C —
  and K103N and G190E are precisely the two legs your Methods already single out for special
  treatment. So the argument does not fully cover the cases that matter most.

The protein-mutation RBFE precedent to cite for the method itself:

> Aldeghi, M.; Gapsys, V.; de Groot, B. L. *Accurate Estimation of Ligand Binding Affinity Changes
> upon Protein Mutation.* **ACS Cent. Sci. 2018**, *4* (12), 1708–1718. DOI: 10.1021/acscentsci.8b00717

**My recommendation.** Rather than argue, run the sensitivity check — it is cheap and it converts
a referee's objection into a result. Re-evaluate MM/GBSA on the *same* 100 snapshots with AM1-BCC
ligand charges and report whether the panel-level conclusions move. No new MD is needed for the
energy re-evaluation, and one sentence ("panel conclusions were unchanged when the MM/GBSA
electrostatics were re-evaluated with AM1-BCC ligand charges; see Supplementary Table X") closes
the question. It cannot rescue the FEP ensembles, which were generated with Gasteiger charges, so
the limitation sentence is still needed — but it would be a bounded, quantified limitation instead
of an open one.

Say the word and I will run it.

---

# B7 — K103N+M230L: I think there is a real answer, and it is in your own data

You said you weren't sure how to answer this. The literature gives a clean one, and it happens to
depend on the very simplification that edit **#6** exposes.

**Met230 is a primer-grip residue.** It sits in the β12–β13 hairpin of p66 (residues 227–235), and
in RT/DNA crystal structures its side chain contacts the ribose of the primer nucleotide at
position −2. The primer grip holds the primer terminus in the orientation required for nucleophilic
attack on the incoming dNTP.

> Ghosh, M.; Jacques, P. S.; Rodgers, D. W.; Ottman, M.; Darlix, J.-L.; Le Grice, S. F. J.
> *Alterations to the Primer Grip of p66 HIV-1 Reverse Transcriptase and Their Consequences for
> Template-Primer Utilization.* **Biochemistry 1996**, *35* (26), 8553–8562. DOI: 10.1021/bi952773j

**And M230L is measurably an enzyme-function mutation, not only an inhibitor-contact mutation.**
Xu et al. showed M230L impairs RT enzymatic function and costs the virus roughly 8-fold in
replication capacity:

> Xu, H.-T.; Quan, Y.; Schader, S. M.; Oliveira, M.; Bar-Magen, T.; Wainberg, M. A. *The M230L
> Nonnucleoside Reverse Transcriptase Inhibitor Resistance Mutation in HIV-1 Reverse Transcriptase
> Impairs Enzymatic Function and Viral Replicative Capacity.* **Antimicrob. Agents Chemother. 2010**,
> *54* (6), 2401–2408. DOI: 10.1128/aac.01795-09

**Our systems contain no template/primer.** 4NCG is a binary RT–DOR complex. A mechanism that acts
through the primer-grip/template-primer interface is structurally absent from every simulation in
this study — so a null result for K103N+M230L is what the model *must* produce, and its 36-fold
phenotype is not evidence that the calculation failed.

Supplementary Table 4 supports this positively rather than by hand-waving: in K103N+M230L the DOR
interface is if anything *tighter* than WT — chlorocyanophenyl burial 22.8 ± 2.4 against 19.5 ± 0.5,
pyridinone contacts 16.6 ± 0.6 against 14.8 ± 1.4, the Lys103 backbone H-bond intact at
3.05 ± 0.03 Å. Nothing at the binding interface is disturbed. That is a clean, defensible finding,
and it is exactly the paper's thesis: affinity and susceptibility are not the same quantity.

**Suggested replacement for the current non sequitur** (which claims the preserved H-bond explains
the combination, when it only explains K103N):

> **Find:** This may explain why K103N+M230L, despite a 36-fold reduction in susceptibility, shows no resolved binding penalty in our simulations (∆∆Gbind = 0.59 ± 0.41 kcal/mol).
>
> **Replace:** This accounts for the absence of a binding penalty at position 103, but not for the 36-fold reduction in susceptibility carried by K103N+M230L, for which our simulations likewise resolve no penalty (∆∆Gbind = 0.59 ± 0.41 kcal/mol) and in which the DOR interface is if anything tighter than WT (chlorocyanophenyl burial 22.8 ± 2.4 against 19.5 ± 0.5 heavy atom contacts; Supplementary Table 4). Met230 lies in the p66 primer grip, where its side chain contacts the primer at position −2,<ref Ghosh 1996> and M230L impairs RT enzymatic function and viral replicative capacity in addition to reducing NNRTI susceptibility.<ref Xu 2010> Our binary RT–DOR systems contain no template/primer, so a resistance mechanism acting through that interface is structurally absent from the model — an instance of the more general point that inhibitor affinity and viral susceptibility are distinct quantities.

---

# B6 — Results subsection

You are right that the Discussion cannot carry this until the numbers appear as results. Now that
Supplementary Table 4 exists, the subsection writes itself from one table. **Say go and I will draft
it in full** — my proposal for scope:

**Heading:** *Structural consequences of DOR-associated substitutions* — third Results subsection,
after "RT-DOR binding free energy analyses."

**Content**, all from Supplementary Table 4, roughly four paragraphs:

1. *The Lys103 anchor is preserved across the panel.* 2.97 ± 0.01 Å in WT; 2.80–3.22 Å in every
   genotype except V106I (3.62 ± 0.30) and V106M (3.56 ± 0.20). One sentence on K103N shortening
   the side chain toward DOR (8.38 → 5.07 Å).
2. *Aromatic packing is set by Tyr188.* Burial 19.5 ± 0.5 in WT, 13.7 ± 1.5 in Y188L — the only
   genotype with a resolved loss. Y181C is 24.0 ± 2.5, i.e. not reduced.
3. *Position 106 supports two opposite responses.* The V106A series slides toward Ser105
   (6.65 → 5.13–5.30 Å) and loses pyridinone contacts (14.8 → 10.3–11.8); V106I/M move away
   (7.04–7.15 Å) and strain the Lys103 anchor.
4. *Position 190 and pocket volume.* Interplanar angle 13.4 ± 0.4° WT, 16.0 ± 1.8° G190A,
   27.4 ± 5.3° G190S, 12.0 ± 0.9° G190E; NNIBP volume 228 ± 12 WT vs 284 ± 3 G190E, 268 ± 18 G190S,
   257 ± 2 Y181C.

**This is also what unblocks B6 in the Discussion.** Paragraph 4 is the honest resolution: G190E's
distortion is *local* — Val179 displaced 3.62 → 6.47 Å, pocket expanded — while the pharmacophore
is untouched, its Tyr188 stack (12.0 ± 0.9°) being if anything more WT-like than WT (13.4 ± 0.4°).
That is why a badly deformed pocket still binds DOR at ∆∆G = −1.02. Once the numbers are in
Results, that becomes a two-sentence addition to the Discussion instead of an assertion.

Two things worth saying plainly in that paragraph, since Supplementary Table 4 now makes them
visible to any reader: **G190S (268 ± 18) and Y181C (257 ± 2) expand the pocket comparably**, so
expansion is not specific to the charge switch; and **WT's own volume SEM is ±12**, so the G190E
expansion is about 2.3σ of the WT spread, not the near-certainty the mutant's ±3 suggests.

---

# Still open from the audit, not addressed this round

Table 1 (**D1** V106A n = 2 vs three raw values; **D2** Y318F; **D3** unsorted Resistant block;
**D4** Feng superscripted 58 instead of 16); **B12/B15/B4/B5/B3** — you've ruled these out;
**B17/B18/B21/B22**; **C1–C11** (front matter, abstract length/structure, Word layout fixes,
Symbol-font λ/γ); **E1–E19** except E1 and E19, which are folded into edits #4 and #9 above.
