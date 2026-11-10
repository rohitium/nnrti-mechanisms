# Coverage audit — every reviewer point, checked against what exists

Audited 21-Oct-2026 against the three reports and both editorial emails.
Status is judged on what is actually in the repository, not on what was planned.

**✅ addressed · ⚠️ partially addressed, limitation stated · ❌ not done**

---

## Reviewer 1

| # | Request | Status | Evidence |
|---|---|---|---|
| General | Sampled pose treated as representative without testing | ✅ | Discussion paragraph on basin-limited sampling; apo pocket measurements; mutation-distance panel |
| General | Is the 4NCG pose representative of RT–DOR at all? | ✅ | **Added 25-Oct-2026.** 4NCG compared with 7Z2G (3.65 Å cryo-EM, RT + DOR + 38-mer DNA, independent group): DOR 0.62 Å heavy-atom RMSD in the NNIBP frame, pocket 0.49 Å Cα, residue-to-drug distances agree to 0.16 Å median with no superposition. New SI section 8, Supplementary Figure 5, Supplementary Table 8. `results/analysis/crystal_pose_comparison/`, `revisions/R1-crystal-representativeness.md` |
| 1 | Weak correlation may reflect mechanism, not only FEP error | ✅ | Confusion matrix + AUC in Methods/Abstract; Discussion attributes the three misses to pre-binding mechanisms, with apo and distance evidence |
| 2a | Test representative systems from independent starting poses | ❌ | **Not done.** Replicates differ only by 0.1 Å jitter. Stated plainly in Discussion and in the response rather than glossed |
| 2b | Longer apo simulations for P225/M230 genotypes vs WT | ⚠️ | Matched 10 ns apo comparison across 8 systems incl. K103N+M230L and V106A+P225H; drift test shows 1–5 Å³ change beyond 10 ns. **No new simulations**, and **K103N+P225H has no apo trajectory at all** |
| 3 | Discuss mismatches mechanistically | ✅ | All three (K103N+M230L, A98G+F227C, G190E) have larger induced-fit gaps than WT; M230L and A98G are the only positions >5 Å from DOR |
| 4 | Figure mapping mutations onto structure | ✅ | `mutation_map.png/pdf` |
| 4 | WT-referenced heatmap of Table 3 | ✅ | `structural_heatmap_wt_referenced.png/pdf` |
| 4 | RMSF | ✅ | `ca_rmsf_per_residue.csv`, `rmsf_delta_by_region.csv` |
| 4 | Secondary-structure changes | ✅ | `secondary_structure_per_residue.csv` |
| 4 | Pocket conformational distributions | ✅ | `pocket_volume_distributions.png/pdf` |

## Reviewer 2

| # | Request | Status | Evidence |
|---|---|---|---|
| 1 | Block averaging | ✅ | Pre-existing `block_averages.csv`, now reported |
| 1 | Autocorrelation analysis | ✅ | `autocorrelation.csv`, τ_int 1.0–3.7 ns, ESS 14–50 |
| 1 | …for the FEP/MM-GBSA estimates too | ✅ | Argued by construction: MM/GBSA snapshots are evenly spaced, FEP work values come from independently seeded switches, so neither is a serial time series |
| 1 | Discuss how limited sampling may bias results | ✅ | Discussion; within-basin caveat stated explicitly |
| 2 | Which mechanisms are inaccessible without nucleic acid/dNTP | ✅ | Discussion names primer grip (M230L) and the ternary-complex route. Strengthened 25-Oct-2026: the 7Z2G comparison shows the NNIBP and the DOR pose are the same with DNA bound, while p51, the fingers and the thumb move by 1.6–2.2 Å — so what nucleic acid changes is outside the drug site |
| 3 | Effect sizes with CIs for V106A vs V106I/M, G190A vs G190E | ✅ | `effect_sizes.csv`; 7 of 42 intervals exclude zero; **V106A vs V106M survives on no metric** and the Discussion is tempered accordingly |
| 3 | Per-atom or torsion-level analysis | ✅ | `dor_torsions_summary.csv`, four dihedrals, circular statistics |
| 4 | Protonation states consistent WT vs mutant | ✅ | Verified: 13 histidines, all HID, identical across all 19 systems |
| 4 | Hybrid residue topologies | ✅ | Methods: pmx single-topology construction described |
| 4 | Soft-core parameters | ✅ | Methods: α = 0.3, σ = 0.25 nm, power 1, sc-coul = yes |
| 4 | Finite-size correction | ✅ | Leading term computed; omitted terms bounded at ≤0.003 kcal/mol |
| 5 | Compare a more accurate charge model, or discuss | ✅ (both) | Charge sets compared (RMS 0.165 e, mean \|q\| 0.228 vs 0.138 e) **and** MM/GBSA recomputed with AM1-BCC: ΔΔE_Total shifts ≤0.81 kcal/mol, signs preserved, shift confined to charge-dependent components |

## Reviewer 3

| # | Request | Status | Evidence |
|---|---|---|---|
| 1 | Define "correctly predicted" quantitatively | ✅ | Abstract states the +0.5 kcal/mol cutoff; Methods state the rule as implemented, (estimate − SEM) > 0.5, with the full confusion matrix and AUC |
| 2 | Typographical issues | ✅ | The three cited strings are absent from our source files; response says the manuscript was proofread in full and explains this honestly |
| 3a | Justify full minimization | ✅ | Methods, with four citations (Kuhn 2005, Rastelli 2010, Weis 2006, Genheden & Ryde 2015) |
| 3b | Report sensitivity of ΔΔE_Total and components | ✅ | Supplementary Table 5; sensitivity is large and the manuscript now says so |

## Editorial

| Item | Status | Note |
|---|---|---|
| TOC graphic (.tif/.pdf, last page, labelled, no caption) | ✅ | `TOC-Graphic.pdf/.tif`, 600 dpi, 3.3 × 1.8 in, and placed on the last page of both manuscript copies, labelled "TOC Graphic" |
| SI cover page | ✅ | `Supporting-Information-Cover-Page.docx` |
| Clean unmarked manuscript | ✅ | `Manuscript-Revision-1-CLEAN.docx`, text verified identical to the marked copy |
| Marked manuscript for review | ✅ | `Manuscript-Revision-1.docx`, all additions in red |
| Clean SI for publication | ✅ | `Supplementary-Text-Revision-1-CLEAN.docx` |
| Point-by-point response | ✅ | `Response-to-Reviewers.docx` |
| Funding (NIAID R24AI136618) | ⚠️ | In the manuscript; **confirm on the submission form** |
| ORCID iDs validated | ❌ | **Each author must do their own** |
| English-language editing | ⚠️ | The ~1,100 words added in this revision were audited and two faults fixed. A pass over the **existing** manuscript text is an author decision; consider ACS Authoring Services |
| Cover art (optional) | ❌ | Not attempted |

---

## What remains, in priority order

1. **ORCID** and the **funding form field** — author actions, cannot be automated.
2. **English-language pass over the pre-existing text** — editor-requested; my additions
   have been audited but the rest of the manuscript has not.
3. **R1.2a (independent starting poses)** — declined. The single remaining
   scientific request not acted on. If a second round pushes, this is what
   they will push on. The 7Z2G comparison does **not** cover it: two
   experimental structures agreeing shows the pose is right, not that the
   simulations explore more than one basin around it.
4. **Fold the 7Z2G work into the documents** — the analysis, figure and draft
   text exist; SI section 8, Supplementary Figure 5, Supplementary Table 8 and
   the Discussion sentence are not yet in the .docx files. Insertion points and
   final text: `revisions/HANDOFF-7Z2G-docx-insertions.md`. Needs an author
   decision on reference numbering (the list is 88 manually numbered entries).
5. **K103N+P225H apo** — no trajectory exists, so this P225 genotype is absent
   from the apo comparison. V106A+P225H covers the P225 class.
