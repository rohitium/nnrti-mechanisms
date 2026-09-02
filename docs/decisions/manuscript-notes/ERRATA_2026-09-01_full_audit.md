# Full audit of DorDRM-MD-09-02-26.docx + Supplementary Text/Tables

Audited 2026-09-01 against `results/` raw data. Every number in the manuscript was
recomputed from `Supplementary-Table-3.xlsx`, `results/analysis/mechanisms/`,
`results/analysis/fep_pmx/panel_ddg.csv`, `results/analysis/modern_md_suite/`, and
`results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv`.

Severity: **[1]** blocks submission · **[2]** a reviewer will raise it · **[3]** editorial.

---

## A. Data integrity — hard contradictions

### A1. [1] Table 2 ∆∆E_SA column is stale and contradicts every other column
The four other ∆∆E columns match the final protocol exactly. The SA column does not
match the final run, the intermediate `even100` run, or `paper/tables/Table-2-energetics.csv`
(which is correct). Three entries have the **wrong sign**:

| genotype | Table 2 | correct |
|---|---:|---:|
| G190A | −0.016 ± 0.024 | **+0.015 ± 0.008** |
| G190S | +0.024 ± 0.057 | **−0.004 ± 0.028** |
| V106A+L234I | −0.070 ± 0.004 | **+0.003 ± 0.031** |
| L100I+K103N | −0.072 ± 0.010 | −0.017 ± 0.014 |
| G190E | −0.024 ± 0.062 | −0.072 ± 0.063 |
| V106A | +0.057 ± 0.011 | +0.016 ± 0.021 |
| A98G+F227C | +0.016 ± 0.013 | +0.051 ± 0.023 |
| V106M | +0.039 ± 0.056 | +0.054 ± 0.017 |
| V106I | −0.025 ± 0.063 | +0.006 ± 0.056 |
| Y318F | −0.070 ± 0.039 | −0.040 ± 0.023 |
| V106A+P225H | −0.044 ± 0.049 | −0.008 ± 0.018 |

Consequence: **the table fails its own additivity check.** vdW+elec+GB+SA should equal
Total. With the printed SA values it misses by more than rounding for 9 of 18 rows —
V106A+L234I by 0.08, L100I+K103N by 0.05, V106A by 0.05, G190E by 0.05. A referee who
adds the columns will find this in minutes.
**Fix:** paste the ∆∆E_SA column from `paper/tables/Table-2-energetics.csv`.

### A2. [1] Supplementary Table 3 (FEP sheet) still holds the pre-resolution G190E
The sheet's `wt_to_G190E` rows are −0.57 / +4.24 / −0.71, summing to
**∆∆G = +0.99 ± 1.63** — the 20 ns-equilibration campaign. Table 2 in the main text says
**−1.02 ± 0.38** (which matches `results/analysis/fep_pmx/panel_ddg.csv`). The supplement
contradicts the paper, in the opposite sign, for the genotype the Discussion argues about.

### A3. [1] Figure 2D in `Figures.pptx` is the old panel
`image22.png` on slide 8 differs by MD5 from the regenerated
`panel_ddg_vs_experiment_by_category.png`; it still plots G190E at ≈ +2.

### A4. [1] Methods describe a superseded MM/GBSA protocol
Methods: *"Each snapshot was energy-minimized for 100 iterations with all atoms free."*
The reported data used **2000 iterations** (`mmgbsa_relaxation_iterations=2000`,
`snapshot_relaxation=unrestrained`, `frame_sampling=even`, `contact_screened=False`,
double precision). Nobody following the Methods will reproduce Table 2 — the 100-iteration
protocol gives materially different numbers (K103N −0.04 vs −0.75; G190E +3.66 vs +1.94).
Also missing: that **no frames were excluded**, and the precision used.
The replacement text is already drafted in `CHANGES_2026-08-30_final_mmgbsa.md` §7.

### A5. [2] Pyridinone contact numbers for the compound V106A genotypes are not reproducible
The sentence quotes five numbers from **two incompatible computations**:

| | manuscript | `dor_moiety_contacts_summary.csv` | `contact_cutoff_sweep.csv` |
|---|---:|---:|---:|
| WT | 14.7 ± 1.4 | **14.7 ± 1.4** | 14.84 ± 1.66 |
| V106A | 11.7 ± 1.8 | **11.7 ± 1.8** | 11.55 ± 1.83 |
| V106A+F227L | 10.3 ± 0.7 | — | 10.12 ± 0.95 |
| V106A+L234I | 10.5 ± 0.8 | — | 10.73 ± 0.51 |
| V106A+P225H | 11.8 ± 0.1 | — | 11.65 ± 0.19 |

WT and V106A come from the moiety table; the other three match **neither** archived file.
Recompute all five from one source.

### A6. [2] The two sheets of Supplementary Table 3 cover different genotypes
The FEP sheet contains a bare **F227C** row; the MMGBSA sheet does not, and F227C appears
in neither Table 1 nor Table 2. A reader has no way to know it is an intermediate
alchemical leg rather than a panel member. Either add a note or drop it.

### A7. [3] G190E pocket volume
Manuscript: 286 ± 3 Å³. Archived: **284.3 ± 3.2** Å³
(`results/analysis/modern_md_suite/tables/pocket_volume_genotype.csv`). WT 230 ± 12 vs
archived 228.4 ± 11.7. Reconcile or state the source.

---

## B. Logical and scientific problems

### B1. [1] Missing negation destroys the last Discussion sentence
> "…perturbing enzyme activity/fitness, etc., these mechanisms **are captured** in our
> calculations due to cost considerations as well as limitations of classical MD."

Must be **"are not captured"**. As written the paper claims the opposite of what it means,
and the "due to cost considerations" clause is left dangling. (Also a comma splice.)

### B2. [2] That same sentence contradicts the Methods
It lists "opening/closing of NNIBP in the apo configuration" among the mechanisms not
captured — but the FEP protocol explicitly simulates the apo leg, and Methods state the
pocket-volume metric "requires no bound ligand and is therefore directly comparable
between holo and apo simulations." Pick one position.

### B3. [1] The ∆∆E-vs-∆∆G paragraph has the causality backwards
> "…these distortions occur in the apo configuration as well, **leading to increased
> binding affinity** between RT and DOR."

Distortion that occurs in both legs **cancels** — it yields ∆∆G ≈ 0, not ∆∆G < 0. A
resolved negative ∆∆G needs the apo penalty to *exceed* the holo penalty, which is a
different and stronger claim. State it that way or drop the causal language.

Separately, **A98G+F227C does not belong in this sentence with G190E.** G190E is
GB/electrostatics-driven (∆∆E_elec −1.37, ∆∆E_GB +2.50); A98G+F227C is almost purely
van der Waals (+1.81) with negligible elec (+0.30) and GB (+0.31). Same sign, different
mechanism.

### B4. [1] k_B T = 0.5 kcal/mol is wrong
k_B T at 300 K is **0.596 kcal/mol** (0.592 at 298 K). The entire strong/weak
classification — and therefore every "3 of 4" and "7 of 9" count in the abstract, results
and discussion — hangs on this threshold. Either correct the value or stop calling 0.5 the
thermal energy and just call it a 0.5 kcal/mol threshold chosen for convenience.

### B5. [2] The strong/weak rule is one-sided, so improvements are scored as "weak"
The rule reads "more than one standard error **above** 0.5 kcal/mol." G190E
(−1.02 ± 0.38) and A98G+F227C (−1.53 ± 0.87) are resolved *strengthenings* of binding —
2.7σ and 1.8σ below zero — yet are counted as "weak impact," i.e. lumped with genotypes
where nothing happened. This is not a rounding matter: it is the difference between "the
calculation says nothing" and "the calculation says the opposite of the phenotype."

### B6. [2] G190E: "major distortions in the NNIBP" is never reconciled
The Discussion asserts major NNIBP distortion, but (a) ∆∆G = −1.02 ± 0.38 — binding
*improves*, and (b) the Tyr188↔chlorocyanophenyl interplanar angle in G190E is
**12.0 ± 0.9°**, i.e. *more* WT-like than WT itself (13.4 ± 0.4°). The aromatic anchor is
undisturbed. Say explicitly that the distortion is local (Val179 displacement, pocket
expansion) and does not reach the pharmacophore — otherwise the reader is left with a
paragraph that argues against its own free energy.

### B7. [2] The K103N+M230L explanation does not explain K103N+M230L
The preserved backbone H-bond explains why **K103N** costs nothing. It says nothing about
why the *combination*, which is 36-fold resistant, also shows no penalty — M230L is the
driver and is not discussed. As written, "This may explain why K103N+M230L … shows no
resolved binding penalty" is a non sequitur.

### B8. [2] Net-negative interface energies: one genotype omitted
"K103N (−0.75 ± 0.13) and K103N+M230L (−1.42 ± 0.67) had net negative interface energy"
— **L100I+K103N is −0.59 ± 1.01**, also negative, also elec-stabilised (−1.82). Three,
not two. (This was flagged in `CHANGES_2026-08-30_final_mmgbsa.md` §3 and not applied.)

### B9. [2] Abstract Conclusion overclaims, and contradicts the Introduction
Abstract: *"Our findings **demonstrate** that in silico relative binding free-energy
calculations **can correlate** with experimentally observed phenotypic data."*
Introduction: *"neither alchemical FEP nor MM/GBSA calculations demonstrate their utility
as a robust proxy."* The strongest result in the paper is R² = 0.26, p = 0.07 on n = 13 —
which by the paper's own Statistical analysis section is a descriptive effect size, not a
demonstration. Soften "demonstrate" to "suggest," or the two framings will be quoted
against each other.

### B10. [2] "Weak correlation" vs "no correlation" is applied inconsistently
R² = 0.26 (p = 0.07) is a "weak correlation"; R² = 0.13 (p = 0.23) is "no correlation";
R² = 0.18 (p = 0.08) is a "weak correlation." At n = 13–18 none of these is
distinguishable from any other. Use one consistent phrase ("no resolved association") or
report all three as effect sizes with confidence intervals and drop the verbal grading.

### B11. [3] Abstract's ∆∆E_Total statistic is the n = 13 subset, but reads as the full panel
"no correlation observed between ∆∆E_Total and experimental susceptibility (R² = 0.13,
p = 0.23)" — that is the established-phenotype subset. The full-panel value is
R² = 0.06, p = 0.31. Say which set.

### B12. [2] Mutations were introduced only into p66
Methods: "introducing the selected p66 substitutions in silico." In the virus, p51 is
cleaved from the same Gag-Pol precursor, so K103N, Y181C, V106A etc. are present in
**both** subunits. Modelling only the p66 copy is a defensible simplification — the NNIBP
is a p66 pocket — but it is a simplification the paper never states, and the pocket
definition itself reaches into p51 (E138). One sentence of justification is needed.

### B13. [2] parmbsc1 DNA parameters are cited for a system with no DNA
4NCG is the RT–doravirine binary complex; there is no nucleic acid. Citing bsc1 (ref 38)
reads as boilerplate carried over from another system and undermines confidence in the
rest of the Methods.

### B14. [2] Gasteiger charges on the ligand, unacknowledged
Gasteiger charges are a topological approximation; the OpenFF standard for binding free
energies is AM1-BCC. The ligand's electrostatics enter both ∆∆E_elec and the FEP
solvation response directly. This will be the first thing a computational referee asks
about. Either justify it (mutations are protein-side, so ligand charges largely cancel in
∆∆) or acknowledge it as a limitation — but do not leave it unmentioned.

### B15. [2] Replicates are near-degenerate
"Replicate starting structures were generated by applying small coordinate jitter (0.1 Å)."
Three trajectories seeded 0.1 Å apart are not independent samples of the conformational
ensemble; SEM across them is a lower bound on the true uncertainty. Every ± in the paper
inherits this. Worth one honest sentence.

### B16. [2] The V106I/V106M contact "rise" is inside the error bars
15.9 ± 0.6 and 16.7 ± 2.3 against **14.7 ± 1.4** in WT. Neither is resolved. "Consistent
with this crowding, pyridinone contacts rise to…" states as fact something the data do not
separate from zero. Reword to "are unchanged or slightly higher," or supply per-frame
statistics.

Related: the paragraph concludes with "steric crowding and **anchor strain**," but the
anchor-strain evidence — the Lys103 backbone H-bond stretching to 3.62 ± 0.30 Å (V106I)
and 3.56 ± 0.20 Å (V106M) against 2.97 ± 0.01 Å in WT, the longest in the panel — was
dropped from the merged text. Either restore it or delete the phrase.

### B17. [2] The G190E pocket-expansion error bars are asymmetric in a misleading way
"230 ± 12 Å³ to 286 ± 3 Å³." The mutant's ±3 sits next to the WT's ±12; the expansion is
~2.3σ of the WT spread, not the near-certainty the tight mutant error implies. Note also
that **G190S (268 ± 18) and Y181C (257 ± 2)** show comparable expansions and are not
discussed — so pocket expansion is not specific to the charge switch.

### B18. [3] "Flips" is the wrong word, and the G190S angle is bimodal
13.4° → 27.4° is a 14° mean shift, not a flip. And the three G190S replicates are
**34.2 / 17.0 / 30.9°** (median 20.5°): one replicate stays near WT. Describe the
distribution rather than the mean. G190A (16.0 ± 1.8°) is intermediate and unmentioned.

### B19. [3] "in vivo" should be "in vitro"
"…indicating that this strain may be tolerated **in vivo**." The 1.1× and 3.4× fold-changes
are from cell-based phenotypic assays. Also "strain" collides with "viral strain" two words
after two genotype names — use "distortion" or "steric penalty."

### B20. [2] No equilibration stage is described
Heating (10→300 K, 25 ps NVT) is followed straight by production. But MM/GBSA samples "the
final 75 ns," implying 25 ns is discarded as equilibration — which is never stated, and
there is no NPT equilibration step between heating and production. Add both.

### B21. [2] Several headline fold-changes rest on n = 1 from non-peer-reviewed abstracts
93× (A98G+F227C), 105× (V106I+F227C), 36× (K103N+M230L) and 18× (G190E) all come from a
single isolate in ref 56, a conference poster hosted on natap.org; ref 57 likewise. These
are the very genotypes the Discussion treats as the interesting failures. The regression
against log10(fold-change) weights them equally with the six-isolate medians. Acknowledge
the asymmetry in data quality, and consider an isolate-count-weighted sensitivity check.

### B22. [2] No Limitations paragraph and no data/code availability statement
The Introduction advertises "open source alternatives … for transparency and
reproducibility," and the paper then ships no repository link, no accession, and no
availability statement. JCIM requires one.

---

## C. Structure and presentation

### C1. [1] Every mechanistic number appears for the first time in the Discussion
Burial (19.5 / 24.0 / 13.7), interplanar angles (13.4° / 27.4°), Val179 distances,
pocket volumes, Ser105 displacements, per-moiety contacts — none of these appear in
Results, in any table, or in any main-text figure. The Methods section
"Trajectory processing and structural features" describes six analyses whose results are
never reported as results. This is the single most likely structural objection from a
referee. Move them into a Results subsection with a supporting table (or promote
Supplementary Figure 3 to a main figure).

### C2. [2] Supplementary Figure 2 is cited for the wrong content
"The results are summarized in Table 2, Supplementary Figure 2, and Supplementary Table 3"
follows the MM/GBSA description — but Supp Fig 2 is the FEP switching-work distributions.
There is no figure anywhere showing the MM/GBSA component decomposition, despite Table 2
carrying five ∆∆E columns.

### C3. [3] Figure 2 caption is embedded mid-paragraph
It runs on directly from "…provided in Supplementary Table 2." with no break. (Word fix.)

### C4. [3] Table 1's title sits inside the table's header row
Should be a paragraph above the table. (Word fix.)

### C5. [3] Figure 2D caption omits the encoding
It does not mention that points are colour-coded by phenotype category or that the x-axis
is logarithmic — both of which the panel uses.

### C6. [3] Supplementary Figure 1 contradicts the Methods on the RMSD reference
Caption: "RMSD of the DOR heavy atoms vs **the minimized structure** within the NNIBP
(at t = 0)." Methods: "relative to **the crystallographic pose**." Different references.

### C7. [3] Supplementary Figure 3A is cited for a passage it does not cover
The caption covers WT and Y188L only; the citation sits in a passage that also quotes
Y181C (24.0 ± 2.6). Either add Y181C to the panel or move the citation.

### C8. [3] Supplementary Tables have no captions, and Table 4 is orphaned
Supp Tables 1–3 are shipped as bare .xlsx with no legends in the Supplementary Text
document. `Supplementary Table 4.xlsx` exists in the folder and is cited nowhere.

### C9. [2] Abstract is structured and 296 words
JCIM wants an unstructured abstract of ~3–4 sentences (and ≤250 words). Drafts already
exist in `paper/abstract-jcim-unstructured.md`.

### C10. [2] Missing front/back matter
No TOC graphic, ORCID, funding statement, author contributions, conflict of interest, or
Supporting Information listing. All are required.

### C11. [3] γ and λ are Symbol-font `w:sym` glyphs, not Unicode
Four λ (Figure 2 caption) and two γ ("Val106 Cγ") are encoded as
`<w:sym w:font="Symbol">`. These vanish in many PDF/XML/HTML conversion pipelines —
including the plain-text extraction used for this audit. Replace with real Unicode
U+03BB and U+03B3.

---

## D. Table 1 specifics

- **D1. [1]** V106A: "Number of isolates = **2**" but **three** raw values are listed
  ([7.1, 9.6, 28]). The median 9.6 is only correct for n = 3.
- **D2. [2]** Y318F: n = 3 but raw values are "Unavailable (see Table 2 in source study)".
  Either extract them or explain why the median is quotable without them.
- **D3. [2]** The Resistant block is not sorted: 11, 9.6, 93, 105, 106, 149, 153, 161, **36**.
  Y318F precedes a smaller value and K103N+M230L is stranded at the bottom. The other two
  blocks are sorted.
- **D4. [2]** V106A+L234I sources read "Feng et al (2014)**58**; Smith et al (2016)**58**".
  Ref 58 is Smith 2016; Feng 2014 is **ref 16**.
- **D5. [3]** Medians round half-down: 1.25→1.2 (K103N), 7.85→7.8 (K103N+P225H),
  161.5→161 (V106A+L234I). Consistent but non-standard; state the convention or round half-up.

---

## E. Language, style, references

**Prose**
- E1. Abstract: "a weak association **between** … ∆∆Gbind, **and with** log10-transformed
  susceptibility fold-change" — delete "with".
- E2. "…perturbing enzyme activity/fitness, etc., these mechanisms…" — comma splice; also
  avoid "etc." in a formal Discussion.
- E3. "we find **the DOR pose** … has room to accommodate it" — the *pocket* has room, not
  the pose.
- E4. "makes a **hydrophobic bond** with Val106" (Introduction) — no such thing; use
  "hydrophobic contact."
- E5. "**MM/GBSA-based absolute energy computation**" (Discussion) — MM/GBSA here yields a
  WT-referenced *interaction energy*, and the supplement's own Notes sheet says so
  explicitly. Contradicts the paper's own definition.
- E6. "the precise identity of residue 103 is **immaterial** to successful binding" —
  overstated; the paper's own K103N ∆∆G is +0.54 and its ∆∆E_elec shifts by −2.10.
- E7. "the Cys181 cavity is **likely absorbed by a modest local relaxation**" — asserted
  with no supporting measurement.
- E8. Genotype naming is inconsistent across the package: **K103N+L100I** (Table 1,
  Results) / **L100I+K103N** (Discussion, Supp Table 3) / **K103N + L100I** (Table 2,
  with spaces). Pick position-ordered, unspaced, and apply everywhere.
- E9. Results say "3 out of 4 Susceptible **mutations**" / "7 out of 9 Resistant
  **mutations**" for entries that are multi-mutation patterns; the rest of the paper
  carefully says "genotypes."
- E10. Unit spacing: "100ns" (×2) vs "100 ns" (×3); "~10Å"; "final 75ns".
- E11. Mixed minus signs: 7 × U+2212 against 30 × ASCII hyphen, including inside Table 2
  (G190E's −1.02 vs A98G+F227C's -1.53).
- E12. Methods heading capitalization alternates between sentence case and Title Case
  ("Molecular dynamics simulations" vs "MM/GBSA Calculations" vs
  "Free Energy Perturbation Calculations").
- E13. Supp Fig 1: "the standard error **this** mean" → "of this mean".

**References**
- E14. Journal name is malformed as "**Journal American Chemical Society**" in refs 23,
  24, 25, 27 and 70, but correct as "J. Am. Chem. Soc." in ref 28. Zotero style artifact.
- E15. Author names malformed: ref 4 "**AD Clark, J.**" (should be Clark, A. D., Jr. — it
  is correct in ref 9); ref 12 "Clark, **Arthur D.**" and "Lichtenstein, **Mark. A.**"
  (full given names plus a stray period).
- E16. Ref 55 (Asante-Appiah 2021): the page field contains a DOI
  ("65 (12), 10.1128/aac.01216-21").
- E17. Non-breaking hyphens (U+2011) appear inconsistently in refs 62, 64, 67, 73
  ("SARS-CoV‑2", "HIV‑1", "MK‑8591").
- E18. Ref 16 (Feng et al.) is given as **2014**, AAC **59** (1), 590–598 — volume 59
  issue 1 is the January **2015** issue (online-first 2014). Table 1 also cites it as
  "Feng et al (2014)". Verify the year.
- E19. Ref 6 (Cilento, a *Chem. Rev.* review) is used as the primary citation for three
  specific structural claims, including "G190A adds a bulge … causing steric conflict with
  HBY097." Cite the primary structural work alongside it.

---

## Verified correct — no action needed

For the record, these were checked against raw data and are right:

- All 18 ∆∆G_bind values in Table 2 reproduce exactly from the Supp Table 3 FEP sheet by
  summing leg means with SEMs in quadrature.
- All ∆∆E_Total / vdW / elec / GB values reproduce from the MMGBSA sheet.
- All six correlation statistics: R² = 0.079/p = 0.258, 0.064/0.310, 0.182/0.078 (n = 18);
  0.261/0.075, 0.126/0.234 (n = 13). Reported values match to the digit.
- The strong/weak counts (3 of 4, 7 of 9, 1 of 4, 13 above threshold, 9 vdW / 3 GB / 1 elec,
  9 of 18 negative elec) are all arithmetically correct **under the stated 0.5 threshold**.
- Every structural number in the Discussion reproduces from
  `results/analysis/mechanisms/mechanism_coordinates.csv`: burial 19.51/23.98/13.66;
  H-bond 2.97/3.08/3.05/3.10/3.20; Lys103 polar 8.38→5.07; Val179 3.62→6.47;
  Ser105 6.65/5.30/7.15/7.04; residue-106 3.36/3.22/3.16; angles 13.37/27.38.
- Supplementary Figure 2's "19 alchemical legs" is correct (10 single + 9 second legs).
- The 16 pocket-lining residues listed in Methods count correctly (15 p66 + 1 p51).
- Reference callouts are in ascending first-appearance order throughout.
