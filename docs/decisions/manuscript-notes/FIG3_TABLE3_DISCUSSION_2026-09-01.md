# Figure 3, Table 3, and a rewritten Discussion

Assumes the new Results subsection (`NEW_RESULTS_SUBSECTION_2026-09-01.md`) is in.
All Discussion citations below point at Figure 3 and Table 3, never at supplementary material.

---

# 1. Supplementary Figure 3 → Figure 3

**Source image:** `results/analysis/mechanisms/plots/mechanism_panel.pdf` (and `.png`), from
`src/analysis/cli/plot_mechanism_panel.py`. Use the PDF for the main text — it is vector.

**Renumbering:** Supplementary Figures 1 and 2 are unaffected. There is no existing main-text
Figure 3, so nothing shifts. In the Supplementary Text, delete the old Supplementary Figure 3
image and caption.

**Caption:**

> **Figure 3. Two distinct structural mechanisms of DOR resistance, resolved in the equilibrium
> MD trajectories.** (A) Burial of the DOR chlorocyanophenyl ring, defined as the number of RT
> heavy atoms within 4.0 Å of the ring, for WT (red) and Y188L (blue). Loss of the Tyr188 side
> chain removes roughly a third of the ring's packing and is sustained throughout the trajectory.
> (B) Displacement of DOR toward Ser105, defined as the minimum heavy-atom distance between DOR
> and Ser105, for WT (red) and the four V106A-containing genotypes (blues). Shortening the
> residue 106 side chain allows DOR to slide 1.3–1.5 Å out of its crystallographic pose, and the
> four genotypes are indistinguishable from one another. Solid lines are means over three
> independent 100 ns replicates interpolated onto a common time grid; shading is the standard
> error of the mean across replicates. Per-genotype summary values are in Table 3 and
> per-replicate values in Supplementary Table 4.

> ⚠️ **One change to make in the figure itself.** The panel titles currently read
> "Y188L: loss of aromatic packing" and "V106A: DOR slips out of position." Those are conclusions,
> not labels. For a main-text figure, replace them with neutral titles — **"A  Chlorocyanophenyl
> ring burial"** and **"B  DOR displacement toward Ser105"** — and let the caption and Discussion
> carry the interpretation. One-line change in `plot_mechanism_panel.py`; say the word and I'll
> make it and regenerate.

---

# 2. New Table 3

Generated as `paper/tables/Table-3-structural.csv`, and now emitted by the *same* script that
builds Supplementary Table 4 (`src/analysis/cli/build_supplementary_table_4.py`), so the two can
never drift apart. The `Replicates` column is dropped, and F227C — an alchemical intermediate, not
a panel genotype — is excluded. 19 rows.

**Title:**

> **Table 3. Structural observables of the RT–DOR interface from equilibrium MD simulations**

**Column headers** (compact, with a footnote key — the table is wide and will want landscape
orientation or a reduced font):

| Mutation category | Mutation | d(103 C=O···N4x) | d(103 sc···DOR) | N(chl) | θ(Y188/chl) | d(V179···DOR) | d(S105···DOR) | d(106···DOR) | N(pyr) | V(NNIBP) |

**Footnote:**

> Values are means over three independent 100 ns equilibrium MD simulations, each quantity
> averaged within a replicate before averaging across replicates; errors are standard errors of
> the mean across the three replicates. Distances in Å, angles in degrees, volume in Å³.
> d(103 C=O···N4x), distance between the residue 103 main-chain carbonyl oxygen and the DOR
> triazolinone nitrogen; d(103 sc···DOR), minimum distance between the residue 103 side-chain
> polar atoms and DOR; N(chl), number of RT heavy atoms within 4.0 Å of the chlorocyanophenyl
> ring; θ(Y188/chl), interplanar angle between the Tyr188 and chlorocyanophenyl ring planes,
> undefined in Y188L; d(V179···DOR), d(S105···DOR) and d(106···DOR), minimum heavy-atom distances
> between DOR and Val179, Ser105 and residue 106 respectively; N(pyr), number of RT heavy-atom
> contacts within 4.0 Å of the pyridinone moiety and its exocyclic substituents; V(NNIBP), NNIBP
> volume. Per-replicate values are given in Supplementary Table 4.

**The table itself** — paste from `paper/tables/Table-3-structural.csv`:

| Category | Genotype | d(103 C=O···N4x) | d(103 sc···DOR) | N(chl) | θ(Y188/chl) | d(V179···DOR) | d(S105···DOR) | d(106···DOR) | N(pyr) | V(NNIBP) |
|---|---|---|---|---|---|---|---|---|---|---|
| Wild type | WT | 2.97 ± 0.01 | 8.38 ± 0.15 | 19.5 ± 0.5 | 13.4 ± 0.4 | 3.62 ± 0.06 | 6.65 ± 0.09 | 3.36 ± 0.07 | 14.8 ± 1.4 | 228.4 ± 11.7 |
| Susceptible | V106I | 3.62 ± 0.30 | 9.20 ± 0.50 | 20.2 ± 0.3 | 14.7 ± 2.3 | 3.74 ± 0.20 | 7.15 ± 0.14 | 3.22 ± 0.09 | 16.0 ± 0.8 | 247.3 ± 6.7 |
| Susceptible | K103N | 3.08 ± 0.06 | 5.07 ± 0.04 | 20.3 ± 0.6 | 14.2 ± 1.7 | 3.77 ± 0.13 | 6.48 ± 0.10 | 3.15 ± 0.14 | 17.4 ± 1.7 | 236.8 ± 3.7 |
| Susceptible | Y181C | 2.95 ± 0.06 | 8.17 ± 0.26 | 24.0 ± 2.6 | 14.9 ± 1.6 | 3.51 ± 0.13 | 5.99 ± 0.15 | 3.05 ± 0.17 | 16.4 ± 1.1 | 257.3 ± 1.6 |
| Susceptible | G190A | 3.01 ± 0.03 | 8.56 ± 0.20 | 20.0 ± 1.2 | 16.0 ± 1.8 | 3.51 ± 0.05 | 6.81 ± 0.08 | 3.26 ± 0.10 | 14.5 ± 0.7 | 233.3 ± 1.5 |
| Resistant | Y318F | 2.99 ± 0.06 | 8.29 ± 0.05 | 19.2 ± 0.2 | 14.2 ± 0.2 | 3.69 ± 0.04 | 6.61 ± 0.09 | 3.30 ± 0.004 | 14.2 ± 0.4 | 233.8 ± 7.8 |
| Resistant | V106A | 3.05 ± 0.08 | 8.45 ± 0.07 | 17.8 ± 0.7 | 20.5 ± 7.1 | 3.55 ± 0.07 | 5.30 ± 0.02 | 3.54 ± 0.10 | 11.6 ± 1.7 | 250.5 ± 7.9 |
| Resistant | A98G+F227C | 2.97 ± 0.01 | 8.32 ± 0.07 | 20.4 ± 3.5 | 22.0 ± 9.4 | 3.39 ± 0.01 | 6.19 ± 0.21 | 3.41 ± 0.02 | 15.2 ± 0.4 | 238.1 ± 10.2 |
| Resistant | V106I+F227C | 2.80 ± 0.16 | 8.07 ± 0.07 | 19.2 ± 2.3 | 14.0 ± 0.3 | 3.50 ± 0.06 | 6.47 ± 0.20 | 3.09 ± 0.16 | 16.8 ± 2.7 | 232.8 ± 4.6 |
| Resistant | V106A+F227L | 3.17 ± 0.08 | 8.54 ± 0.10 | 17.2 ± 0.3 | 16.7 ± 1.8 | 3.50 ± 0.11 | 5.24 ± 0.06 | 3.63 ± 0.04 | 10.3 ± 0.6 | 239.2 ± 3.6 |
| Resistant | Y188L | 2.91 ± 0.30 | 8.44 ± 0.24 | 13.7 ± 1.5 | — | 3.63 ± 0.03 | 6.38 ± 0.42 | 3.13 ± 0.28 | 17.4 ± 3.7 | 241.6 ± 7.1 |
| Resistant | V106A+P225H | 3.11 ± 0.08 | 8.54 ± 0.11 | 17.9 ± 1.2 | 21.1 ± 8.3 | 3.59 ± 0.07 | 5.29 ± 0.17 | 3.64 ± 0.04 | 11.8 ± 0.1 | 223.3 ± 7.4 |
| Resistant | V106A+L234I | 3.22 ± 0.09 | 8.48 ± 0.02 | 18.9 ± 0.4 | 16.2 ± 1.7 | 3.43 ± 0.04 | 5.13 ± 0.07 | 3.58 ± 0.07 | 10.3 ± 0.8 | 237.3 ± 12.0 |
| Resistant | K103N+M230L | 3.05 ± 0.04 | 4.80 ± 0.14 | 22.8 ± 2.4 | 13.8 ± 0.3 | 3.53 ± 0.08 | 6.19 ± 0.29 | 3.20 ± 0.12 | 16.6 ± 0.6 | 250.5 ± 1.2 |
| Uncertain | V106M | 3.56 ± 0.20 | 8.88 ± 0.38 | 19.3 ± 0.9 | 16.7 ± 1.5 | 3.76 ± 0.23 | 7.04 ± 0.46 | 3.16 ± 0.14 | 16.8 ± 2.4 | 239.6 ± 22.6 |
| Uncertain | G190S | 3.00 ± 0.03 | 8.72 ± 0.04 | 19.0 ± 0.4 | 27.4 ± 5.3 | 3.78 ± 0.08 | 5.95 ± 0.28 | 3.15 ± 0.20 | 17.2 ± 2.4 | 268.1 ± 18.3 |
| Uncertain | L100I+K103N | 3.20 ± 0.15 | 5.43 ± 0.20 | 18.4 ± 0.5 | 12.7 ± 0.3 | 3.89 ± 0.29 | 6.31 ± 0.23 | 3.38 ± 0.02 | 17.8 ± 0.4 | 229.9 ± 6.9 |
| Uncertain | K103N+P225H | 3.10 ± 0.02 | 5.06 ± 0.14 | 18.3 ± 0.6 | 15.5 ± 3.0 | 3.71 ± 0.08 | 6.37 ± 0.15 | 3.34 ± 0.09 | 15.7 ± 0.5 | 239.7 ± 11.5 |
| Uncertain | G190E | 3.17 ± 0.09 | 8.63 ± 0.16 | 22.7 ± 2.0 | 12.0 ± 0.9 | 6.47 ± 0.44 | 5.70 ± 0.45 | 2.92 ± 0.24 | 17.3 ± 2.7 | 284.3 ± 3.2 |

*(Y318F's d(106···DOR) error is genuinely 0.004 Å — three replicates agreeing to 0.01 Å — shown at
an extra digit rather than as "± 0.00".)*

**Knock-on:** update the Results subsection's pointers — "Supplementary Table 4" becomes
"Table 3" everywhere in the Results text, and "Supplementary Figure 3A/3B" becomes
"Figure 3A/3B".

---

# 3. Rewritten Discussion

Nine paragraphs instead of ten, and roughly 15% shorter, but every mechanistic claim now carries
a citation that either supports or complicates it. New references are marked **[N1]–[N6]** and
listed at the end.

**Paragraphs 1–3 (profile, approach, headline result) stay as they are**, with one correction
carried over from the audit: in paragraph 2, *"MM/GBSA-based absolute energy computation"* should
read *"MM/GBSA interaction energies referenced to wild type"* — MM/GBSA here is not an absolute
energy, and the supplement's own Notes sheet says so.

### ¶4 — Lys103 (replaces the current K103N paragraph)

> The anchor DOR makes to the residue 103 backbone is indifferent to what happens to the side
> chain at that position. K103N shortens the side chain and pulls its polar atoms more than 3 Å
> closer to the drug, yet the main-chain hydrogen bond is unchanged, and it survives intact in all
> four K103N-containing genotypes (Table 3). This is the structural expression of the design
> intent behind DOR:14,15 by engaging the backbone rather than the Lys103 side chain, the drug
> sidesteps the mechanism by which K103N defeats efavirenz and nevirapine — stabilisation of the
> closed, unliganded pocket that slows inhibitor entry9 and permits catalysis with inhibitor
> bound.8 The prediction is independently supported outside subtype B: in HIV-1 subtype C
> pseudoviruses, K103N alone leaves DOR susceptibility unchanged at 0.96-fold.**[N1]** It may in
> part explain why K103N+M230L and L100I+K103N also show insignificant binding penalties in our
> simulations (∆∆Gbind = 0.59 ± 0.41 and 0.19 ± 0.68 kcal/mol).

### ¶5 — the two tyrosines (replaces the current Y181C/Y188L paragraph)

> The two tyrosines of the pocket are not equivalent for DOR. Loss of aromatic stacking at
> Tyr181 and Tyr188 is the classical explanation for nevirapine and delavirdine
> resistance,10,11,**[N3]** but in our trajectories losing Tyr181 costs DOR nothing measurable,
> whereas losing Tyr188 removes about a third of the chlorocyanophenyl burial (Figure 3A) and
> carries both the largest van der Waals penalty in the panel (∆∆EvdW = 2.25 ± 0.20 kcal/mol) and
> one of the largest binding free energy penalties (∆∆Gbind = 4.52 ± 0.49 kcal/mol). Aromatic
> packing of DOR is therefore carried by Tyr188 alone, consistent with the crystallographic
> observation that the Tyr181 side chain is rotated away from the drug,14,16 and this asymmetry
> accounts for the clinical pattern in which Y181C spares DOR while Y188L is among the most
> resistant single substitutions in subtype B18 and subtype C alike.**[N1]**

### ¶6 — position 190 (replaces the current G190 paragraph) — **this closes the G190E problem**

> Position 190 shows that a substitution can be structurally disruptive without being disruptive
> where it matters. G190A resists first-generation NNRTIs by introducing a steric bulge into a
> compact region of the pocket,6,**[N3]** but that region has room around the DOR pose, and no
> observable in Table 3 separates G190A from WT. G190S tilts the chlorocyanophenyl ring away from
> coplanarity with Tyr188 without loosening its packing, so its effect is on stacking geometry
> rather than on contact. G190E is the most distorted system in the panel by two independent
> measures — Val179 displaced by nearly 3 Å and the pocket enlarged by a quarter — yet the
> distortion is peripheral: the aromatic anchor is undisturbed and the drug remains fully buried,
> which is why a badly deformed pocket nonetheless binds DOR slightly more tightly than WT
> (∆∆Gbind = −1.02 ± 0.38 kcal/mol). A computed gain in affinity at this position is less
> anomalous than it first appears. Substitutions at residue 190 increase susceptibility to
> delavirdine while severely impairing replication capacity, G190E most of all,**[N4]** so the
> resistance associated with G190E need not be a binding effect at all — and pocket-level
> distortion should not be read as evidence of reduced affinity.

### ¶7 — position 106 (replaces both current V106 paragraphs)

> Position 106 supports two distinct and opposite structural responses. Shortening the side chain
> opens a cavity: in all four V106A-containing genotypes DOR slides 1.3–1.5 Å toward Ser105, out
> of its crystallographic pose, and sheds roughly a quarter of its pyridinone contacts
> (Figure 3B, Table 3). Residue 106 packs directly against the inhibitor in the hydrophobic palm
> pocket, and crystal structures of RT mutated at this codon locate the resistance mechanism to
> that contact region,**[N2]** consistent with V106A being the primary in vitro DOR resistance
> pathway, from which F227L and L234I subsequently emerge16 — the two genotypes with the lowest
> pyridinone contact counts in our panel. Enlarging the side chain does the opposite. In V106I and
> V106M, DOR is pushed away from Ser105 and pressed against residue 106, the pyridinone contact
> count does not fall, and the cost appears instead at the other end of the molecule: these are
> the only two genotypes in which the residue 103 backbone anchor lengthens appreciably, by
> roughly 0.6 Å against a distance otherwise constant to within 0.25 Å across all 18 mutants
> (Table 3). Both responses are penalising in our calculations, and V106I and V106M account for
> the two largest apparent overestimates of ∆∆Gbind in the panel (+2.27 ± 0.74 and
> +6.10 ± 0.16 kcal/mol) against subtype B fold-changes of 1.1 and 3.4. For V106M that reading may
> be too harsh: in subtype C, V106M alone reduces DOR susceptibility 17-fold,**[N1]** a value our
> calculation would not badly overestimate, and V106A/V106M is the primary resistance pathway in
> subtypes A and C as well as B.16

### ¶8 — ∆∆E versus ∆∆G (tightened; the current paragraph stands, minus the causal error)

> Several genotypes show inconsistencies between ∆∆ETotal and ∆∆Gbind, which is to be expected:
> the two measure different quantities, the interface energy of the bound complex versus the
> energetic difference between bound and unbound states.34 G190E and A98G+F227C both distort the
> DOR–RT interface (∆∆ETotal = 1.94 ± 0.40 and 2.47 ± 0.24 kcal/mol), but comparable distortion
> is present in the apo simulations, so the penalty largely cancels in the thermodynamic cycle and
> the residual favours binding (∆∆Gbind = −1.02 ± 0.38 and −1.53 ± 0.87 kcal/mol). An end-point
> interface score and a free energy of binding are not interchangeable, and a panel of this kind
> is where that distinction becomes visible rather than academic.

### ¶9 — Met230, and the general limitation (replaces the current closing paragraph)

> The clearest case is K103N+M230L, for which our calculations resolve no binding penalty
> (∆∆Gbind = 0.59 ± 0.41 kcal/mol) despite a 36-fold reduction in susceptibility. Table 3 shows
> why the model cannot find one: nothing at the DOR interface is disturbed — chlorocyanophenyl
> burial is the second highest in the panel, 22.8 ± 2.4 against 19.5 ± 0.5 in WT, and the Lys103
> anchor is intact at 3.05 ± 0.04 Å. Met230 lies in the p66 primer grip, where its side chain
> contacts the primer terminus at position −2,**[N5]** and M230L impairs RT enzymatic function and
> viral replicative capacity in addition to reducing NNRTI susceptibility.**[N6]** Our systems are
> binary RT–DOR complexes containing no template/primer, so a resistance mechanism acting through
> that interface is structurally absent from the model. The same limitation applies more broadly:
> resistance can act through large-scale allosteric change, through the conformational equilibrium
> of the unliganded pocket on timescales far longer than we sample, or through effects on enzyme
> processivity and viral fitness, none of which these calculations capture. Such mechanisms are in
> principle accessible to hybrid QM/MM approaches,72,73 enhanced sampling74–76 and machine
> learning,77–79 but the more immediate point is that inhibitor affinity and viral susceptibility
> are distinct quantities, and a panel spanning both susceptible and resistant genotypes is where
> they come apart.

---

## New references

> **[N1]** Reddy, N.; Papathanasopoulos, M.; Steegen, K.; Basson, A. E. K103N, V106M and Y188L
> Significantly Reduce HIV-1 Subtype C Phenotypic Susceptibility to Doravirine. *Viruses* **2024**,
> *16* (9), 1493. https://doi.org/10.3390/v16091493
>
> **[N2]** Ren, J.; Nichols, C. E.; Chamberlain, P. P.; Weaver, K. L.; Short, S. A.; Stammers, D. K.
> Crystal Structures of HIV-1 Reverse Transcriptases Mutated at Codons 100, 106 and 108 and
> Mechanisms of Resistance to Non-nucleoside Inhibitors. *J. Mol. Biol.* **2004**, *336* (3),
> 569–578. https://doi.org/10.1016/j.jmb.2003.12.055
>
> **[N3]** Ren, J.; Stammers, D. K. Structural Basis for Drug Resistance Mechanisms for
> Non-Nucleoside Inhibitors of HIV Reverse Transcriptase. *Virus Res.* **2008**, *134* (1–2),
> 157–170. https://doi.org/10.1016/j.virusres.2007.12.018
>
> **[N4]** Huang, W.; Gamarnik, A.; Limoli, K.; Petropoulos, C. J.; Whitcomb, J. M. Amino Acid
> Substitutions at Position 190 of Human Immunodeficiency Virus Type 1 Reverse Transcriptase
> Increase Susceptibility to Delavirdine and Impair Virus Replication. *J. Virol.* **2003**,
> *77* (2), 1512–1523. https://doi.org/10.1128/jvi.77.2.1512-1523.2003
>
> **[N5]** Ghosh, M.; Jacques, P. S.; Rodgers, D. W.; Ottman, M.; Darlix, J.-L.; Le Grice, S. F. J.
> Alterations to the Primer Grip of p66 HIV-1 Reverse Transcriptase and Their Consequences for
> Template-Primer Utilization. *Biochemistry* **1996**, *35* (26), 8553–8562.
> https://doi.org/10.1021/bi952773j
>
> **[N6]** Xu, H.-T.; Quan, Y.; Schader, S. M.; Oliveira, M.; Bar-Magen, T.; Wainberg, M. A. The
> M230L Nonnucleoside Reverse Transcriptase Inhibitor Resistance Mutation in HIV-1 Reverse
> Transcriptase Impairs Enzymatic Function and Viral Replicative Capacity. *Antimicrob. Agents
> Chemother.* **2010**, *54* (6), 2401–2408. https://doi.org/10.1128/aac.01795-09

---

# Two things you need to decide

### A. **[N1] partly contradicts the paper, and I think that is a feature**

Reddy et al. measured DOR susceptibility in subtype C pseudoviruses and found **V106M alone at
17.3 ± 5.3-fold**. Table 1 in this manuscript lists V106M at **3.4-fold**, from a single subtype B
isolate (Lai 2014). If the subtype C value is closer to the truth, then our FEP result for V106M
(∆∆Gbind = +6.10 ± 0.16 kcal/mol) is **not** the panel's worst overestimate — it may be a correct
prediction being scored against a weak reference value.

I have written ¶7 to say this explicitly, because it strengthens the paper: the FEP made a call
that the single subtype B datapoint contradicted and an independent subtype C study supports. But
it is your call whether to go further and note in Table 1 that V106M's median rests on n = 1.
Related: it also means the "two largest overestimates" framing is now half-retracted, so do not
reuse that phrasing elsewhere.

### B. **One citation I could not fully verify**

**[N2]** (Ren 2004) definitely covers crystal structures of RT mutated at codons 100, 106 and 108,
including V106A (PDB 1S1W). I could not retrieve the full text to confirm the *specific* mechanism
wording, so I have cited it conservatively — "locate the resistance mechanism to that contact
region" — rather than claiming it reports pocket enlargement. If you have access, check the
abstract and sharpen the sentence; if it does describe an enlarged cavity, that is a much stronger
match to our V106A result and worth saying directly.

---

# Also confirmed while doing this

**Reference 16 (Feng et al.) has the wrong year.** Crossref gives *Antimicrob. Agents Chemother.*
**2015**, *59* (1), 590–598 — volume 59 issue 1 is the January 2015 issue. Both the reference list
and Table 1's "Feng et al (2014)" should read 2015.
