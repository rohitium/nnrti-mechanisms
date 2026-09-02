# New Results subsection + the Discussion trims it forces

Insert after the correlation paragraph that closes **"RT-DOR binding free energy analyses"**
(the one ending "…∆∆ETotal correlation with log10(Fold-change) is weaker still (R2 = 0.13,
p = 0.23)."), immediately before the **Discussion** heading.

Every number below is from `Supplementary-Table-4.xlsx` and was checked against it. The
resolution criterion is stated once in the first paragraph and applied consistently; I verified
each "resolved" and "not resolved" claim arithmetically rather than by eye.

---

## Structural consequences of DOR-associated substitutions

To identify the structural changes underlying these energetic shifts, we measured a set of
geometric observables in every equilibrium trajectory: the hydrogen bond between the residue 103
main-chain carbonyl and DOR, the burial and stacking geometry of the chlorocyanophenyl ring, the
minimum distances from DOR to Ser105, Val179 and residue 106, the number of contacts made by the
pyridinone moiety, and the NNIBP volume. Each quantity was averaged within a replicate before
averaging across the three replicates, so the reported uncertainty is the replicate-to-replicate
standard error; per-replicate values for every simulated system are given in Supplementary Table 4. In
what follows we describe a difference from WT as resolved when it exceeds twice the two standard
errors combined in quadrature.

The hydrogen bond between the residue 103 main-chain carbonyl oxygen and the DOR triazolinone
nitrogen is the most conserved feature of the interface. It measures 2.97 ± 0.01 Å in WT and lies
between 2.80 and 3.22 Å in 16 of the 18 mutant genotypes, including all four genotypes carrying
K103N (3.05–3.20 Å); the two exceptions are V106I (3.62 ± 0.30 Å) and V106M (3.56 ± 0.20 Å).
Substituting lysine by asparagine at position 103 shortens the side chain and brings its polar
atoms much closer to DOR, from 8.38 ± 0.15 Å in WT to between 4.80 and 5.43 Å in every
K103N-containing genotype, against 8.07–9.20 Å in every genotype that retains Lys103, while
leaving the backbone contact itself unchanged.

Burial of the chlorocyanophenyl ring, counted as RT heavy-atom contacts within 4.0 Å, is
19.5 ± 0.5 in WT. The largest reduction in the panel, by a wide margin, occurs in Y188L
(13.7 ± 1.5; Supplementary Figure 3A); the next lowest value, in V106A+F227L, is 17.2 ± 0.3.
Y181C is not reduced (24.0 ± 2.6), and the three highest values in the panel are those of Y181C,
K103N+M230L (22.8 ± 2.4) and G190E (22.7 ± 2.0). The interplanar angle between the Tyr188 and
chlorocyanophenyl ring planes is 13.4 ± 0.4° in WT and is undefined in Y188L, which has no Tyr188
ring. G190S shows the largest departure from WT (27.4 ± 5.3°) and is the only genotype in the
panel whose mean angle differs from WT by more than 10°; several genotypes with intermediate mean
angles (V106A, A98G+F227C and V106A+P225H, 20–22°) carry replicate errors of 7–9° and are not
resolved.

Substitutions at position 106 displace DOR in opposite directions depending on whether the side
chain is shortened or enlarged. In all four V106A-containing genotypes DOR moves toward Ser105,
from 6.65 ± 0.09 Å in WT to between 5.13 and 5.30 Å (Supplementary Figure 3B), moves away from
residue 106 (3.36 ± 0.07 Å in WT to 3.54–3.64 Å), and loses pyridinone contacts, from 14.8 ± 1.4
in WT to between 10.3 and 11.8. In V106I and V106M the displacement is reversed: DOR moves away
from Ser105 (7.15 ± 0.14 and 7.04 ± 0.46 Å) and closer to residue 106 (3.22 ± 0.09 and
3.16 ± 0.14 Å), while pyridinone contacts are unchanged or slightly higher (16.0 ± 0.8 and
16.8 ± 2.4), a shift that is not resolved against replicate variation. These are also the two
genotypes in which the residue 103 anchor lengthens.

The three substitutions at position 190 behave differently from one another. G190A leaves the DOR
pose essentially unperturbed: burial (20.0 ± 1.2), stacking angle (16.0 ± 1.8°), the distances to
Ser105, Val179 and residue 106, and the NNIBP volume (233 ± 2 Å³) are all within replicate error
of the corresponding WT values. G190S produces the largest stacking distortion in the panel while
leaving burial unchanged (19.0 ± 0.4). G190E is distinct from every other genotype in two
respects: the minimum distance from Val179 to DOR increases from 3.62 ± 0.06 Å in WT to
6.47 ± 0.44 Å, more than 2.5 Å beyond the next largest value in the panel (3.89 ± 0.29 Å, in
L100I+K103N), and the NNIBP volume expands from 228 ± 12 Å³ in WT to 284 ± 3 Å³, the largest of
any genotype. Y181C (257 ± 2 Å³) is the only other genotype whose pocket volume is resolved above
WT; G190S (268 ± 18 Å³) has a comparable mean but a replicate error that spans the WT value. The
G190E pocket expansion is not accompanied by any loss of aromatic packing: its stacking angle
(12.0 ± 0.9°) is marginally closer to coplanar than WT and its chlorocyanophenyl burial
(22.7 ± 2.0) is above the WT value.

---

# Discussion trims this forces

The point of the subsection is that these numbers stop being new in the Discussion. Four
Discussion paragraphs currently *report* them; they should now *interpret* them. Rewrites below —
each is shorter, and each says something the Results paragraph does not.

### Trim 1 — K103N paragraph

**Find:**
> In our equilibrium MD trajectories, the residue 103 mainchain carbonyl oxygen accepts a hydrogen bond from a triazolinone nitrogen in DOR at 2.97 ± 0.01 Å in WT, and this distance is almost unchanged at 3.08 ± 0.06 Å in K103N. In fact, substituting lysine by asparagine shortens the side chain and brings its polar atoms closer to DOR (8.38 ± 0.15 Å to 5.07 ± 0.04 Å). Moreover, the same backbone hydrogen bond is preserved in every K103N-containing genotype we simulated (3.05–3.20 Å across K103N+M230L, K103N+P225H and L100I+K103N). This may (in part) explain why K103N+M230L and L100I+K103N also show an insignificant free energy binding penalty, ∆∆Gbind, in our simulations (0.59 ± 0.41 and 0.19 ± 0.68 kcal/mol).

**Replace:**
> The anchor DOR makes to the residue 103 backbone is indifferent to what happens to the side chain at that position. K103N shortens the side chain and pulls its polar atoms more than 3 Å closer to the drug, yet the main-chain hydrogen bond is unchanged, and it survives intact in all four K103N-containing genotypes we simulated. This is the structural expression of the design intent behind DOR: by binding the backbone rather than the Lys103 side chain, the drug is insensitive to substitution at a position that abolishes efavirenz and nevirapine activity. It may in part explain why K103N+M230L and L100I+K103N also show an insignificant binding free energy penalty in our simulations (∆∆Gbind = 0.59 ± 0.41 and 0.19 ± 0.68 kcal/mol).

### Trim 2 — Y181C / Y188L paragraph

**Find:**
> In simulations of the Y181C genotype, DOR packing in the NNIBP is not reduced relative to WT: 24.0 ± 2.6 RT heavy atom contacts within 4.0 Å of the chlorocyanophenyl ring against 19.5 ± 0.5 in WT. In contrast, the most prominent van der Waals penalty on the RT-DOR interface (∆∆EvdW = 2.25 ± 0.20 kcal/mol) and one of the largest binding free energy penalties (∆∆Gbind = 4.52 ± 0.49 kcal/mol) occurs in the well-known DOR Resistant mutation, Y188L. Correspondingly, equilibrium MD simulations reveal a remarkable drop in chlorocyanophenyl ring packing efficiency falling from 19.5 ± 0.5 heavy atom contacts in WT to 13.7 ± 1.5 in Y188L (Supplementary Figure 3A).

**Replace:**
> The two tyrosines of the NNIBP are not equivalent for DOR. Losing Tyr181 costs the drug nothing measurable — chlorocyanophenyl burial in Y181C is if anything above the WT value — whereas losing Tyr188 removes roughly a third of that burial, and Y188L carries both the largest van der Waals penalty on the interface in the whole panel (∆∆EvdW = 2.25 ± 0.20 kcal/mol) and one of the largest binding free energy penalties (∆∆Gbind = 4.52 ± 0.49 kcal/mol).

*(The following sentence, "Aromatic packing of DOR is therefore largely determined by Tyr188…",
already interprets and should stay as it is.)*

### Trim 3 — position 190 paragraph

**Find:**
> While G190A confers resistance to first-generation NNRTIs by introducing a steric bulge into a compact region of the pocket,6 we find the DOR pose in MD simulations of G190A RT-DOR complexes has room to accommodate it. In contrast, the interplanar angle between the Tyr188 and chlorocyanophenyl ring flips from 13.4 ± 0.4° in WT to 27.4 ± 5.3° in G190S simulations, suggesting distorted aromatic packing. Notably, the charge switch in case of G190E simulations causes major distortions in the NNIBP: minimum distance from Val179 to DOR increases from 3.62 ± 0.06 Å to 6.47 ± 0.44 Å, and the NNIBP volume expands from 228 ± 12 Å³ to 284 ± 3 Å³ (Supplementary Table 4).

**Replace:**
> Position 190 illustrates that a substitution can be structurally disruptive without being disruptive where it matters. G190A confers resistance to first-generation NNRTIs by introducing a steric bulge into a compact region of the pocket,6 but the pocket has room to accommodate that bulge around the DOR pose, and none of our observables separates G190A from WT. G190S tilts the chlorocyanophenyl ring away from coplanarity with Tyr188 without loosening its packing, so its effect is on stacking geometry rather than on contact. G190E is the most distorted system in the panel by two independent measures — Val179 displaced by nearly 3 Å and the pocket enlarged by a quarter — yet that distortion is confined to the periphery: the aromatic anchor is undisturbed and the drug remains fully buried. This is why a badly deformed pocket nonetheless binds DOR slightly more tightly than WT (∆∆Gbind = −1.02 ± 0.38 kcal/mol), and it is a caution against reading pocket-level distortion as evidence of reduced affinity.

**This closes B6** — the "major distortions in the NNIBP" claim was previously left unreconciled
against a negative ∆∆G. The reconciliation is that the distortion is peripheral, and the Results
subsection now supplies the numbers that show it.

### Trim 4 — V106A and V106I/V106M paragraphs

**Find (V106A paragraph):**
> Interestingly, all V106A containing genotypes considered here demonstrate DOR sliding 1.3–1.5 Å toward Ser105 (Supplementary Figure 3B), which suggests DOR "slipping" out of its crystallographic pose. In addition, we observe a drop in the number of pyridinone moiety contacts within 4.0 Å from 14.8 ± 1.4 in WT to 11.6 ± 1.7 in V106A, 10.3 ± 0.6 in V106A+F227L, 10.3 ± 0.8 in V106A+L234I, and 11.8 ± 0.1 in V106A+P225H (Supplementary Table 4).

**Replace:**
> The four V106A-containing genotypes share one structural signature: DOR slides 1.3–1.5 Å toward Ser105, out of its crystallographic pose, and sheds roughly a quarter of its pyridinone contacts in the process.

**Find (opening of the V106I/V106M paragraph):**
> In V106I and V106M simulations DOR slides in the opposite direction, i.e. away from Ser105 rather than toward it, to 7.15 ± 0.14 Å and 7.04 ± 0.46 Å against 6.65 ± 0.09 Å in WT, and is pressed closer to residue 106 (3.22 ± 0.09 Å and 3.16 ± 0.14 Å against 3.36 ± 0.07 Å). Pyridinone contacts are correspondingly unchanged or slightly higher, 16.0 ± 0.8 in V106I and 16.8 ± 2.4 in V106M against 14.8 ± 1.4 in WT, a shift that is not resolved against replicate variation. Importantly, the residue 103 main-chain carbonyl to triazolinone nitrogen distance, 2.97 ± 0.01 Å in WT and within 0.25 Å of that value in every other genotype we simulated, lengthens to 3.62 ± 0.30 Å in V106I and 3.56 ± 0.20 Å in V106M — the two longest in the panel (Supplementary Table 4).

**Replace:**
> Enlarging the same side chain does the opposite. In V106I and V106M, DOR is pushed away from Ser105 and pressed against residue 106, and the pyridinone contact count does not fall. The cost instead appears at the other end of the molecule: these are the only two genotypes in the panel in which the residue 103 backbone anchor lengthens appreciably, by roughly 0.6 Å against a distance that is otherwise constant to within 0.25 Å across all 18 mutants.

---

# Two knock-on items

**1. Supplementary Figure 3 citations move.** Panels 3A and 3B are now first cited in Results
(paragraphs 3 and 4), which is where they belong. Remove the duplicate citations from the
Discussion — I have already dropped them from the replacement text above.

**2. Supplementary Table 4 needs a caption.** Supp Tables 1–3 currently ship as bare .xlsx with no
legends anywhere (audit item C8). At minimum add to the Supplementary Text:

> **Supplementary Table 4.** Per-replicate structural observables for every simulated system.
> Nine geometric quantities measured in every production trajectory, averaged within each
> replicate. The Summary sheet gives the mean and replicate-to-replicate standard error per
> genotype; the Per-replicate sheet gives the underlying values with the number of frames
> analysed; the Definitions sheet documents each quantity, the aggregation convention and the
> periodic-boundary treatment.

---

# Still pending

**Met230 / primer grip in the Discussion** — you said you want this, and it is not in this file.
The drafted replacement (with the Ghosh 1996 and Xu 2010 citations) is in
`EDITS_2026-09-01_round2.md` under **B7**. It fits naturally at the end of Trim 1's paragraph or
as its own paragraph after it. Say the word and I will place it against whatever the Discussion
looks like once these trims are in.

One thing that makes it land better now: the Results subsection reports that K103N+M230L has the
second-highest chlorocyanophenyl burial in the panel (22.8 ± 2.4) and an intact anchor
(3.05 ± 0.04 Å). So the Discussion can state that *nothing at the binding interface is disturbed
in a genotype that is 36-fold resistant* as a fact already established in Results, rather than
introducing it as new evidence.
