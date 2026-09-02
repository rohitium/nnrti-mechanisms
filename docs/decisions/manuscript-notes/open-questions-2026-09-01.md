# Two open observations, answered from the data

## 1. Why do V106I and V106M carry large ∆∆G_bind penalties?

**Because they act by the opposite mechanism to V106A — steric crowding rather
than cavity creation — and the strain lands on the Lys103 anchor.**

| | ring burial | Lys103 backbone H-bond | DOR → Ser105 | DOR → residue 106 | ∆∆G_bind | fold |
|---|---:|---:|---:|---:|---:|---:|
| WT | 19.5 ± 0.5 | **2.97 ± 0.01** | 6.65 ± 0.09 | 3.36 ± 0.07 | — | 1 |
| V106A | 17.8 ± 0.7 | 3.05 ± 0.08 | **5.30 ± 0.02** | 3.54 ± 0.10 | +1.76 | 9.6 |
| V106I | 20.2 ± 0.3 | **3.62 ± 0.30** | **7.15 ± 0.14** | 3.22 ± 0.09 | +2.27 | 1.1 |
| V106M | 19.3 ± 0.9 | **3.56 ± 0.20** | **7.04 ± 0.46** | 3.16 ± 0.14 | **+6.10** | 3.4 |

**DOR moves in opposite directions.** Removing the two γ-methyls (V106A) opens a
cavity and the drug *falls into it*, sliding 1.35 Å toward Ser105 and away from
residue 106. Adding bulk (Ile, Met) does the reverse: DOR is pushed 0.4–0.5 Å
*away* from Ser105 and pressed 0.14–0.20 Å closer to residue 106. Same site,
two mechanisms, opposite displacement vectors.

**The penalty is anchor strain, not lost packing.** V106I and V106M do not lose
burial — V106I actually gains slightly (20.2 vs 19.5) — and the NNIBP volume is
barely changed. What changes is the anchoring hydrogen bond: the Lys103
main-chain carbonyl to triazolinone N–H distance stretches from 2.97 ± 0.01 Å to
**3.62 ± 0.30** (V106I) and **3.56 ± 0.20** (V106M), the two longest values in
the entire panel. The bulkier side chain displaces the drug along the pocket
axis, and the anchor takes the strain.

**This also explains the classification failure.** V106I is FEP's only false
positive (∆∆G = +2.27, measured 1.1-fold) and V106M is the panel's largest
∆∆G at 6.10 against only 3.4-fold. The calculation is not wrong about the
physics — it detects a real, reproducible strain on the anchoring hydrogen
bond. The virus simply tolerates it: a stretched hydrogen bond is not a broken
one, and DOR remains packed. This is a concrete instance of binding affinity
and susceptibility diverging, and it is worth stating as such rather than
treating V106I as noise.

Suggested text:

> The two bulkier substitutions at position 106 act by the opposite mechanism to
> V106A. Where the loss of both γ-methyls in V106A opens a cavity into which DOR
> slides toward Ser105, the larger isoleucine and methionine side chains displace
> the drug in the reverse direction, 0.4–0.5 Å further from Ser105 and closer to
> residue 106. Neither loses packing — ring burial is unchanged or slightly
> increased — but both strain the anchoring interaction: the Lys103 main-chain
> hydrogen bond lengthens from 2.97 ± 0.01 Å in WT to 3.62 ± 0.30 Å in V106I and
> 3.56 ± 0.20 Å in V106M, the longest values in the panel. That both genotypes
> nonetheless remain phenotypically susceptible or near-susceptible (1.1- and
> 3.4-fold) indicates that this strain is tolerated in vivo, and accounts for the
> two largest overestimates of ∆∆G_bind in our panel.

---

## 2. Why do ∆∆G_bind and ∆∆E_Total disagree qualitatively?

**Because they are not the same quantity. ∆∆E_Total is a bound-state
interaction energy; ∆∆G_bind is a difference between the bound and unbound
states.** A mutation that costs the same in both states has ∆∆G_bind ≈ 0 no
matter what it does to the interface.

The FEP decomposes into the two legs it actually computes:

| leg | ∆G holo | ∆G apo | **∆∆G = holo − apo** |
|---|---:|---:|---:|
| wt → G190E | −53.84 ± 1.22 | −52.82 ± 1.51 | **−1.02** |
| K103N → K103N+M230L | −14.09 ± 1.94 | −14.14 ± 2.21 | **+0.05** |
| F227C → A98G+F227C | −8.97 ± 0.27 | −7.75 ± 0.19 | **−1.22** |
| wt → Y188L | −0.40 ± 0.35 | −4.92 ± 0.32 | **+4.52** |

**G190E is the clearest case.** Introducing the glutamate charge costs ~53
kcal/mol — in *both* states. The binding effect is the 1 kcal/mol residual, about
2% of either leg. Almost the entire energetic consequence of the mutation is an
intrinsic charging/desolvation cost that the protein pays whether or not DOR is
bound, and which therefore cancels. MM/GBSA never sees the apo state, so it
reports only the bound-state consequence: ∆∆E_Total = +1.94, dominated by
∆∆E_GB = +2.50 — the new charge is poorly solvated at a buried interface.
Both numbers are right about different questions.

Contrast Y188L, where the legs do *not* cancel: the mutation is 4.5 kcal/mol
cheaper in apo than in holo, because removing the tyrosine costs little in the
empty pocket but a great deal when it is stacking against the drug. That is what
a genuine binding effect looks like, and both methods agree on it.

**A second reason, specific to the charged and polar cases.** For G190E and
K103N+M230L the MM/GBSA total is a small residual of two large, opposing terms:

| genotype | ∆∆E_elec | ∆∆E_GB | ∆∆E_Total |
|---|---:|---:|---:|
| G190E | −1.37 ± 0.43 | +2.50 ± 0.32 | +1.94 ± 0.40 |
| K103N+M230L | −2.74 ± 0.56 | +2.33 ± 0.57 | −1.42 ± 0.67 |
| A98G+F227C | +0.30 ± 0.22 | +0.31 ± 0.15 | +2.47 ± 0.24 |

Coulombic attraction and generalised-Born desolvation nearly cancel, so the sign
of ∆∆E_Total is set by an implicit-solvent approximation operating in exactly
the regime where it is least reliable. K103N+M230L's favourable total is a
−2.74 electrostatic gain surviving a +2.33 desolvation penalty; a modest error in
either term flips it.

**A98G+F227C is a third case again.** Its electrostatic and GB terms are both
negligible; the +2.47 is driven by van der Waals (+1.81), i.e. genuinely worse
packing in the bound complex. Yet the FEP finds the mutation 1.22 kcal/mol
*more* favourable in holo than in apo. The two are reconcilable: MM/GBSA compares
interaction energies between two separately equilibrated ensembles and ignores
how the mutation changes the protein's own internal energy, which the alchemical
calculation includes. Neither is measuring interface packing alone.

**What to say in the Discussion.** The disagreements are not a failure of either
method but a consequence of what each omits. ∆∆E_Total is blind to the apo
state, so it cannot distinguish a mutation that damages the interface from one
that damages the protein equally whether or not the drug is present. ∆∆G_bind
sees both but pays for it: it is a small difference between two large numbers,
which is why its error bars are dominated by the legs rather than by the result.
The genotypes where they agree — Y188L, the V106A series — are those where the
mutation acts specifically at the bound interface.

---

## 3. A charge sign error in the current draft

> "charge-changing transitions, e.g. Lys+1 to Asn0 in K103N or **Gly0 to Glu+1**
> in G190E"

Glutamate is **negatively** charged. `ops/slurm/fep/charge_correction.py`
records the legs correctly as `wt_to_K103N: Lys+ -> Asn0` and
`wt_to_G190E: Gly0 -> Glu-`, both with Δq = −1. The draft should read
**"Gly0 to Glu−1"**.

The same paragraph also states that longer switching "was not enough to mitigate
the errors" for G190E. That is now false — G190E's SEM is 0.38, the best in the
panel. See `paper/CHANGES_2026-09-01_g190e_resolved.md`; the corrected
point is that poor overlap did not determine the error, one unconverged
replicate did.
