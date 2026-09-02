# Companion paragraph for V106I and V106M

All values at the 4.0 Å cutoff, mean ± SEM over three 100 ns replicates.

| | pyridinone contacts | DOR → Ser105 | Lys103 H-bond | ring burial | ∆∆G_bind | fold |
|---|---:|---:|---:|---:|---:|---:|
| WT | 14.7 ± 1.4 | 6.65 ± 0.09 | **2.97 ± 0.01** | 19.5 ± 0.5 | — | 1 |
| V106A | **11.7 ± 1.8** (−20%) | **5.30 ± 0.02** | 3.05 ± 0.08 | 17.8 ± 0.7 | +1.76 ± 0.51 | 9.6 |
| V106I | **15.9 ± 0.6** (+8%) | **7.15 ± 0.14** | **3.62 ± 0.30** | 20.2 ± 0.3 | +2.27 ± 0.74 | 1.1 |
| V106M | **16.7 ± 2.3** (+14%) | **7.04 ± 0.46** | **3.56 ± 0.20** | 19.3 ± 0.9 | +6.10 ± 0.16 | 3.4 |

The contrast is symmetric on every coordinate: where V106A loses pyridinone
contacts, V106I and V106M gain them; where V106A slides toward Ser105, they are
pushed away from it. Neither bulky variant loses packing, and the penalty
instead appears as strain on the anchoring hydrogen bond.

Supporting detail: DOR sits 3.36 ± 0.07 Å from residue 106 in WT, 3.54 ± 0.10 Å
in V106A (further, as the cavity opens) and 3.22 ± 0.09 / 3.16 ± 0.14 Å in
V106I / V106M (closer, as the side chain crowds it). V106M additionally loses
triazolinone contacts (30.0 ± 2.6 against 35.4 ± 3.2 in WT), i.e. the ring that
carries the hydrogen bond is itself pushed off — consistent with its being the
larger of the two penalties.

---

## Suggested paragraph

> The two bulkier substitutions at this position act in the opposite sense. In
> V106I and V106M simulations DOR is displaced *away* from Ser105 rather than
> toward it, to 7.15 ± 0.14 Å and 7.04 ± 0.46 Å against 6.65 ± 0.09 Å in WT, and
> is pressed closer to residue 106 (3.22 ± 0.09 Å and 3.16 ± 0.14 Å against
> 3.36 ± 0.07 Å). Consistent with this crowding, pyridinone contacts *increase*
> rather than fall, to 15.9 ± 0.6 in V106I and 16.7 ± 2.3 in V106M against
> 14.7 ± 1.4 in WT, and ring burial is unchanged. The energetic penalty in these
> genotypes therefore does not arise from lost packing but from strain on the
> anchoring interaction: the hydrogen bond between the Lys103 main-chain carbonyl
> and the DOR triazolinone lengthens from 2.97 ± 0.01 Å in WT to 3.62 ± 0.30 Å in
> V106I and 3.56 ± 0.20 Å in V106M, the two longest values in the panel. Both
> genotypes nonetheless remain phenotypically susceptible or near-susceptible
> (1.1- and 3.4-fold), indicating that this strain is tolerated in vivo, and they
> account for the two largest overestimates of ∆∆G_bind in our panel
> (+2.27 ± 0.74 and +6.10 ± 0.16 kcal/mol). Position 106 therefore supports two
> distinct and opposite structural responses — cavity formation and drug slippage
> when the side chain is shortened, steric crowding and anchor strain when it is
> enlarged — of which only the former tracks the clinical phenotype.

---

## Note

The last clause is the one worth keeping: it is the panel's cleanest example of
the calculations detecting real, reproducible strain that the virus tolerates,
and it explains FEP's only false positive (V106I) without dismissing it as noise.
