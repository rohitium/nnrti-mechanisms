# Contact cutoff: what the literature uses, and what changes at 4.0 Å

## 1. Is there a standard?

No. Several conventions coexist, and the right one depends on what is being
counted:

| cutoff | typical use | examples |
|---|---|---|
| **4.0 Å** | protein–ligand **hydrophobic contact** | PLIP, LigPlot+ (default max 3.9 Å) |
| **4.5 Å** | heavy-atom contact in MD contact analyses | widely used for contact maps and occupancy |
| **5.0 Å** | defining "binding-site residues" | pocket/shell definitions |
| 6–8 Å | coarse residue-level contact maps | usually Cβ–Cβ, not heavy-atom |

The physical anchor: two carbons in van der Waals contact sit at the sum of
their vdW radii, ≈ **3.4 Å**. So 4.0 Å is "vdW contact plus a small tolerance",
and 4.5 Å is roughly the first coordination shell. Both are defensible; **4.0 Å
is the more common convention for ligand–protein contacts specifically**, while
4.5 Å is more common in MD trajectory analysis, which is what this paper does.

Because there is no standard, the question that matters is not which value is
"right" but **which conclusions survive the choice**. That is measured below.

Reproduce with:

```bash
PYTHONPATH=. python -m src.analysis.cli.sweep_contact_cutoff
```

Output: `results/analysis/mechanisms/contact_cutoff_sweep.csv`.

---

## 2. Sensitivity, % change vs WT

### Pyridinone — the headline result. **Robust.**

| genotype | 3.5 Å | **4.0 Å** | **4.5 Å** | 5.0 Å |
|---|---:|---:|---:|---:|
| V106A | −29.5% | **−22.2%** | **−22.6%** | −16.3% |
| V106A+F227L | −56.7% | **−31.8%** | **−26.2%** | −18.8% |
| V106A+L234I | −47.1% | **−27.7%** | **−25.0%** | −17.1% |
| V106A+P225H | −35.6% | **−21.5%** | **−19.7%** | −14.2% |
| Y188L | +61.4% | +16.7% | +5.4% | +1.3% |

The pyridinone loses ~20–30% of its contacts in every V106A genotype at every
cutoff from 4.0 to 5.0 Å, and it is always the largest loss of any moiety.
**4.0 Å slightly strengthens the result** (−21.5 to −31.8% against −19.7 to
−26.2%). The claim does not depend on the threshold.

Note Y188L is *unchanged or slightly increased* at the pyridinone — confirming
the V106A signature is specific and not a general loosening.

### Chlorocyanophenyl. **Robust.**

| genotype | 3.5 Å | **4.0 Å** | **4.5 Å** | 5.0 Å |
|---|---:|---:|---:|---:|
| V106A | +11.0% | −6.8% | −7.1% | −4.2% |
| V106A+F227L | −43.8% | −13.2% | −8.6% | −7.5% |
| V106A+L234I | −9.4% | −7.3% | −5.4% | −5.1% |
| V106A+P225H | −37.2% | −10.5% | −3.2% | −2.0% |
| **Y188L** | +29.9% | **−20.2%** | **−20.0%** | −16.7% |

Y188L's loss is essentially identical at 4.0 and 4.5 Å (−20.2% vs −20.0%), so
**Figure 3A and the "loss of roughly a quarter" claim are cutoff-independent**.
The V106A genotypes lose consistently less here than at the pyridinone at every
cutoff ≥ 4.0 Å.

### Triazolinone. **NOT robust — the claim must be dropped.**

| genotype | 3.5 Å | **4.0 Å** | **4.5 Å** | 5.0 Å |
|---|---:|---:|---:|---:|
| V106A | −13.0% | **−2.5%** | **+1.2%** | +3.4% |
| V106A+F227L | −43.1% | **−9.7%** | **−1.3%** | +2.8% |
| V106A+L234I | −30.2% | **−3.9%** | **+0.4%** | +3.9% |
| V106A+P225H | −34.5% | **−7.9%** | **+0.1%** | +6.3% |

**The sign flips between 4.0 and 4.5 Å.** My earlier statement that "the distal
triazolinone ring gains contacts, indicating that DOR pivots about its distal
end rather than withdrawing bodily" is therefore **not supportable and should be
removed**. At 4.5 Å the change is within noise of zero (+1.2 to −1.3%); at
4.0 Å it is a small loss. Neither supports a pivot interpretation.

The surviving statement is the one that matters and is fully robust: *the loss
is concentrated at the pyridinone*.

### 3.5 Å is unusable

Every moiety changes sign somewhere in the 3.5 Å column (chlorocyanophenyl is
+11.0% for V106A but −43.8% for V106A+F227L; Y188L is +29.9%). At this cutoff
the counts fall to 2.6–5.9 pairs, so relative changes are dominated by counting
noise. Do not go below 4.0 Å.

---

## 3. Recommendation

**Keep 4.5 Å**, and add the 4.0 Å column to the SI as a sensitivity check.

Reasons: every conclusion the paper draws is unchanged at 4.0 Å, so switching
buys nothing scientifically; 4.5 Å is the more common convention in MD
trajectory analysis specifically; and it is already used consistently throughout
the draft, in Figure 3A, and in the 224-contact figure, so a change would mean
re-deriving numbers you have already merged for no gain.

**But drop the triazolinone claim regardless of which cutoff is chosen** — that
one is genuinely threshold-dependent.

If you would rather standardise on 4.0 Å, everything can be regenerated from the
scripts; say the word and I will re-derive Figure 3A, the burial values, the
Figure 1B contact map and the moiety table together so nothing is left mixed.

---

## 4. A definition inconsistency found while doing this

`compute_dor_moiety_contacts.py` added exocyclic substituents (Cl, C≡N) to the
**chlorocyanophenyl** ring but treated the pyridinone and triazolinone as
**ring atoms only** — so the CF₃ and carbonyl of the pyridinone, and the methyl
and carbonyl of the triazolinone, were excluded. That is ad hoc, and it is why
two sets of absolute numbers appear in my notes:

| moiety | ring atoms only | ring + own substituents |
|---|---:|---:|
| pyridinone, WT, 4.5 Å | 24.0 | 38.5 |
| triazolinone, WT, 4.5 Å | 43.3 | 75.0 |

`sweep_contact_cutoff.py` treats all three rings consistently. **The pyridinone
conclusion holds under both definitions** (ring-only −20.8 to −27.5%;
ring+substituents −19.7 to −26.2%), so nothing scientific turns on it — but the
two scripts should agree before any of these numbers reach the manuscript.

The triazolinone "gain" was partly an artifact of this inconsistency, on top of
being cutoff-dependent, which is a further reason to drop it.

**Decision needed:** ring-only or ring+substituents. I would use
**ring + own substituents** — it is symmetric, and the CF₃ group is a large part
of the pyridinone's actual contact surface — but that changes the absolute
numbers in the suggested V106A replacement text from "24.0 ± 2.1 → 17.8 ± 1.3"
to "38.5 → 29.8". The percentages barely move either way.
