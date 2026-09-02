# DOR moiety terminology — and a correction

Answers the question of which names to standardise on, and corrects an error I
introduced into the Discussion draft that has since been merged.

---

## 1. The chemistry

Doravirine is
*3-chloro-5-({1-[(4-methyl-5-oxo-4,5-dihydro-1H-1,2,4-triazol-3-yl)methyl]-2-oxo-4-(trifluoromethyl)-1,2-dihydropyridin-3-yl}oxy)benzonitrile*.

It has three rings, joined by an ether oxygen and a methylene:

| ring | chemistry | role |
|---|---|---|
| **chlorocyanophenyl** | benzonitrile bearing Cl — a genuine benzene | stacks with Tyr188 |
| **pyridinone** | 2-oxo-1,2-dihydropyridine bearing CF₃ — **a nitrogen heterocycle, not a benzene** | central ring; packs against Val106 |
| **triazolinone** | 4-methyl-5-oxo-1,2,4-triazole, carries the N–H | donates the H-bond to the Lys103 main chain |

### Recommendation

Use **chlorocyanophenyl**, **pyridinone** (or "central pyridinone ring" on first
use), and **triazolinone** throughout.

**The Introduction's "central phenyl ring in DOR makes a hydrophobic bond with
Val106" is chemically incorrect.** The central ring is a pyridinone, not a
phenyl. The *interaction* described is right — Val106 does contact that
ring — only the name is wrong. Suggested fix:

> …the central **pyridinone** ring in DOR makes hydrophobic contact with Val106…

Calling it "phenyl" is also actively confusing in this paper, because the
molecule contains a real phenyl ring (the chlorocyanophenyl) that the same
paragraph discusses two sentences earlier.

---

## 2. CORRECTION — the Lys103 hydrogen bond is donated by the **triazolinone**, not the pyridinone

The merged Discussion text reads:

> "the residue 103 mainchain carbonyl oxygen accepts a hydrogen bond from the
> **pyridinone nitrogen** in DOR at 2.97 ± 0.01 Å"

**This is wrong, and the error is mine** — it entered in ¶2 of
`discussion-expansion-2026-08-31.md`. The donor is the **triazolinone N–H**.

Three independent lines of evidence:

**(a) Distance.** Lys103 main-chain O to each DOR nitrogen, averaged over the WT
trajectories:

| DOR nitrogen | moiety | distance |
|---|---|---:|
| **N4x** | **triazolinone** | **3.05 Å** |
| N3x | triazolinone | 3.70 Å |
| N5x | triazolinone | 5.11 Å |
| N2x | pyridinone | 6.52 Å |
| N1x | nitrile | 14.06 Å |

The pyridinone nitrogen is 6.5 Å away — not hydrogen bonded to anything.

**(b) Valence.** The pyridinone nitrogen is tertiary: it bears the methylene
linker to the triazolinone. It has no hydrogen and therefore *cannot* donate a
hydrogen bond. This is a structural impossibility, not a matter of geometry.

**(c) The code.** The quantity plotted and quoted is
`res103_bb_to_triazN` in `compute_mechanism_coordinates.py`, computed against
`out["triaz_N"]` — the **triazolinone** nitrogens. The number 2.97 ± 0.01 Å is
correct; only its description was wrong.

### Fix

> …the residue 103 main-chain carbonyl oxygen accepts a hydrogen bond from the
> **triazolinone N–H** of DOR at 2.97 ± 0.01 Å in WT, and this distance is
> unchanged at 3.08 ± 0.06 Å in K103N.

Nothing else in that paragraph changes: the distances, their invariance across
K103N backgrounds, and the conclusion all stand.

### Knock-on: the V106A framing

In the moiety analysis I described the pyridinone as *"the ring that carries the
Lys103 backbone hydrogen bond"* and framed the V106A slide as *"displacement of
the anchored end of the ligand."* Both follow from the same error and should be
dropped. The correct reading is simpler and more direct: **Val106 packs against
the pyridinone**, so V106A removes packing from precisely the ring it contacts.
Use:

> …localises the loss to the central pyridinone ring — the ring that Val106
> packs against directly — whose contacts fall from 24.0 ± 2.1 in WT to …

---

## 3. Supporting contact map (WT)

Atom-pair contacts within 4.5 Å between key residues and each moiety, averaged
over three 100 ns WT replicates:

| residue | chlorocyanophenyl | pyridinone | triazolinone |
|---|---:|---:|---:|
| Tyr188 | **27.5** | 0.0 | 0.0 |
| Trp229 | 16.9 | 0.0 | 0.0 |
| Phe227 | 9.0 | 1.3 | 1.0 |
| Tyr318 | 0.0 | 8.9 | **18.6** |
| Lys103 | 0.0 | 1.9 | **10.9** |
| Val106 | 0.4 | **8.9** | **8.9** |
| Tyr181 | 0.3 | 0.1 | 0.0 |
| Gly190 | 0.0 | 0.1 | 0.0 |
| Ser105 | 0.0 | 0.0 | 0.6 |

Three things worth using in the text:

- **Tyr188's engagement is exclusively with the chlorocyanophenyl ring** — 27.5
  contacts there and exactly zero with either other ring. The "Tyr188 alone
  stacks with the chlorocyanophenyl ring" claim in the Introduction is confirmed
  quantitatively.
- **Tyr181 makes 0.4 contacts with the entire ligand.** This is stronger than
  what the Discussion currently claims. Y181C removes a residue that was never
  engaged, which is the cleanest possible statement of why DOR is unaffected.
- **Gly190 makes 0.1 contacts.** Position 190 has essentially no direct contact
  with DOR, so G190S and G190E must act indirectly — which is exactly what the
  position-190 analysis found (G190S tilts the Tyr188 ring; G190E displaces
  Val179 and expands the pocket). This is a useful consistency check to state.

---

## 4. Figure 1B

`results/plots/figure1B_dor_schematic.pdf` (and `.png`, 300 dpi). 2D structure
with the three moieties shaded and labelled, each key residue placed beside the
moiety it contacts and annotated with its measured contact count, the
non-contacting residues (Tyr181, Gly190) shown greyed with dotted leaders, and
the anchoring N–H···O=C Lys103 hydrogen bond drawn as a dashed red line from the
donor nitrogen, which is circled.

Regenerate with:

```bash
PYTHONPATH=. python -m src.analysis.cli.plot_dor_schematic
```

Suggested caption:

> **Figure 1B.** Doravirine, its three ring moieties, and the NNIBP residues that
> engage them. Contacts are protein–ligand heavy-atom pairs within 4.5 Å averaged
> over three independent 100 ns wild-type simulations. Tyr188 engages the
> chlorocyanophenyl ring exclusively; Val106 packs against the central pyridinone
> ring; the triazolinone N–H donates a hydrogen bond to the Lys103 main-chain
> carbonyl (dashed), an interaction independent of the residue-103 side chain.
> Tyr181 and Gly190 (grey) make no appreciable direct contact with the drug.

Note the figure uses the atom-pair convention throughout, consistent with the
rest of the manuscript — see the "protein heavy atoms are atom pairs" correction
in `discussion-expansion-2026-08-31.md`.
