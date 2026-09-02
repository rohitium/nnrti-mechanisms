# Paragraph-level corrections against the current data (2026-09-01)

Four Results/Discussion paragraphs checked number by number against the
regenerated analyses. **Every contact count in the draft is stale**: the cutoff
moved 4.5 Å → 4.0 Å, and the moiety definition became symmetric (each ring now
carries its own exocyclic substituents, where before only the phenyl did).
Distances, angles, volumes and energies are unaffected.

---

## ¶1 — K103N backbone hydrogen bond

**Factual error.** "…accepts a hydrogen bond from the **pyridinone nitrogen**"
— the donor is the **triazolinone N–H**. Lys103:O to triazolinone N4x is
2.90 Å; to the pyridinone N it is 6.52 Å, and that nitrogen is tertiary (it
carries the methylene linker) so it has no hydrogen to donate. The quantity
quoted is `res103_bb_to_triazN`.

**Overreach.** "…none of these backgrounds produces a large binding penalty" —
K103N+P225H does: ∆∆G_bind = 1.76 ± 0.79, i.e. 0.97 above zero after one SEM,
clearing the manuscript's own 0.5 kcal/mol rule. Also only K103N+M230L among
these is clinically Resistant; the other two are Uncertain.

Suggested ending: *"…This may explain why K103N+M230L, despite a 36-fold
reduction in susceptibility, shows no resolved binding penalty in our
simulations (∆∆G_bind = 0.59 ± 0.41 kcal/mol)."*

All distances quoted (2.97 ± 0.01, 3.08 ± 0.06, 8.38 ± 0.15 → 5.07 ± 0.04,
3.05–3.20 Å) are exact.

---

## ¶2 — Y181C / Y188L

- Burial numbers stale and internally mixed. Current at 4.0 Å: **WT 19.5 ± 0.5**,
  **Y181C 24.0 ± 2.6**, **Y188L 13.7 ± 1.5**; the Y188L loss is **30%**.
  Supplementary Figure 3A has been regenerated at 4.0 Å, so the text as written
  no longer matches its own figure.
- "nearly identical to WT" → Y181C is now **1.23× WT** (1.7 σ, not significant).
  Use **"not reduced relative to WT"**.
- "RT heavy atoms" → **"RT heavy-atom contacts"** (these are atom pairs).
- "Aromatic packing … entirely determined by Tyr188" → narrow to **face-to-face
  stacking**. Trp229 also contacts the ring (6.4 contacts) but meets it edge-on
  at **82.3 ± 0.5°**, against Tyr188's **13.4 ± 0.4°**.
- Correct as written: ∆∆E_vdW = 2.25 ± 0.20 **is** the panel's largest;
  "one of the largest" ∆∆G_bind is correctly hedged (V106M is larger at 6.10).

---

## ¶3 — Position 190

All numbers verified: Tyr188 angle 13.4 ± 0.4° → 27.4 ± 5.3°; Val179→DOR
3.62 ± 0.06 → 6.47 ± 0.44 Å; NNIBP volume 230 ± 11 → 286 ± 3 Å³ (draft says
± 12; trivial). G190A is WT-like on every coordinate, supporting "has room to
accommodate it".

Two additions worth making:

**G190S is bimodal, not tilted.** 37% of frames remain below 15° and 32% exceed
40°, with per-replicate means 34.2° / 17.0° / 30.9°. The rotated population also
withdraws (centroid 4.23 → 5.05 Å). The stack is intermittently *broken*, and
the ± 5.3 is real replicate heterogeneity rather than noise.

**G190E's distortion is confined to the pocket, not the drug.** Its Tyr188 stack
is *tighter* than WT (12.0 ± 0.9°) and its ring burial higher (22.7 ± 2.0),
which is what reconciles "major distortions in the NNIBP" with the new
∆∆G_bind = −1.02 ± 0.38. Suggested closing sentence:

> Notably these distortions are confined to the pocket rather than the drug:
> DOR's aromatic anchoring is unaffected in G190E (Tyr188 interplanar angle
> 12.0 ± 0.9°, ring burial 22.7 ± 2.0 contacts, both comparable to or tighter
> than WT), consistent with the slightly favourable ∆∆G_bind computed for this
> genotype and indicating that its 18-fold reduction in susceptibility is not
> mediated by loss of binding affinity.

---

## ¶4 — V106A slide  ← **replacement text below**

Stale on both counts (cutoff and moiety definition):

| genotype | draft (4.5 Å, ring only) | **current (4.0 Å, symmetric)** | change vs WT |
|---|---:|---:|---:|
| WT | 24.0 ± 2.1 | **14.7 ± 1.4** | — |
| V106A | 17.8 ± 1.3 | **11.7 ± 1.8** | **−20%** |
| V106A+F227L | 18.2 ± 1.3 | **10.3 ± 0.7** | **−30%** |
| V106A+L234I | 17.4 ± 1.7 | **10.5 ± 0.8** | **−29%** |
| V106A+P225H | 19.0 ± 0.2 | **11.8 ± 0.1** | **−20%** |

Also "~1.3 Å" understates the slide slightly: it is 1.35, 1.41, 1.52 and 1.36 Å
— **1.3–1.5 Å**.

**The loss is specific to the pyridinone**, which is the point the paragraph is
really making and currently leaves implicit:

| moiety | V106A | +F227L | +L234I | +P225H |
|---|---:|---:|---:|---:|
| **pyridinone** | **−20%** | **−30%** | **−29%** | **−20%** |
| chlorocyanophenyl | −4% | −11% | −4% | −8% |
| triazolinone | +1% | −7% | −3% | −6% |
| whole ligand | −4% | −9% | −4% | −5% |

The whole-ligand figure (−4 to −9%) understates the effect three- to sevenfold.

### Suggested replacement paragraph

> Interestingly, all V106A-containing genotypes considered here demonstrate DOR
> sliding 1.3–1.5 Å toward Ser105 (Supplementary Figure 3B), which suggests DOR
> "slipping" out of its crystallographic pose. Resolving the interface by ligand
> moiety shows that this displacement is not a uniform loosening but is
> concentrated on the central pyridinone ring — the ring that Val106 packs
> against directly — whose heavy-atom contacts within 4.0 Å fall from 14.7 ± 1.4
> in WT to 11.7 ± 1.8 in V106A, 10.3 ± 0.7 in V106A+F227L, 10.5 ± 0.8 in
> V106A+L234I and 11.8 ± 0.1 in V106A+P225H, a loss of 20–30% in every case
> against 1–11% for the other two ring systems. Notably, Val106 contacts the
> pyridinone through its two γ-methyl carbons (3.75 and 3.80 Å, against 4.49 Å
> for Cβ), which alanine does not possess, so V106A removes the very atoms
> making this contact. That a comparable loss is common to all four genotypes,
> which span a 16-fold range of measured susceptibility, suggests a shared
> structural mechanism of DOR resistance in these viruses.

**Do not** restore the earlier claim that the triazolinone *gains* contacts —
it does not survive at 4.0 Å (+1% to −7%) and was partly an artifact of the
asymmetric moiety definition.
