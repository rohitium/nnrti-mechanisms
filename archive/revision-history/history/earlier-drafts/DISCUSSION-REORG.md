# Discussion: reorganisation, and one claim that needs your decision

Updated 02-Nov-2026 after the apo tables were regenerated. Two parts: a claim
the new numbers do not support, and the rearrangement you asked for.

---

## Part 1 — the cancellation argument in [114] is not supported

Paragraph [114] currently says:

> G190E and A98G+F227C both distort the DOR–RT interface (∆∆ETotal = 1.94 ±
> 0.40 and 2.47 ± 0.24 kcal/mol), **but comparable distortion is present in the
> apo simulations, so the penalty largely cancels in the thermodynamic cycle**
> and the residual favors binding.

Counting features whose 95% interval excludes zero, in each state:

| | holo | apo |
|---|---|---|
| **G190E** | **8** resolved — pocket volume +55.8 ± 12.1 Å³, Val179 +2.84 ± 0.44 Å, and RMSF in the pocket lining, primer grip, β6, fingers and palm | **0** |
| **A98G+F227C** | 2 resolved — Tyr181 χ1 −75.8 ± 12.0°, primer-grip RMSF +0.09 ± 0.03 Å | 2 resolved, but **different features** — palm helix content, connection RMSF |

For G190E the distortion is a bound-state phenomenon and is absent from the
unliganded enzyme, so there is nothing for it to cancel against. For
A98G+F227C the two states are each perturbed, but not in the same way.

Figure 4 says the same thing from the other direction: G190E's induced fit
differs from wild type by **+71 ± 26 Å³**, which is a *discordant* change —
precisely the case that should NOT cancel.

So the mechanism and the observation now point opposite ways. G190E distorts
the bound complex more than any genotype in the panel, that distortion does
not cancel, and ΔΔGbind still comes out favourable (−1.02 ± 0.38). That is a
sharper puzzle than the one the paragraph currently poses, and resolving it is
an interpretation — yours, not mine. Three options, with no recommendation
from me:

1. Drop the cancellation sentence and state the puzzle plainly.
2. Keep cancellation as a hypothesis, explicitly flagged as unsupported by the
   apo data at three replicates.
3. Point instead at the G190E-specific confounds already in [112] — the
   replication-capacity problem and the single-isolate phenotype — and treat
   the ΔΔGbind value as untestable for this genotype.

The same sentence is quoted in the Reviewer 1 comment 3 and Reviewer 2 comment
2 responses, so whatever you choose has to propagate to both.

---

## Part 2 — the rearrangement

A note on how to apply it: **these paragraphs hold 26 live Zotero fields**
(five in [107] alone), and moving them with a script risks orphaning
citations. Cut-and-paste inside Word preserves fields. The list below is
therefore ordered so you can do it top to bottom in a few minutes.

### Opening — unchanged
1. **[106]** context. 2. **[107]** approach. 3. **[108]** headline result.

### Mechanisms — one per paragraph, limitations removed
4. **[109]** K103N and the backbone anchor. Unchanged.
5. **[111]** Tyr188 carries the aromatic packing. Unchanged.
6. **[112]** the G190 series — **delete the last two sentences** (replication
   capacity; "might affect the primer grip, which is not modeled here"). They
   move to §11 and §12.
7. **[113]** the V106 series — **delete the last 90 words**, from "We
   quantified the uncertainty…" to "…which they do not", and replace with:
   > The displacement is well resolved (−1.35 Å [−1.73, −0.97] for V106A
   > relative to wild type, Hedges' g = −9.5); the change in pyridinone contact
   > count is not, and the argument rests on the displacement.
8. **[114]** ΔΔETotal against ΔΔGbind — **delete the final sentence** on
   implicit solvent, which moves to §11, and settle Part 1 above.
9. **[115]** K103N+M230L — **delete the last two sentences** on the binary
   complex, which move to §11.

### Limitations — one per paragraph, each stated once
10. **Sampling.** Three replicates from one seeded pose sample one basin;
    convergence within a basin is not evidence the basin is the only relevant
    one; the unliganded equilibrium is slower than 100 ns. Add: of 1,287
    genotype-by-feature structural comparisons, 66 reach 95%, which is the
    null expectation — most structural differences in the panel are
    unresolved at three replicates.
11. **A simplified model system.** The binary RT–DOR complex omits
    template/primer and dNTP, so primer-grip and assembly mechanisms are
    structurally absent (merged from [112] and [115], currently stated twice);
    MM/GBSA uses implicit solvent, which may overestimate desolvation for a
    buried glutamate (from [114]); the ligand charge model shifts ΔΔETotal by
    up to 1.33 kcal/mol without changing its sign.
12. **Mechanisms we cannot reach.** [116] as written, plus the
    replication-capacity sentence moved out of [112].

### What leaves the Discussion
- The 90-word statistical passage in [113] → Results or SI.
- "comparable distortion is present in the apo simulations" → see Part 1.
- [112]'s "pocket enlarged by ~25%" reads as a new result here; it is in
  Table 3 and can be cited rather than restated.

---

## Verified while doing this

All 17 energetic and structural numbers quoted in the Discussion check out
against Tables 2 and 3: every ΔΔGbind, ΔΔETotal, ΔΔEvdW and ΔΔEGB value, the
chlorocyanophenyl burial counts for G190E, Y188L and wild type, and the Lys103
anchor distance. Only the apo-derived claims needed attention.
