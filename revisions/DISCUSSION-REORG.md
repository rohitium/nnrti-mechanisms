# Discussion, reorganised

Mechanisms separated from limitations, as you asked. This is a **rearrangement
of your existing text**, not new prose: each block below names the paragraph it
comes from and what was taken out of it. New connective sentences are marked
`[NEW]` and there are only four of them.

Current: 11 paragraphs, 1,436 words, with limitations embedded in five of them.
Proposed: 12 paragraphs, roughly the same length, limitations in three.

## The problem being fixed

The same limitation is currently stated more than once because it is attached
to whichever mechanism raised it:

| limitation | currently appears in |
|---|---|
| no nucleic acid / primer grip | **[112]** (G190E "might affect the primer grip, which is not modeled here") and **[115]** (K103N+M230L "no template/primer… absent from our simulations") |
| sampling / timescales | **[116]**, and implicitly in [114] via "comparable distortion is present in the apo simulations" |
| implicit solvent | **[114]** only, buried at the end of a mechanism paragraph |
| three replicates not resolving a metric | **[113]**, as a 90-word statistical aside inside the V106 mechanism |

## Proposed order

### Opening — unchanged
1. **[106]** DOR's distinct resistance profile; sparse public phenotype data.
2. **[107]** What we did and why; prior free-energy work in HIV and cancer.
3. **[108]** Headline result: 3 of 4 Susceptible, 7 of 9 Resistant.

### Mechanisms — one per paragraph, limitations removed
4. **K103N and the backbone anchor** — [109] as written. Already clean.
5. **Tyr188 carries the aromatic packing alone** — [111] as written. Already
   clean, and the strongest literature-anchored paragraph in the section.
6. **The G190 series** — [112] **minus its final two sentences** (the
   replication-capacity point and "might affect the primer grip, which is not
   modeled here"). Those move to §11 and §12.
7. **The V106 series** — [113] **minus the last 90 words** (from "We quantified
   the uncertainty…" to "…which they do not"). The interval itself is worth
   keeping in the mechanism; the methodological hedging is not. Replace with
   one clause:
   > …loses roughly a quarter of its pyridinone contacts. The displacement is
   > well resolved (−1.35 Å [−1.73, −0.97] for V106A relative to wild type);
   > the change in contact count is not, and the argument rests on the
   > displacement. `[NEW — compressed from the deleted block]`
8. **Where ΔΔETotal and ΔΔGbind disagree** — [114] **minus its final sentence**
   on implicit solvent, which moves to §11.
9. **K103N+M230L** — [115] **minus its last two sentences** on the binary
   complex, which move to §11. What remains is the mechanism: Met230 sits in
   the primer grip, M230L impairs replicative capacity, and the DOR interface
   is unperturbed.

### Limitations — one per paragraph, each stated once
10. **Sampling.** `[NEW — one opening sentence]` Then the sampling material
    from [116], plus the apo-equilibrium point. Covers: three replicates seeded
    from one pose sample one basin; convergence within a basin is not evidence
    that the basin is the only relevant one; the unliganded pocket equilibrium
    is slower than we sample; and which structural features are unresolved at
    three replicates.
11. **A simplified model system.** `[NEW — one opening sentence]` Then, each
    stated once: the binary RT–DOR complex omits template/primer and dNTP,
    so primer-grip and assembly mechanisms are structurally absent (from [112]
    and [115], now merged); MM/GBSA uses implicit solvent, which may
    overestimate the desolvation penalty for a buried glutamate (from [114]);
    the ligand charge model shifts ΔΔETotal by up to 1.33 kcal/mol without
    changing its sign.
12. **Mechanisms we cannot reach.** [116] as written, which already does this
    job — allosteric change, unliganded equilibrium, processivity and fitness —
    plus the replication-capacity sentence moved out of [112], and the closing
    pointer to QM/MM, enhanced sampling and machine learning.

## What this costs

Three things leave the Discussion and should land elsewhere if you want to keep
them:

- The 90-word statistical passage in [113]. It belongs in Results or the
  Supporting Information, not in interpretation.
- "comparable distortion is present in the apo simulations" in [114] is a
  **result**, not an interpretation, and its number is stale — see the numbers
  memo. It needs a Results sentence to point at.
- [112]'s "pocket enlarged by ~25%" is also a result. The current figure is
  from the holo pocket volumes and is fine, but it reads as new data here.

## Numbers in the Discussion that need checking

These are the ones I have not yet been able to confirm against regenerated
tables — the apo analysis is still rerunning after a stale PBC sidecar was
found and rebuilt:

- [114] "comparable distortion is present in the apo simulations" — the apo
  numbers behind this all predate the 100 ns runs.
- [112] "pocket enlarged by ~25%" for G190E.

Everything else in the Discussion cites ΔΔGbind, ΔΔETotal or Table 3 values,
which are unaffected by the apo rerun.
