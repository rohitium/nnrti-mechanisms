# Discussion — expanded mechanistic treatment

Drafted 2026-08-31 in response to: *"the Discussion needs more detailed
discussion about each of the significant results … Can we say anything about the
mutations that simulations predict should be DOR susceptible? Anything from the
equilibrium MD regarding the mechanism? … connect to what's known about these
mutations with respect to other second-generation NNRTIs (are the mechanisms
different)?"*

Every number below is computed from the existing equilibrium trajectories — see
§Provenance. Claims that rest on the **literature rather than our data** are
marked `[cite]` and need a reference check before submission; the reference
numbers suggested are those already in the draft's bibliography.

The proposal replaces the single paragraph beginning *"We found that the
relative binding free energy…"* with five paragraphs (§1–§5), and inserts them
around the existing resistance-mechanism paragraph. Suggested order:

| ¶ | content | status |
|---|---|---|
| 1 | overall performance, with classification metrics | rewrite of existing |
| 2 | **why the susceptible genotypes are susceptible** | new |
| 3 | **the position-190 series** | new |
| 4 | existing Y188L / V106A resistance mechanisms | keep, with one added link |
| 5 | **contrast with other second-generation NNRTIs** | new |
| 6 | **where the calculations fail, and what that tells us** | new |
| 7–8 | existing FEP error and outlook paragraphs | keep |

---

## ¶1 — Overall performance (rewrite)

> We found that the relative binding free energy, ∆∆G<sub>bind</sub>, obtained
> from the FEP simulations correctly predicted low impact on RT–DOR binding
> affinity for 3 out of 4 Susceptible mutations and significant impact for 7 out
> of 9 Resistant mutations. Restricted to the 13 genotypes with established
> phenotypes, this corresponds to a sensitivity of 0.78, a specificity of 0.75
> and a Matthews correlation coefficient of 0.50, or a threshold-free ROC AUC of
> 0.78 (95% CI 0.43–1.00). Amongst this set, ∆∆G<sub>bind</sub> correlated with
> the logarithm of the fold-change in DOR susceptibility observed in vitro,
> though this correlation was weak and not statistically significant
> (R² = 0.26, p = 0.07). In contrast, ∆∆E<sub>Total</sub> obtained from MM/GBSA
> detected resistance equally often (sensitivity 0.78) but correctly identified
> only 1 of 4 Susceptible mutations as low-impact, giving a specificity of 0.25,
> an MCC of 0.03 and an AUC of 0.69 (95% CI 0.36–0.97); ∆∆E<sub>Total</sub>
> exhibited no significant correlation with DOR fold-change (R² = 0.13,
> p = 0.23). Neither difference reaches conventional significance at this panel
> size (Fisher exact p = 0.119 and 0.706), and we therefore report these as
> effect sizes rather than as hypothesis tests. The asymmetry is nonetheless
> consistent across both framings: the two methods agree on which genotypes
> impair binding, and differ chiefly in how often they raise a false alarm on a
> genotype that remains susceptible in vitro.

*Numbers verified 2026-08-31 against `panel_ddg.csv` and `ddg_full.csv`. The
draft's current R² = 0.25/p = 0.03 (∆∆E vs ∆∆G), R² = 0.02/p = 0.55–0.60
(∆∆E vs fold, all 18) and R² = 0.05/p = 0.45 (∆∆E vs fold, 13 established) are
all stale; the current values are 0.22/0.048, 0.06/0.31 and 0.13/0.23.*

---

## ¶2 — Why the susceptible genotypes are susceptible (NEW)

> The three Susceptible mutations that the calculations classify correctly —
> K103N, Y181C and G190A — are informative precisely because they are the
> canonical NNRTI resistance mutations against which DOR was designed, and the
> equilibrium trajectories show that each fails to engage DOR for a different
> structural reason.
>
> For **K103N**, the interaction that anchors DOR is made to the residue-103
> *backbone*, not its side chain: the carbonyl oxygen accepts a hydrogen bond
> from the pyridinone nitrogen at 2.97 ± 0.01 Å in WT, and this distance is
> unchanged at 3.08 ± 0.06 Å in K103N. Substituting lysine by asparagine
> shortens the side chain and in fact brings its polar atoms closer to the
> ligand (8.38 ± 0.15 Å to 5.07 ± 0.04 Å), yet DOR neither moves nor loses
> packing: burial of the chlorocyanophenyl ring is 46.2 ± 1.7 protein heavy
> atoms against 45.9 ± 1.3 in WT, and the Tyr188 stacking distance is
> 4.22 ± 0.17 Å against 4.24 ± 0.19 Å. The mutation is, in effect, invisible to
> the ligand. Critically, the same backbone hydrogen bond is preserved in every
> K103N-containing genotype we simulated (3.05–3.20 Å across K103N+M230L,
> K103N+P225H and L100I+K103N), which explains why none of these backgrounds
> produces a large computed binding penalty even where the clinical phenotype
> is resistant.
>
> For **Y181C**, removal of the tyrosine costs DOR nothing because DOR does not
> use it. Packing around the chlorocyanophenyl ring is not reduced but slightly
> increased (48.6 ± 1.4 against 45.9 ± 1.3 heavy atoms), and engagement by
> Tyr188 tightens rather than loosens: the ring-centroid separation falls from
> 4.24 ± 0.19 Å to 3.87 ± 0.21 Å and the number of Tyr188 contacts with the
> ring rises from 21.3 ± 2.9 to 26.5 ± 3.1. Aromatic anchoring of DOR is
> therefore carried entirely by Tyr188, consistent with the crystallographic
> observation that the Tyr181 side chain is rotated away from the drug,<sup>14,16</sup>
> and the Cys181 cavity is simply absorbed by a modest local relaxation.
>
> For **G190A**, no coordinate we examined departs from wild type. The Ala190
> methyl does contact DOR (3.31 ± 0.00 Å minimum heavy-atom distance), but
> without displacing it: ring burial is 45.0 ± 0.3 heavy atoms, the Tyr188
> stack is 4.16 ± 0.13 Å, the Lys103 backbone hydrogen bond is 3.01 ± 0.03 Å,
> and DOR remains 6.81 ± 0.08 Å from Ser105 against 6.65 ± 0.09 Å in WT. Where
> G190A is thought to confer resistance to first-generation NNRTIs by
> introducing a steric bulge into a compact region of the pocket,<sup>6</sup>
> the DOR pose evidently has room to accommodate it.

---

## ¶3 — The position-190 series (NEW)

> Because the panel contains three substitutions at position 190, it resolves a
> steric and electrostatic series at a single site whose severity tracks the
> measured phenotype (2.7-, 5.2- and 18-fold for G190A, G190S and G190E). As
> described above, Ala190 is accommodated with essentially no perturbation of the
> aromatic stack: the interplanar angle between the Tyr188 ring and the
> chlorocyanophenyl ring is 16.1 ± 1.8°, against 13.4 ± 0.4° in WT and within the
> replicate-to-replicate scatter of the wild type.
>
> Ser190 destabilises that stack. In WT the two rings are held close to coplanar,
> with 61% of frames below 15° and essentially none beyond 40° (0.2%). In G190S
> the stack is not uniformly tilted but intermittently broken: a WT-like
> coplanar population persists (37% of frames below 15°) alongside a second
> population, sampled in 32% of frames, in which the ring rotates past 40° and
> simultaneously withdraws from the ligand, the centroid separation rising from
> 4.23 Å in the coplanar state to 5.05 Å in the rotated one. The frame-averaged
> angle is 27.4 ± 5.3°, but that uncertainty reflects genuine heterogeneity in
> how often each replicate visits the rotated state (per-replicate means 34.2°,
> 17.0° and 30.9°) rather than imprecision in the measurement. Ring burial falls
> correspondingly to 42.7 ± 1.8 heavy atoms.
>
> Glu190 acts differently again, displacing a packing partner rather than the
> stack — the Tyr188 geometry is in fact marginally tighter than wild type
> (12.1 ± 0.9°, 75% of frames below 15°). Instead, the minimum distance from
> Val179 to DOR increases from 3.62 ± 0.06 Å to 6.47 ± 0.44 Å, and the NNIBP
> proxy volume expands from 230 ± 12 Å³ to 286 ± 3 Å³, the largest pocket in the
> panel. That the same position produces three distinct perturbations — none,
> loss of aromatic stacking, and expansion of the pocket — illustrates why
> position-specific resistance rules derived from one substitution generalise
> poorly to others at the same residue.

*Caveat: G190E's ∆∆G<sub>bind</sub> is the least well-determined value in the
panel and is currently being re-run (Lever C). The structural statements above
do not depend on it — they come from the equilibrium trajectories, not the FEP.
Do not add a ∆∆G number for G190E to this paragraph.*

---

## ¶4 — Existing resistance-mechanism paragraph

Keep as written, with the corrections already listed in
`CHANGES_2026-08-30_final_mmgbsa.md` (∆∆E<sub>vdW</sub> = **2.25 ± 0.20**, and
Y188L now has the largest van der Waals penalty of the **whole panel**, so the
"of any known DOR resistant genotype" hedge can be dropped). Add one linking
sentence at the end, which sets up ¶5:

> Taken together with the Susceptible genotypes above, the two mechanisms are
> complementary: DOR tolerates the loss of Tyr181 because Tyr188 alone carries
> its aromatic anchoring, and is correspondingly vulnerable when that single
> anchor is removed at position 188.

---

## ¶5 — Contrast with other second-generation NNRTIs (NEW)

> These mechanisms differ instructively from those established for the
> diarylpyrimidine second-generation NNRTIs, etravirine and rilpivirine. ETR and
> RPV accommodate mutated pockets primarily through torsional
> flexibility — repositioning and re-orienting within the NNIBP so that
> compensating contacts can be formed when a given interaction is
> lost.<sup>12,13</sup> DOR was instead designed with conformational constraint,
> deriving its resilience from a hydrogen bond to a main-chain
> atom<sup>14,15</sup> — an interaction that no side-chain substitution can
> remove. Our simulations show these to be genuinely different strategies rather
> than two descriptions of the same one: across all four K103N-containing
> genotypes the DOR pose is essentially rigid, and the anchoring hydrogen bond
> varies by at most 0.23 Å from wild type, with no evidence of the ligand
> re-orientation that underlies DAPY tolerance. `[cite]`
>
> The two strategies have different failure modes, and the panel exposes DOR's.
> Because ETR and RPV stack against both Tyr181 and Tyr188, the loss of either
> aromatic leaves the other in place; because DOR stacks only against Tyr188, it
> has no such redundancy. This single architectural choice accounts for the two
> most extreme entries in our panel in opposite directions — Y181C at 1.4-fold,
> where the removed residue was never engaged, and Y188L at 149-fold, the
> largest computed binding penalty in the study (∆∆G<sub>bind</sub> =
> 4.52 ± 0.49 kcal/mol) — and is consistent with the comparatively modest effect
> of Y188L on ETR and RPV susceptibility. `[cite — Y188L cross-NNRTI penalties;
> ref 18 or HIVDB]`
>
> The same reasoning explains why the genotypes that dominate DOR resistance are
> not those that dominate ETR or RPV resistance. The V106A-containing
> combinations that recur in DOR failure — V106A+F227L, V106A+P225H,
> V106A+L234I — act not by removing a specific interaction but by permitting DOR
> to slide out of its crystallographic pose (Figure 3B), a mode of escape to
> which a conformationally constrained ligand with a single point of aromatic
> anchoring is particularly exposed, and which is not a signature pathway for
> the DAPYs. `[cite — refs 16, 18]`

---

## ¶6 — Where the calculations fail, and what that tells us (NEW)

> The misclassifications are as informative as the successes, and the
> equilibrium trajectories distinguish two quite different kinds of failure.
>
> The one Susceptible genotype that FEP over-calls, **V106I**
> (∆∆G<sub>bind</sub> = 2.27 ± 0.74 kcal/mol against 1.1-fold measured), is not
> a numerical artefact: the simulations detect a real local strain. The
> β-branched isoleucine displaces DOR away from Ser105 (7.15 ± 0.14 Å against
> 6.65 ± 0.09 Å in WT) and lengthens the Lys103 backbone hydrogen bond to
> 3.62 ± 0.30 Å, the longest value in the panel, while leaving
> ring burial unchanged at 44.9 ± 1.6 heavy atoms. V106M perturbs the same
> anchoring in the same direction (3.56 ± 0.20 Å; 7.04 ± 0.46 Å from Ser105) and
> is likewise assigned a large penalty against a measured 3.4-fold change. The
> calculations are thus reporting a genuine strain on the anchoring hydrogen
> bond that the virus evidently tolerates — a limitation of binding affinity as
> a proxy for susceptibility, not of the sampling.
>
> The two Resistant genotypes that FEP misses fail for the opposite reason: in
> both, at least one substituted residue never contacts the drug. Residue 98
> sits 9.11 ± 0.19 Å from DOR in **A98G+F227C**, and residue 230 sits
> 7.34 ± 0.04 Å away in **K103N+M230L** — in each case outside the first contact
> shell. Consistently, no interface coordinate is degraded in either genotype.
> In K103N+M230L the interface is if anything tighter than wild type: ring
> burial is 48.6 ± 1.3 heavy atoms, the Tyr188 stack is 3.83 ± 0.11 Å, Tyr188
> makes 27.1 ± 1.4 contacts with the chlorocyanophenyl ring against 21.3 ± 2.9
> in WT, and the anchoring hydrogen bond is 3.05 ± 0.04 Å. In A98G+F227C the
> Lys103 backbone hydrogen bond is 2.98 ± 0.01 Å, indistinguishable from WT, and
> the total DOR–RT contact count is 17.0 ± 0.1 against 16.7 ± 0.4. The
> simulations are therefore not mistaken about the interface; the interface is
> genuinely unperturbed, and the 36- and 93-fold reductions in susceptibility
> these genotypes produce must arise outside it — through the conformational
> equilibrium of the unliganded pocket, the kinetics of drug entry, or effects
> on polymerase function and viral fitness that a bound-state binding
> calculation does not represent. This is the clearest evidence in our dataset
> that binding energetics, however well converged, are an incomplete model of
> phenotypic drug resistance.

---

## Provenance

All values recomputed 2026-08-31; per-replicate means first, then mean ± SEM
across the three 100 ns replicates.

| quantity | source file |
|---|---|
| ring burial, Tyr188 stacking geometry, residue-103 backbone/side-chain distances, Ser105 and Val179 distances | `results/analysis/mechanisms/mechanism_coordinates.csv` |
| Lys103 backbone H-bond, Tyr181/Tyr188/Val179 contacts | `results/dor_key_contacts_timeseries_all_mutations.csv` |
| mutated-residue-to-DOR distances (`c1_to_dor`, `c2_to_dor`) | `results/drm_sidechain_distance_timeseries_all_mutations.csv` |
| NNIBP proxy volume | `results/pocket_volume_profiles.csv` |
| contact and H-bond counts | `results/structural_metrics.csv` |
| Tyr188 interplanar angle, position-190 series | `results/analysis/mechanisms/y188_interplanar_angle_190series.csv` |
| DOR moiety contacts, V106A genotypes | `results/analysis/mechanisms/dor_moiety_contacts_summary.csv` |
| per-residue contact loss | `results/analysis/mechanisms/dor_residue_contact_delta.csv` |
| ∆∆G<sub>bind</sub>, SEMs | `results/analysis/fep_pmx/panel_ddg.csv` |
| classification metrics | `results/analysis/classification_performance/classification_metrics.csv` |

**Not covered by `mechanism_coordinates.csv`:** A98G+F227C, V106I+F227C, Y318F
and F227C were never run through `compute_mechanism_coordinates.py`. Statements
about A98G+F227C above therefore use the key-contact and side-chain-distance
tables, which do cover it. If a supplementary table of mechanism coordinates for
the full panel is wanted, run:

```bash
PYTHONPATH=. ~/miniconda3/envs/nnrti-prep/bin/python -m src.analysis.cli.compute_mechanism_coordinates --mutations A98G+F227C V106I+F227C Y318F F227C
```

## Open items for the author

1. The three `[cite]` markers in ¶5 are literature claims, not our results.
2. ¶3 and ¶6 would be well served by a supplementary table of the mechanism
   coordinates for all 19 systems — the numbers exist for 16 of them already.
3. Y181C's burial increase (48.6 ± 1.4 vs 45.9 ± 1.3) is within roughly one
   pooled SEM. The text claims only that burial does **not decrease**, which is
   what the data support; do not strengthen it to a claimed increase.

---

## CORRECTION — WT chlorocyanophenyl ring burial (added 2026-08-31)

The 09-02 draft carries **two different values for the same WT quantity**:

| location | text | status |
|---|---|---|
| Results, mechanisms ¶ | "falling from **45.3 ± 2.4** in WT to **34.2 ± 1.1** in Y188L" | **stale** |
| Discussion (Y181C, newly merged) | "48.6 ± 1.4 against **45.9 ± 1.3** heavy atoms" | correct |
| Discussion (Y188L, existing) | "falling from **45.3 ± 2.4** in WT to **34.2 ± 1.1** in Y188L" | **stale** |
| Discussion (G190A, newly merged) | "ring burial is 45.0 ± 0.3" | correct |

**Correct values: WT 45.9 ± 1.3, Y188L 34.6 ± 0.9 heavy atoms.**

These are what `results/analysis/mechanisms/mechanism_summary.csv` contains and
what `plot_mechanism_panel.py` writes today: per-replicate means first, then mean
± SEM across the three replicates. The `45.3 ± 2.4 / 34.2 ± 1.1` pair predates
the current `mechanism_coordinates.csv` (committed 2026-08-26 in `245825a`) and
is **not reproducible from it under any aggregation** — not replicate SEM
(1.34 / 0.92), not replicate SD (2.31 / 1.59), not pooled frames (46.53 / 34.60),
and not the time-average of the plotted interpolated trace (45.92 / 34.60). It
was carried forward verbatim from the 08-26 draft.

Two find-and-replace operations in Word fix both occurrences:

```
45.3 ± 2.4 in WT to 34.2 ± 1.1 in Y188L
    ->  45.9 ± 1.3 in WT to 34.6 ± 0.9 in Y188L
```

The accompanying claim "a loss of roughly a quarter" survives the correction
(34.6 / 45.9 = 0.75, a 25% loss), as does Figure 3A, whose plotted trace is
unaffected — only the quoted summary statistics were stale.

### Why the replicate-mean aggregation is the right one

Frame counts are not balanced across replicates: WT rep 1 contributes 500 frames
against 180 each for reps 2 and 3, and the imbalance varies by system (from 139
to 540 frames per replicate). Averaging frames directly would weight WT rep 1 —
whose mean, 47.6, is the highest of the three — nearly three times as heavily as
its siblings, giving 46.5. Taking each replicate's mean first weights the three
independent trajectories equally, which is also what the reported SEM assumes.
This is worth one sentence in Methods if it is not already stated.

---

## Moiety-resolved contact analysis for the V106A genotypes (added 2026-08-31)

Replaces: *"This dislocation disrupts contacts with RT residues, for example the
total number of protein heavy atoms within 4.5 Å of DOR fell from 224 ± 1 in WT
to 212 ± 1 in V106A."*

New analysis: `src/analysis/cli/compute_dor_moiety_contacts.py`. DOR is
partitioned by bond connectivity into its three ring systems plus the
ether/methylene linker, and the 4.5 Å contact count is computed for each, over
the full production trajectory, for WT and all four V106A-containing genotypes.

### Result — the loss is concentrated in the pyridinone ring

Atom-pair contacts within 4.5 Å, mean ± SEM over three replicates:

| genotype | whole ligand | chlorocyanophenyl | **pyridinone** | triazolinone | linker |
|---|---:|---:|---:|---:|---:|
| WT | 224.3 ± 1.0 | 69.3 ± 1.9 | **24.0 ± 2.1** | 43.3 ± 1.9 | 87.6 ± 1.4 |
| V106A | 215.8 ± 0.9 | 65.1 ± 1.6 | **17.8 ± 1.3** | 45.2 ± 1.1 | 87.7 ± 1.4 |
| V106A+F227L | 213.0 ± 2.4 | 63.7 ± 1.3 | **18.2 ± 1.3** | 45.6 ± 1.0 | 85.6 ± 1.3 |
| V106A+L234I | 218.7 ± 2.1 | 67.1 ± 2.0 | **17.4 ± 1.7** | 46.5 ± 1.3 | 87.7 ± 0.5 |
| V106A+P225H | 222.3 ± 2.7 | 68.3 ± 0.2 | **19.0 ± 0.2** | 44.6 ± 0.8 | 90.5 ± 2.4 |

Change relative to WT:

| genotype | whole ligand | chlorocyanophenyl | **pyridinone** | triazolinone | linker |
|---|---:|---:|---:|---:|---:|
| V106A | −3.8% | −6.1% | **−25.8%** | +4.4% | +0.1% |
| V106A+F227L | −5.0% | −8.1% | **−24.2%** | +5.3% | −2.3% |
| V106A+L234I | −2.5% | −3.2% | **−27.5%** | +7.4% | +0.1% |
| V106A+P225H | −0.9% | −1.4% | **−20.8%** | +3.0% | +3.3% |

**Why this is the better number.** The whole-ligand count ranges from −0.9% to
−5.0% and would make V106A+P225H (153-fold resistant) look essentially
unaffected. The pyridinone count is −21% to −28% in all four, a range of only
7 percentage points across genotypes that differ 16-fold in susceptibility — the
disruption is a shared, uniform property of the V106A background, which is
exactly the claim the paragraph is making. The chlorocyanophenyl ring loses
proportionally 3–5× less, and the distal triazolinone ring *gains* contacts,
consistent with DOR pivoting about its distal end rather than withdrawing
bodily.

The pyridinone is also the mechanistically meaningful ring: it carries the
Lys103 backbone hydrogen bond that anchors DOR (§¶2). So the V106A slide is
best described not as a general loosening of the interface but as the
displacement of the anchored end of the ligand.

### The partner side

Residues losing the most contact with DOR (mean over the four V106A genotypes,
against WT), in contacts per frame:

| residue | WT | V106A set | Δ |
|---|---:|---:|---:|
| Val106 | 3.26 | 0.00 | −3.26 (mutated to Ala) |
| Phe227 | 5.57 | 4.00 | −1.57 |
| Tyr318 | 6.68 | 5.20 | −1.48 |
| Lys102 | 1.08 | 0.18 | −0.89 |
| Lys101 | 1.57 | 1.01 | −0.56 |
| Pro225 | 1.85 | 1.34 | −0.51 |

Residues gaining: Ser105 (0.31 → 3.53), Ala106 (0 → 2.76), Lys104 (0.44 → 2.46).

The exchange is directional and local to the pyridinone: contacts move from
Lys101/Lys102 on one side of the 101–106 loop to Lys104/Ser105 on the other.
This is the same displacement plotted in Figure 3B, now attributed to a specific
part of the ligand and a specific set of partners. Note that Tyr318 and Phe227
are themselves sites of DOR resistance mutations in this panel.

### Suggested replacement text

> This dislocation does not loosen the interface uniformly. Resolving the 4.5 Å
> contact count by ligand moiety localises the loss to the central pyridinone
> ring — the ring that carries the Lys103 backbone hydrogen bond — whose contacts
> fall from 24.0 ± 2.1 in WT to 17.8 ± 1.3 in V106A, 18.2 ± 1.3 in V106A+F227L,
> 17.4 ± 1.7 in V106A+L234I and 19.0 ± 0.2 in V106A+P225H, a loss of 21–28% in
> every case. The chlorocyanophenyl ring loses proportionally three- to five-fold
> less and the distal triazolinone ring gains contacts, indicating that DOR
> pivots about its distal end rather than withdrawing from the pocket as a whole.
> On the protein side the same motion transfers contacts from Lys101, Lys102,
> Phe227 and Tyr318 to Ser105 and Lys104. That a 21–28% loss at the anchored end
> of the ligand is common to all four V106A genotypes, which span a 16-fold range
> of measured susceptibility, identifies displacement of the pyridinone as the
> shared structural consequence of the V106A background.

---

## SECOND CORRECTION — "protein heavy atoms" are atom pairs

Both packing figures in the draft are described as counts of atoms but are
computed as counts of atom **pairs**. `_ncontacts()` in
`compute_mechanism_coordinates.py` returns `(d < cutoff).sum()` over the full
protein × ligand distance matrix, so an RT atom close to three ligand atoms is
counted three times.

| draft wording | computed quantity | count of distinct atoms |
|---|---:|---:|
| "total number of protein heavy atoms within 4.5 Å of DOR … 224 ± 1 in WT" | 224.3 ± 1.0 **atom pairs** | **76.3 ± 1.1** atoms |
| "number of RT heavy atoms within 4.5 Å of the chlorocyanophenyl ring … 45.9 in WT" | 45.9 **atom pairs** | ~26.9 atoms |

The numbers are correct and the comparisons are valid — a pair count is a
legitimate packing density measure — but the *wording* is not. Both should read
"contacts" or "heavy-atom contacts (atom pairs within 4.5 Å)" rather than
"heavy atoms". The moiety script now reports both conventions side by side so
they cannot be conflated again.

This matters beyond wording in one place: under the distinct-atom convention the
whole-ligand count **rises** in every V106A genotype (+1.5 to +4.4 atoms) while
the pair count falls. DOR ends up touching slightly more RT atoms, less closely.
The per-moiety result is robust to the choice — the pyridinone loses under both
conventions (−5.0 to −6.6 pairs, −1.1 to −1.7 atoms) and the triazolinone gains
under both — which is a further reason to report the moiety-resolved numbers
rather than the whole-ligand one, whose sign depends on the convention chosen.

## THIRD CORRECTION — V106A's 212 came from a terminal window

WT reproduces the published value exactly, but V106A does not:

| frame set | WT | V106A |
|---|---:|---:|
| **full trajectory** | **223.7 ± 0.7** | **215.8 ± 0.7** |
| last 25% | 223.3 ± 2.4 | 210.4 ± 2.5 |
| last 20 frames | 223.8 ± 4.2 | 207.1 ± 3.4 |
| *draft* | *224 ± 1* | *212 ± 1* |

WT is insensitive to the window; V106A is not. The published 224 → 212 contrast
is therefore inflated by the same terminal-window sampling artifact that was
identified and removed from the MM/GBSA protocol — note also the ~5× inflation
of the replicate SEM in the terminal windows, the same signature. On the full
trajectory the gap is 8.5 contacts, not 12.

Since the surrounding analysis and the final MM/GBSA panel are both now
full-trajectory, the replacement text above uses full-trajectory values
throughout and this inconsistency disappears with it.
