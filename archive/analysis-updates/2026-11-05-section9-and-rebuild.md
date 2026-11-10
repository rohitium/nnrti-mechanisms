# Section 9 rewrite, and the Table 4 / Figure 5 rebuild behind it

2026-11-05. Addresses Reviewer 1 comments 2 and 3.

## What changed in the submission package

| File | Change |
|---|---|
| `Supporting-Information-Clean.docx` | Section 9 rewritten (1 paragraph -> 5); Supplementary Figure 6 replaced; its caption rewritten |
| `Supporting-Information-Clean.pdf` | **STALE — must be re-exported from the DOCX.** No LibreOffice on this machine |
| `Supplementary-Table-4.xlsx` | Rebuilt. 68 cells changed, all in the unliganded state |
| `Manuscript-Clean.docx`, `Manuscript-Marked.docx` | Figure 5 art replaced; four numbers in the Tyr181/Tyr188 paragraph updated |

Backup of every published table this touched:
`revisions/analysis-updates/published-tables-backup-20261105/`.

## Why the rebuild was needed

Two upstream tables had fallen behind the trajectories on disk.

1. **`chi1_per_replicate.csv` averaged as few as 7 frames** for some apo
   replicates, against the 200-380 now present; 18 of 114 apo entries used
   fewer than 50 frames. Every holo entry was correct.
2. **V106A's apo rows predated the V106A regeneration** (commit `abe566c`).
   Replicates 2 and 3 disagreed with the current trajectories by up to 0.32 A
   in fingers RMSF; replicate 1 matched.

Both are now regenerated. The apo chi1 minimum frame count is 200.

## What the rebuild did and did not change

Figure 5 carries **11 resolved cells before and 11 after, and not one cell
changed status.** Every claim in the manuscript survives. The movement is
confined to cells that were never significant, the largest being K103N Tyr188
(+7 -> +57 deg) and K103N+M230L Tyr181 (+8 -> +46 deg). V106A pocket volume
moved from -30 +/- 18 to -7 +/- 29 A^3, also unresolved either way.

Four numbers in the manuscript needed updating, all in the Tyr181/Tyr188
paragraph of the Results:

| Was | Now |
|---|---|
| G190S unliganded Tyr181 `-52 deg` | `-53 deg` |
| G190S induced fit `-83 +/- 28 deg ([-148, -18])` | `-82 +/- 28 deg ([-147, -17])` |
| A98G+F227C `6.9 standard errors` | `7.0 standard errors` |
| Y318F interval `[22, 169]` | `[21, 170]` |

Every RMSF value quoted in the Results is unchanged, including the
K103N+M230L primer-grip numbers and the A98G+F227C pocket-lining numbers.

## The new Supplementary Figure 6

`src/nnrti/cli/plot_rt_metric_violins.py`. Ten panels, one per Figure 5
metric, nineteen clinical genotypes, bound and unliganded side by side.

The unit of observation differs by metric, and each panel says which:

| Metric | Unit | Reproduces Supplementary Table 4? |
|---|---|---|
| NNIBP volume | per frame | apo exact; holo within 0.42 of the frame-to-frame SD |
| Tyr181 / Tyr188 chi1 | per frame | all 60 holo replicates exact |
| 6 Ca RMSF regions, b9-b10 strand | per 10 ns block | 118 of 120 replicates exact |

An RMSF has no per-frame value, so those panels are blocked; a block RMSF is
smaller than the 100 ns RMSF of the same replicate, and the published value is
overlaid as an open marker so the two are not confused.

chi1 keeps the signed (-180, +180] convention of Supplementary Table 4, but
each chi1 panel's axis is rotated to wrap at its least populated angle. A
straight violin on the signed axis draws a trans rotamer as two half violins
pinned to opposite edges with the kernel smeared through the middle, which
asserts a bimodality the data do not contain.

## Result

Of 178 genotype-by-metric comparisons, **167 are concordant** between the two
states: the mutation moves the unliganded and bound ensembles by the same
amount within error, which is what the thermodynamic cycle cancels. The eleven
exceptions are the resolved cells of Figure 5.

G190E has the largest bound pocket in the panel, 56 A^3 above WT, but its
unliganded pocket is 15 A^3 *smaller* than WT's: the expansion is induced fit,
not a property of the resting enzyme. K103N+M230L is the opposite, concordant
on six of seven fluctuation metrics and discordant only in the primer grip.

For comment 2 the answer is negative and stated as such: the P225- and
M230-containing genotypes show no distinct unliganded state. They do have the
widest block-to-block scatter in the panel, but it does not survive correction
for the eighteen comparisons made (smallest Holm-adjusted p = 0.84). Three
100 ns replicates constrain where the unliganded ensemble sits, not how broad
it is.

## New and changed code

- `src/nnrti/cli/plot_rt_metric_violins.py` — the figure
- `src/nnrti/cli/compute_block_metrics.py` — 10 ns block RMSF and strand content
- `src/nnrti/cli/compute_chi1_frame_series.py` — per-frame chi1, both states
- `src/nnrti/cli/compute_holo_volume_frames.py` — per-frame bound pocket volume
- `src/nnrti/cli/diff_structural_rebuild.py` — what a rebuild changes, cell by cell
- `src/nnrti/cli/compute_pocket_volume.py` — also writes a per-frame companion
- `ops/docx/replace_figure.py` — swap art above or inside a caption paragraph
- `ops/docx/insert_paragraphs.py` — grow a section past one paragraph
- `ops/docx/edit_docx_text.py`, `ops/docx/insert_paragraphs.py` — preserve the
  document's namespace prefixes across an ElementTree round trip

## Two traps found on the way

**`compare_apo_ensembles.frame_times_ns` is an apo-only time axis.** It
reconstructs the apo two-stride schedule (200 frames at 0.05 ns, then 0.5 ns).
20 of the 60 holo replicates hold more than 200 frames, for which it returns
anything from 29 to 180 ns against a 100 ns step record. The holo axis is
uniform and lives in each run's own JSON (`analysis_n_frames`,
`analysis_dt_ps`, repaired 2026-08-13), cross-checked against the step record;
`compute_block_metrics.frame_times` does this. Published Table 4 and Figure 5
are unaffected because they average over all frames and never need the axis.

**No holo `*_analysis_pbcfix.dcd` sidecar exists on this machine**, so
`compute_pocket_volume` would skip all 60 holo replicates. The bound volumes
are therefore recomputed with in-memory imaging on the convention
`compare_apo_ensembles` uses for the apo volumes, which reproduces the
published per-replicate means within 0.42 of their own SD.

**G190E holo rep_03** reached 50,000,000 steps but its analysis DCD holds 139
frames at the same 558.659 ps stride as its siblings, so it covers 77.6 ns of
a 100 ns run. The DCD is truncated, not re-strided.

## Cleanup, and a naming collision left in place

`src/nnrti/cli/plot_rt_metric_distributions.py` and its five-page
`supplementary_figure_6_all_rt_metrics_*` outputs are removed: they drew
violins over three replicate means and are superseded by
`plot_rt_metric_violins.py`. The script was untracked, so git could not have
recovered it; a copy is in `published-tables-backup-20261105/superseded/`.

Not touched, but worth knowing about:

- **Two tracked scripts write the same file.**
  `plot_pocket_volumes_from_si.py` and `plot_pocket_volume_violins.py` both
  save `supplementary_figure_6_pocket_volume.{png,pdf}`. Whichever runs last
  wins and nothing warns. Neither produced the submitted art.
- **`supplementary_figure_5_pocket_volume.png` is the art that was actually
  submitted as Supplementary Figure 6** (149,019 bytes, byte-identical to
  `word/media/image6.png` in the pre-2026-11-05 SI). The name is one
  renumbering behind. Keep it: it is the provenance of the submitted figure.
- The new figure is deliberately named
  `supplementary_figure_rt_metric_violins.*` rather than
  `supplementary_figure_6_*`, so it cannot join that collision.

## Crystallographic chi1 anchors (added 2026-11-05)

The two chi1 panels now carry the experimental values from the five wild-type
RT-NNRTI structures of Supplementary Section 1
(`results/analysis/nnrti_pose_comparison/tables/nnrti_tyrosine_comparison.csv`):

| Structure | Ligand | Tyr181 chi1 | Tyr188 chi1 |
|---|---|---|---|
| 4NCG, 2.58 A X-ray | doravirine | -96.0 deg | -49.4 deg |
| 7Z2G, 3.65 A cryo-EM | doravirine | -86.4 deg | -46.1 deg |
| 2ZD1, 1.80 A X-ray | rilpivirine | **+170.8 deg** | -63.2 deg |
| 1FK9, 2.50 A X-ray | efavirenz | **+168.9 deg** | -62.4 deg |
| 1VRT, 2.20 A X-ray | nevirapine | **+173.8 deg** | -69.1 deg |

Three things follow.

**The simulation reproduces the crystallographic rotamers.** Wild-type bound
MD gives Tyr181 -87.0 deg against -96.0 and -86.4 in the two doravirine
structures, and Tyr188 -44.2 deg against -49.4 and -46.1. Both within 9 deg.

**Doravirine is the only one of the four NNRTIs that holds Tyr181
gauche-minus.** Rilpivirine, efavirenz and nevirapine all hold it trans, about
100 deg away. A bound-state Tyr181 trans excursion in the MD is therefore the
side chain moving into the rotamer the other NNRTIs enforce and doravirine
does not -- a concrete structural statement in place of "the chi1 mean shifted
by -76 deg".

**Tyr188 trans appears in no structure and in no bound frame.** All five
crystals are gauche-minus, between -46 and -69 deg, and the trans rotamer is
absent from all 14,229 bound MD frames across 57 replicates. The crystals
corroborate what the simulations show.

Scope, honestly: bound Tyr181 trans occupancy does **not** track resistance
(Spearman rho = 0.24, p = 0.34 over 18 genotypes). A98G+F227C is the one
genotype that is both extreme and consistent, 0.851 +/- 0.075 across all three
replicates against 0.026 for wild type, at 93x. G190S (5.2x) and V106A (9.6x)
are also high; V106A+P225H (153x) and K103N+M230L (36x) are at baseline. This
is a single-genotype observation, not a panel trend.

## Candidate alternate mechanisms (reframed 2026-11-05)

The question is not which metric correlates with resistance. It is whether
these distributions suggest *why* ddG_bind reads near zero for a genotype the
clinic calls resistant, and the reverse.

Using the FEP ddG_bind (the numbers the reviewer's comment quotes), a
regression of log10(fold) on ddG_bind over 18 genotypes is weak (r = 0.28,
p = 0.26). The residual of that fit is "resistance the binding free energy
does not account for".

### The one pattern that survives

Genotypes carrying a substitution in the primer grip (p66 225-235) are
systematically more resistant than their ddG_bind predicts:

| | mean residual | n |
|---|---|---|
| primer-grip substitution | **+0.63** | 7 |
| no primer-grip substitution | **-0.40** | 11 |

Welch t = 4.16, p = 0.0008; Mann-Whitney p = 0.0028. In an OLS on the residual
with both primer-grip membership and induced-fit volume change as predictors,
primer grip carries it (+1.20, p = 0.003) and induced fit does not (p = 0.46);
R^2 = 0.51.

Held against its own background, every available pair agrees:

| background | with primer-grip substitution | delta |
|---|---|---|
| V106A | V106A+F227L | +0.77 |
| V106A | V106A+L234I | +1.02 |
| V106A | V106A+P225H | +1.00 |
| K103N | K103N+M230L | +1.47 |
| K103N | K103N+P225H | +0.68 |
| V106I | V106I+F227C | +2.03 |

All six positive (sign test p = 0.031), mean +1.16 log10 units, about 15-fold
of resistance unaccounted for. The effect survives dropping F227 and L234,
which line the NNIBP as well as the primer grip: P225H and M230L alone give
+1.05.

**The structural correlate moves the same way.** In all six pairs the
unliganded primer-grip Ca RMSF rises relative to its background: V106A 0.676 ->
0.698 / 0.774 / 0.702; K103N 0.602 -> 0.772 / 0.636; V106I 0.590 -> 0.661.

**The mechanism this suggests** is the one M230L already implies: a primer-grip
substitution acts on template/primer engagement, and the binary RT-DOR complex
simulated here contains no template/primer, so that route is invisible to
ddG_bind by construction rather than by error. K103N+M230L is the clean case --
+0.59 +/- 0.41 kcal/mol at 36x, with the most mobile unliganded primer grip in
the panel and the largest binding-associated suppression of it.

### What stops this being a conclusion

The primer-grip genotypes are all double mutants, and the panel holds exactly
one double mutant without a primer-grip component: L100I+K103N, whose delta
from K103N is **+0.79** -- inside the +0.68 to +2.03 range of the primer-grip
pairs. With one control pair, "a second mutation adds resistance ddG_bind does
not see" explains the data as well as "a primer-grip mutation does". The two
cannot be separated here, and the primer-grip reading is preferred on
mechanism, not on the statistics.

### What did not hold

ddG_bind also *over*-predicts for V106I (+2.27 kcal/mol, 1.1x) and V106M
(+6.10, 3.4x). Both start from a collapsed unliganded pocket (55.6 and
54.0 A^3 against 76.4 for wild type) and expand on binding (+39.6, +33.5 A^3
of induced fit above wild type), which suggests the cycle charges a steric
cost at a fixed pose that the real pocket relieves by opening. Across the
panel that does not hold up: the residual does not track induced-fit volume
(Spearman rho = -0.26, p = 0.31), and induced fit is not significant beside
primer grip in the joint model. Reported as an observation about two
genotypes, not a mechanism.

## The NNIBP entrance gates (added 2026-11-05)

`src/nnrti/cli/compute_nnibp_gates.py`. Minimum heavy-atom distance between
the residue pairs lining the mouth of the pocket, per frame, both states:
103-188 (the K103N latch), 101-188, 106-181, 100-227, 103-225.

**Why this was needed.** The literature mechanism for K103N is not a change in
cavity volume -- it is the Asn103 side chain reaching Tyr188 and latching the
pocket mouth shut, stabilising the closed unliganded form and slowing
inhibitor entry. The pocket-volume proxy measures the wrong coordinate for
that, which is why K103N looks identical to wild type in Supplementary
Figure 6 panel A (74.2 +/- 10.2 A^3 against 76.4 +/- 9.4, open-state fraction
0.098 against 0.085).

**On the right coordinate the mechanism appears.** Unliganded 103-188 distance:

| | mean | n |
|---|---|---|
| K103N-containing (4 genotypes) | **6.61 +/- 0.51 A** | 12 reps |
| all other mutants (15) | 8.21 +/- 0.19 A | 45 reps |
| wild type | 8.17 +/- 0.81 A | 3 reps |

Welch t = -2.94, p = 0.011; Mann-Whitney p = 0.0017. The fraction of
unliganded frames with the latch closed below 4 A is 0.10 for the K103N
family against 0.007 for other mutants and 0.012 for wild type.

**And it is specific.** Of the five gates, 103-188 is the only one that
separates the K103N family from the rest (p = 0.011; the other four give
p = 0.25, 0.79, 0.95, 0.93). That is the gate the mechanism names.

In the bound state the latch is displaced in every genotype (K103N family
8.31 A, others 9.34 A, wild type 9.00 A), as it must be with the drug in the
pocket. This is a kinetic mechanism -- an entry rate -- and ddG_bind is an
equilibrium quantity, so the cycle cannot represent it even in principle.

A sanity filter drops any frame with a gate wider than 20 A. It fires on 132
frames, all G190E unliganded replicate 2 from 20 ns on: the same unfolding
documented elsewhere, which otherwise reported that genotype's channel floor
as 30.5 A and its primer-grip rim as 29.1 A. Before 20 ns that replicate reads
6.6 and 10.3 A, against 7.7 and 9.6 for wild type.

## Metrics computed but never surfaced

`mechanism_coordinates.csv` holds 17 per-frame columns that reach no figure or
table: `y188_cent_dist`, `y188_mindist`, `y188_ringcent_dist`,
`chl_trp229_interplanar_deg`, `chl_trp229_cent_dist`, `dor_tau2_deg`,
`n_contacts_188_dor`, `n_contacts_188_chl`, `res227_to_dor_mindist`,
`res229_to_dor_mindist`, `dor_cent_to_S105OG`, `dor_cent_to_S105CA`,
`dor_cent_to_res106`, `res103_sc_to_dor`, `res190_to_dor_mindist`,
`res190CA_to_dor_mindist`, `res190_to_K101NZ`. The moiety table adds
`triazolinone`, `chlorocyanophenyl`, `linker` and `whole_ligand` per frame.

Three were checked here:

- **`res190_to_K101NZ`** is computed only for G190E. A Glu190-Lys101 salt
  bridge would be an obvious explanation for that genotype's favourable
  ddG_bind. It does not form: 7.49 +/- 1.33 A, with 0.2% of frames inside
  4 A. The hypothesis is ruled out.
- **`chl_trp229_interplanar_deg`** gives doravirine a second aromatic partner
  besides Tyr188. Wild type 82.6 deg; the largest departures are
  V106A+F227L -6.5 and V106M -5.3 deg. Small.
- **`res227_to_dor_mindist` / `res229_to_dor_mindist`** put two primer-grip
  residues in direct contact with the drug, 3.5 A in wild type. Relevant to
  the primer-grip observation above and not yet used.

Also unused as analyses: `dor_torsions`, `dor_preorganization`,
`apo_collapse`, `pose_coordinates`, `effect_sizes`.

## Heatmap redesign (2026-11-05)

`plot_structural_heatmap.py` now writes four panels instead of two, with two
changes to how every panel reads.

**Columns are blocked by unit.** Angstroms, cubic Angstroms, degrees, contact
counts and fractions each get their own block, separated by a rule and headed
with the unit. One colour scale across mixed units invited a comparison that
cannot be made: +0.1 A and +0.1 contacts are not the same quantity.

**Colour is ternary, not continuous.** Red where the 95% interval is entirely
above wild type, blue where entirely below, white where it includes zero. The
value and its error still print in the cell, so no magnitude is lost; the fill
now carries exactly one claim, and it is one three replicates can support. The
old per-column normalisation shaded cells by relative magnitude, which ranked
differences the data cannot resolve.

**Two new panels, referenced to wild type within a single state** rather than
as induced fit:

| panel | content | resolved cells |
|---|---|---|
| A | RT–DOR contacts, bound, vs WT | 21 / 126 |
| B | RT conformation, induced fit (bound − unliganded), vs WT | 11 / 88 |
| C | RT conformation, **unliganded**, vs WT | 12 / 89 |
| D | RT conformation, **bound**, vs WT | 18 / 106 |

### What panel C shows that the induced-fit panel cancels

Resolved cells per genotype, split by whether the genotype carries a
primer-grip substitution:

| panel | primer-grip genotypes | others |
|---|---|---|
| C, unliganded | **1.14** | 0.27 |
| D, bound | 0.57 | 1.27 |
| B, induced fit | 0.57 | 0.64 |

Primer-grip genotypes depart from wild type about four times as often as the
rest in the unliganded state, slightly less often than the rest in the bound
state, and not at all distinguishably in the induced-fit panel -- which is the
difference of the two, and so cancels exactly this signal. The two genotypes
with a resolved elevation of unliganded primer-grip RMSF are K103N+M230L
(+0.20 +/- 0.06 A) and V106A+L234I (+0.20 +/- 0.07 A), both primer-grip
mutants.

This is the same observation as the resistance-residual analysis above,
reached from the structures alone, and it is the reason a WT-referenced panel
was needed alongside the induced-fit one.
