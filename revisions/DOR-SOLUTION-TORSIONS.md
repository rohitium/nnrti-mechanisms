# DOR alone in water: what the free-ligand ensemble is good for

Computed 01-Nov-2026 from 2 x 100 ns of doravirine in explicit solvent
(rep_01 still running, ETA ~23:15). Bound ensembles are the existing holo runs,
read only.

Sources: `results/analysis/dor_torsions_solution/tables/`,
`results/analysis/dor_preorganization/tables/`.
Scripts: `nnrti.cli.compute_dor_torsions` (pointed at a separate output dir so
the reported `results/analysis/dor_torsions` is untouched),
`nnrti.cli.compute_dor_preorganization`.

System check: 2,537 atoms, cubic box 2.92-2.96 nm over 25,000 logged frames, so
the half-box stays at 1.46 nm against the 1.0 nm cutoff. The periodic-image
worry the padding comment flagged does not materialise.

## The lead result: binding does not narrow the distributions, it stops the clock

Figures: `results/analysis/dor_torsions_solution/plots/dor_torsion_trajectories.png`
and `..._distributions.png`. Table: `tables/dor_torsion_kinetics.csv`.
Script: `nnrti.cli.plot_dor_torsion_kinetics`.

Counting how often each torsion leaves its bound rotamer, rather than how wide
it sits:

| | free, escapes/100 ns | bound, escapes in 60 replicates | fraction of time in the bound rotamer, free | bound |
|---|---|---|---|---|
| tau1 pyridinone-O | 7.5 | 1 | 66% | 99.97% |
| tau2 O-aryl | 4.5 | 1 | 31% | 99.96% |
| tau3 pyridinone-N-CH2 | 2.5 | **0** | 32% | 100.00% |
| tau4 CH2-triazolinone | 3.5 | 1 | 33% | 99.79% |

**Across 60 bound replicates x 4 torsions x 100 ns -- 24 microseconds of
cumulative torsion-time -- doravirine leaves its bound rotamer three times.**
Free in water, each torsion leaves several times per 100 ns per replicate. At
matched time resolution that is a 150-450x reduction in exchange rate, and the
true factor is larger: the bound runs write frames 4-11x more coarsely than the
solvent runs, and coarse sampling can only undercount transitions.

Two methodological points the figures and script carry, because both would have
produced a wrong answer if left alone:

- **Frame spacing is not matched.** Solvent runs write 2000 frames per 100 ns
  (50 ps); bound runs write 139-540 (200-720 ps), recorded per replicate as
  `analysis_dt_ps` and cross-checked against `production_ns_completed`. Free
  trajectories are subsampled to the bound resolution before any rate is
  quoted. At full resolution the free rates are 3-50 per 100 ns instead of
  2.5-7.5.
- **The first binning was wrong.** Fixed 120 degree sectors put a well edge at
  -60, which is 2.3 degrees from where the bound tau1 distribution is centred,
  so tau1 rattling about its own mean counted as rotamer hopping -- it appeared
  to move as freely bound as free. Occupancy is now a basin centred on the
  bound circular mean, which has no edge near the populated region, and tau1 is
  as arrested as the rest.

Plain-language version: in water, doravirine is a floppy molecule that
interconverts between rotamers every few nanoseconds. In the pocket it stops
doing that entirely -- not slowed, arrested on the timescale of every
simulation in this study, and in all 20 genotypes alike, resistant ones
included. The pocket does not merely prefer one conformer; it removes the
barrier crossings.

## The same thing seen as amplitude

Circular SD of each rotatable dihedral, degrees:

| | bound mean | bound SD (WT) | free SD | rigidification |
|---|---|---|---|---|
| tau1 pyridinone-O ether | -62.3 | 10.0 | 69.1 | 6.9x |
| tau2 O-aryl | -43.1 | 12.3 | 116.0 | 9.5x |
| tau3 pyridinone-N-CH2 | 78.1 | 7.9 | 81.7 | 10.3x |
| tau4 CH2-triazolinone | 94.2 | 10.6 | 67.9 | 6.4x |

This is block-independent: the free SD is 40-120 degrees in every 25 ns quarter
of both replicates, against bound SDs of 8-12. The statement survives any
reasonable slicing of the data.

## What it does to the existing torsion claim

This is the part worth having. The manuscript reports mutant-vs-WT torsion
widths as evidence about ligand flexibility. The free ligand supplies the
reference state those numbers were missing, and it rescales them:

| | bound -> free gap | largest mutational effect | ratio |
|---|---|---|---|
| tau1 | 59.2 | 2.5 | 24x smaller |
| tau2 | 103.8 | 8.1 (V106A+F227L) | 13x smaller |
| tau3 | 73.8 | 2.9 | 26x smaller |
| tau4 | 57.2 | 13.0 (G190S) | 4.4x smaller |

So **every genotype, resistant ones included, holds DOR in essentially the same
internal conformation.** The largest effect any mutation has on ligand rigidity
is a quarter of the bound-to-free gap, and most are a few percent of it.

Plain-language consequence: resistance does not work by letting the ligand
relax. The pocket still clamps DOR's internal geometry in V106A+P225H and
G190E as firmly as in wild type. Whatever these mutations do, they do it to how
well the clamped ligand is accommodated -- its position, its contacts, the
pocket it sits in -- not to the ligand's own conformational freedom.

The G190S tau4 observation survives with context rather than being overturned:
it roughly doubles the linker's spread (10.6 -> 23.6 degrees), which is real,
but even G190S holds the linker about three times more rigidly than water does.

## The order-of-magnitude result: the pre-organisation cost

ddG_bind omits the cost of pre-organising the ligand by construction, because
single-trajectory MM/GBSA extracts the free ligand from the complex frame -- it
is already in its bound conformation and pays nothing. Measuring that term by
counting avoids the minimiser artefact that makes absolute MM/GBSA energies
unusable (~2,200 kcal/mol of drift between 100 and 1,000 relaxation iterations):
a population ratio has no minimiser in it.

Free DOR occupies its bound conformation in 0.1-0.5% of frames, giving
**dG_preorg of roughly 3-4 kcal/mol**. The window-width sweep spans 3.1 to 4.1,
and the per-torsion product agrees with direct counting to within 11-17% where
counts are adequate, so the method is self-consistent.

**Do not quote this to one decimal.** tau3 crosses its barrier only 3-6 times
per 100 ns replicate (~17-33 ns per crossing), and the joint population
multiplies through it. Note that this is not purely a sampling failure: tau3 is
the slowest torsion in the free ligand too (3 escapes per 100 ns against 18-50
for the others), so part of what limits the estimate is a real property of the
molecule. It limits the number nonetheless. Within-trajectory 25 ns blocks give populations from
0.00% to 1.13%; pooled halves differ fivefold. The other three torsions cross
every ~3 ns and are fine. A converged number needs enhanced sampling of the
torsions, not more plain MD -- at this crossing rate it would take microseconds.
The script prints this check and flags the undersampled torsion every run.

rep_01 will not fix it. A third replicate narrows the average; it does not
converge a coordinate with a ~25 ns correlation time.

## What this answers, and what it does not

Answers:
- R2.3's request for torsion-level support, now with a reference state,
  and with a kinetic observable rather than only a width.
- The size of one term ddG_bind leaves out: a few kcal/mol, against a 7.21
  kcal/mol panel range. Worth stating that the omission is not negligible.
- Whether resistance acts through ligand strain. It does not.

Does **not** answer:
- **Why ddG_bind misses G190E.** dG_preorg is a property of the free ligand, so
  it is identical for every genotype and **cancels exactly** in a mutant-minus-WT
  difference. It cannot explain a genotype-specific discrepancy, and should not
  be offered as if it could.

## Caveats

- n=2 pending rep_01. Replicate-to-replicate joint populations differ 2.5-fold
  (0.10% vs 0.25%), which is the honest error bar on the population and is
  smaller than the tau3 sampling problem above.
- The four dihedrals are treated as the coordinate. Strain in bond angles, ring
  puckers and the trifluoromethyl rotor is not counted, so a few kcal/mol is a
  lower bound on total ligand strain.
- The whole result is a torsional population, so it inherits the ligand force
  field's torsion profiles. A 1 kcal/mol error in an OpenFF barrier moves
  dG_preorg directly. Worth one sentence in Methods rather than silence.
- Both replicates were prepared from the same starting structure and
  decorrelate only through initial velocities, so tau3's agreement between them
  (68/2/30 vs 66/2/32% well occupancy) is weaker evidence than it looks.
- Literature placement is worth a check I have not done: ligand strain on
  binding is generally reported at a few kcal/mol for drug-sized molecules
  (Perola & Charifson 2004 is the usual reference), which would put DOR in the
  ordinary range rather than making it unusual. Verify before citing.
