# MM/GBSA Method and Recompute Notes

## Purpose

This note documents the binding-energy workflow after discovery of stale cached
MM/GBSA results and later promotion of the last-20-frame analysis to the default
manuscript-facing output.

## Inputs

Each replicate is read from `manifests/md_manifest.csv` through the repository result collector. For each run, the MM/GBSA calculation uses:

- the solute-only analysis topology PDB,
- the corresponding analysis DCD trajectory,
- the minimized PDB as a fallback structure input,
- doravirine SDF parameters from `data/ligands/dor.sdf`,
- ligand residue name `2KW`.

The calculation is run in the `nnrti-prep` conda environment with `PYTHONPATH=.`.

## Snapshot Sampling

The current default analysis uses the final 20 saved trajectory frames from each
replicate. This terminal-frame protocol is stored under:

```text
results/analysis/binding_energy/last20frames/
```

and has been promoted to the default top-level binding-energy plots and tables:

```text
results/analysis/binding_energy/plots/
results/analysis/binding_energy/tables/
```

Earlier alternate analyses used 100 trajectory snapshots after discarding the
first 25% of frames, final 1 ns, or final 5 ns windows. These non-default outputs
were archived on 2026-05-14 under:

```text
results/archive/2026-05-14_binding_energy_nondefault/
```

## Per-Snapshot Energy Calculation

For each sampled snapshot:

1. The solute coordinates are loaded from the analysis trajectory and converted from angstrom to nanometer.
2. Hydrogen atoms are locally relaxed while all non-hydrogen atoms are harmonically restrained to the MD snapshot coordinates. This step reduces finite-timestep hydrogen overlap artifacts without moving the heavy-atom trajectory geometry.
3. Separate OpenMM systems are evaluated for the complex, receptor, and ligand.
4. Component binding terms are calculated as:

```text
vdW binding term = E_vdW(complex) - E_vdW(receptor) - E_vdW(ligand)
electrostatic binding term = E_elec(complex) - E_elec(receptor) - E_elec(ligand)
GB polar solvation term = G_GB,polar(complex) - G_GB,polar(receptor) - G_GB,polar(ligand)
SA nonpolar term = G_SA(complex) - G_SA(receptor) - G_SA(ligand)
total MM/GBSA score = vdW + electrostatic + GB polar + SA
```

The GB polar-solvation term is evaluated with OpenMM `GBSAOBCForce`, solvent dielectric 80, solute dielectric 2, charges from the nonbonded force, generic element-based radii, and surface-area energy disabled for the polar-only term.

## Units

OpenMM reports energies in kJ/mol, so `src/md/openmm/mmgbsa.py` and the raw
per-replicate checkpoints are kJ/mol. **All canonical outputs are kcal/mol**,
matching the pmx non-equilibrium FEP results in `scripts/fep_pmx` (which use the
same 4.184 factor; see `analyze_neq.py`).

Conversion happens once, at the canonical-table rebuild boundary, via
`src/analysis/units.py`:

- `KJ_PER_KCAL = 4.184`, applied to the 15 `binding_dg*` energy columns
  (value, `_std`, `_sem`) before WT referencing, so the derived `ddg*` and
  `wt_*` columns are already in kcal/mol.
- The result is stamped in an `energy_units` column, which makes the conversion
  **idempotent** — re-running a rebuild cannot double-divide.
- Unstamped legacy tables are assumed kJ/mol. `read_energy_table()` converts on
  read, so stale top-level CSVs cannot silently reach a kcal-labelled plot.
- Structural metrics (distances, volumes) are merged *after* conversion and are
  never touched.

Select with `--energy-units {kcal/mol,kJ/mol}` on `rebuild_binding_energy_sources`;
the source and target units are both recorded in the rebuild config JSON.

## WT-Referenced Shifts

**Current default (2026-08-18): unmatched WT reference.** Every mutant replicate
is referenced to the *mean* of the three WT production replicates:

```text
WT_ref(component) = mean over WT replicates of component(WT, replicate i)
shift(mutant, replicate i) = component(mutant, replicate i) - WT_ref(component)
```

Positive shifts indicate a less favorable MM/GBSA score relative to WT; negative
shifts indicate a more favorable score.

The reference mode is selected with `--wt-reference {unmatched,matched}` on
`rebuild_binding_energy_sources` and recorded in the rebuild config JSON and in
the `wt_reference_mode` column of `ddg_full.csv`.

### Why the reference changed

The previous default paired mutant and WT by replicate index:

```text
shift(mutant, replicate i) = component(mutant, replicate i) - component(WT, replicate i)
```

WT replicate 2 is a large outlier in the van der Waals term (-5.6 kcal/mol, versus
-64.9 and -64.9 kcal/mol for replicates 1 and 3; total -15.9 +/- 22.0 kcal/mol
across the three). Under index matching, that single trajectory injected a common
~-60 kcal/mol offset into the replicate-2 shift of *every* mutation, and the
resulting between-replicate spread was reported as per-mutation uncertainty.
Mean SEM on the total shift was 26 kcal/mol.

Switching to the WT replicate mean leaves every reported mean **exactly
unchanged** — the mean of `mut_i - WT_i` and the mean of `mut_i - mean(WT)` are
algebraically identical for balanced replicate counts — and it leaves all
fold-change correlation statistics unchanged. What changes is the error bars:
mean SEM on the total shift drops from 26 to 12 kcal/mol, and the spread now
reflects mutant-replicate variability rather than one WT trajectory.

SEMs do not shrink uniformly. They rise for mutations whose own replicate 2 was
being incidentally cancelled by WT replicate 2 (K103N+M230L 9.7 -> 30.7 kcal/mol,
G190E 20.5 -> 24.4 kcal/mol; V106I+F227C 27.8 -> 23.3 kcal/mol is a mild drop). This is the
honest accounting: those genotypes were never as well determined as index
matching made them look.

### Reference uncertainty

The WT reference carries its own uncertainty, which is a **common offset shared
by every mutation** and is therefore excluded from the per-mutation SEM:

```text
total              -15.9 +/- 22.0 kcal/mol
van der Waals      -45.1 +/- 19.8 kcal/mol
electrostatic       -5.4 +/-  0.8 kcal/mol
GB polar solvation  35.2 +/-  1.8 kcal/mol
nonpolar SA         -0.5 +/-  0.0 kcal/mol
```

Because it is common to all rows it does not affect *comparisons between*
genotypes, but it does mean the absolute zero of the total and vdW columns is
poorly determined. The electrostatic column is the one component whose reference
is tightly determined, which is consistent with it being the only component with
a non-negligible (if still insignificant) association with phenotype.

Matched-reference outputs, plots, and the corresponding workbook are archived
under `results/archive/2026-08-18_binding_energy_matched_wt_reference/`.

## Cache Issue Found

The stale-cache issue arose because previous output generation could reuse a complete top-level `results/mmgbsa_replicate_metrics.csv` without validating that it matched current trajectory/topology files. The incremental plot/output path also preferred that top-level file over a newer checkpoint. Historical affected outputs were archived under:

```text
results/archive/2026-05-13_binding_energy_pre_recompute/
```

## Fresh Recompute

The fresh recompute writes to a new checkpoint path first:

```text
results/.checkpoints/.checkpoint_mmgbsa_replicate_metrics_fresh_2026-05-13.csv
```

The fresh file will be compared against the current corrected source before any promotion to canonical outputs.

## Current Default Promotion

On 2026-05-14, the last-20-frame workbook was copied to the manuscript-facing
Supplementary Table 3 path:

```text
manuscript/Supplementary-Table-3.xlsx
```

The previous top-level default plots/tables/config were moved into the
non-default archive before the last-20-frame plots/tables/config were copied into
the top-level default locations.

## Manuscript role (Atanu revision)

**Demote MM/GBSA in the main text.** Per the Atanu feedback plan (§2B):

- Collapse the current MM/GBSA Results block (Table 2 + Figure 2 style content) into
  **one paragraph + SI table** (Supplementary Table 3 / last-20-frame workbook).
- Do **not** lead the resistance story with MM/GBSA totals.
- Delete any hybrid-topology / 5 ns endpoint language from the MM/GBSA section —
  that belongs to the FEP protocol walkthrough, not end-point scoring.
- FEP (tight-SEM subset) and the modern-MD descriptor suite carry the quantitative
  claims; MM/GBSA is supporting context only.

This note is the Methods/SI home for how the scores were computed; the draft
rewrite should cite this path rather than expand MM/GBSA in Results.
