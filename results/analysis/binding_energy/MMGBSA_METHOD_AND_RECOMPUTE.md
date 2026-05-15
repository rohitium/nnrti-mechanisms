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

## WT-Referenced Shifts

Replicate-level WT-referenced shifts are calculated by matched replicate:

```text
shift(mutant, replicate i) = component(mutant, replicate i) - component(WT, replicate i)
```

Positive shifts indicate a less favorable MM/GBSA score relative to matched WT; negative shifts indicate a more favorable score.

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
