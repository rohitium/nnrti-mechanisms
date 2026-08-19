# Binding Energy

This folder contains the default MM/GBSA summary outputs. As of 2026-05-14,
the default binding-energy analysis uses the final 20 saved trajectory frames
from each replicate.

## Main plots

- `plots/mmgbsa_ddg_components_vs_wt.png`
  WT-referenced component shifts as bar charts with replicate SEM error bars
- `plots/mmgbsa_total_vs_fold_change.png`
  total MM/GBSA energy vs fold reduction with labeled mutations
- `plots/mmgbsa_vdw_vs_fold_change.png`
  vdW component vs fold reduction with labeled mutations
- `plots/mmgbsa_electrostatic_vs_fold_change.png`
  electrostatic component vs fold reduction with labeled mutations
- `plots/mmgbsa_gb_vs_fold_change.png`
  GB polar solvation component vs fold reduction with labeled mutations
- `plots/mmgbsa_sa_vs_fold_change.png`
  SA nonpolar component vs fold reduction with labeled mutations

## Main tables

- `tables/mutation_component_summary.csv`
- `tables/mutation_ddg_summary.csv`
- `tables/component_vs_fold_change_stats.csv`
- `tables/ddg_vs_fold_change_stats.csv`

## Current readout

These outputs were promoted on 2026-05-14 from:

```text
results/analysis/binding_energy/last20frames/
```

The previous top-level defaults and alternate window/diagnostic analyses were
archived under:

```text
results/archive/2026-05-14_binding_energy_nondefault/
```

Earlier stale cached binding-energy outputs remain archived under:

```text
results/archive/2026-05-13_binding_energy_pre_recompute/
```

After contact screening and the GB parameterisation fix (both 2026-08-18), every
component is null against susceptibility, and now with enough precision for that
to be an informative result rather than an inability to resolve one:

- total `ddg` vs fold change: `R^2 = 0.005`, `p = 0.77`
- `ddg_vdw`: `R^2 = 0.012`, `p = 0.67`
- `ddg_electrostatic`: `R^2 = 0.026`, `p = 0.52`
- `ddg_gb`: `R^2 = 0.0001`, `p = 0.97`
- `ddg_sa`: `R^2 = 0.005`, `p = 0.79`

Mean SEM on the total shift is 0.92 kcal/mol, and no component dominates: mean
|ddG| is vdW 1.33, GB 1.42, elec 1.11, SA 0.04 kcal/mol.

Earlier readouts from this folder are superseded. In particular the electrostatic
term was once reported as the largest association (`R^2 = 0.106`, `p = 0.187`);
that signal did not survive contact screening. The GB term's former dominance did
not survive the screening-factor fix. See `MMGBSA_METHOD_AND_RECOMPUTE.md`.

So this folder is best treated as a mechanistic comparison package, not a
standalone predictive explanation of susceptibility.
