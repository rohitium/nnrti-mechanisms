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

The total MM/GBSA signal remains weak against susceptibility in the last-20-frame
analysis:

- total `binding_dg` vs fold change: `R^2 = 0.010`, `p = 0.697`
- total `ddg` vs fold change: `R^2 = 0.010`, `p = 0.697`

Among the components, the electrostatic term is still the largest association,
but it remains weak:

- `binding_dg_electrostatic` vs fold change: `R^2 = 0.106`, `p = 0.187`
- `ddg_electrostatic` vs fold change: `R^2 = 0.106`, `p = 0.187`

So this folder is best treated as a mechanistic comparison package, not a standalone predictive explanation of susceptibility.
