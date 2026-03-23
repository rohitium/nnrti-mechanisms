# Binding Energy

This folder contains the final MM/GBSA summary outputs.

## Main plots

- `plots/mmgbsa_components_by_mutation.png`
  mutation-level MM/GBSA component means with replicate SEM error bars
- `plots/mmgbsa_ddg_components_vs_wt.png`
  WT-referenced component shifts with replicate SEM error bars
- `plots/mmgbsa_components_vs_fold_change.png`
  mutation-level component means vs fold reduction on a log-scaled x-axis
- `plots/mmgbsa_ddg_vs_fold_change.png`
  WT-referenced component shifts vs fold reduction on a log-scaled x-axis

## Main tables

- `tables/mutation_component_summary.csv`
- `tables/mutation_ddg_summary.csv`
- `tables/component_vs_fold_change_stats.csv`
- `tables/ddg_vs_fold_change_stats.csv`

## Current readout

The total MM/GBSA signal remains weak against susceptibility:

- total `binding_dg` vs fold change: `R^2 = 0.018`, `p = 0.586`
- total `ddg` vs fold change: `R^2 = 0.018`, `p = 0.586`

Among the components, the clearest remaining association is the electrostatic term:

- `binding_dg_electrostatic` vs fold change: `R^2 = 0.254`, `p = 0.028`
- `ddg_electrostatic` vs fold change: `R^2 = 0.254`, `p = 0.028`

So this folder is best treated as a mechanistic comparison package, not a standalone predictive explanation of susceptibility.
