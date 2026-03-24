# Binding Energy

This folder contains the final MM/GBSA summary outputs.

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

The total MM/GBSA signal remains weak against susceptibility:

- total `binding_dg` vs fold change: `R^2 = 0.018`, `p = 0.586`
- total `ddg` vs fold change: `R^2 = 0.018`, `p = 0.586`

Among the components, the clearest remaining association is the electrostatic term:

- `binding_dg_electrostatic` vs fold change: `R^2 = 0.254`, `p = 0.028`
- `ddg_electrostatic` vs fold change: `R^2 = 0.254`, `p = 0.028`

So this folder is best treated as a mechanistic comparison package, not a standalone predictive explanation of susceptibility.
