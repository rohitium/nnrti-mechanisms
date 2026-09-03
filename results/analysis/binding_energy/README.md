# MM/GBSA binding energy

WT-referenced MM/GBSA interaction energies for the genotype panel, and the
source for Table 2's energy components and Supplementary Table 3.

Each replicate is scored over 100 frames spaced evenly across the final 75 ns of
its production trajectory, with every frame minimized to convergence and none
excluded. Mutants are referenced to the mean of the three WT replicates.

Rebuild with the two commands in [REPRODUCE.md](../../../REPRODUCE.md).

## Tables

| File | Content |
| --- | --- |
| `tables/ddg_full.csv` | Per-replicate components; read by Table 2 and Supp. Table 3 |
| `tables/mutation_ddg_summary.csv` | Per-genotype total, mean ± SEM |
| `tables/mutation_component_summary.csv` | Per-genotype components, mean ± SEM |
| `tables/ddg_vs_fold_change_stats.csv` | Regression of ΔΔE_Total against log₁₀ fold-change |
| `tables/component_vs_fold_change_stats.csv` | The same, per component |

## Plots

`plots/mmgbsa_ddg_components_vs_wt.png` gives the component shifts as bars with
replicate SEM. The `mmgbsa_*_vs_fold_change.png` series plots the total and each
component against fold reduction.
