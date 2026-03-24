# Feature Screening Including Energy

This folder contains the mutation-level feature matrix and exploratory association plots for the structural classifier with MM/GBSA energy terms added.

Compared with the main [logistic_regression feature screening](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/analysis/logistic_regression/feature_screening), the extra features are:

- `binding_dg_mean`
- `binding_dg_repstd`
- `binding_dg_vdw_mean`
- `binding_dg_vdw_repstd`
- `binding_dg_electrostatic_mean`
- `binding_dg_electrostatic_repstd`
- `binding_dg_gb_mean`
- `binding_dg_gb_repstd`
- `binding_dg_sa_mean`
- `binding_dg_sa_repstd`

These are mutation-level MM/GBSA component means and between-replicate SD terms aggregated from [mmgbsa_replicate_metrics.csv](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/mmgbsa_replicate_metrics.csv).
