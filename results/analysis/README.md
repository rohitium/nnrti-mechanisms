# Analysis Outputs

Active analysis folders:

- `pocket_volume/`
  NNIBP descriptor suite (H-bonds, ligand RMSF, pose clusters, pocket volume, PCA, DCCM, contact networks)
- `md_convergence/`
  100 ns coordinate-vs-time convergence traces (RMSD / COM / reporter distances)
- `fep_pmx/`
  NEQ FEP panel, protocol walkthrough figures, discussion tiers
- `triplet_contact_story_100ns/`
  final triplet-contact figures and supporting trace tables
- `ligand_pocket_features/`
  aligned frame-level structural features used as reusable inputs
- `logistic_regression/`
  consolidated low-vs-high resistance classifier, interpretability outputs, and feature triplets
- `logistic_regression_allsd/`
  parallel low-vs-high classifier using pooled all-frame SD features instead of replicate-SD features
- `logistic_regression_lasso/`
  parallel low-vs-high classifier with the same structural inputs as the main model but `L1` logistic regularization
- `logistic_regression_including_energy/`
  parallel low-vs-high classifier with MM/GBSA component means and replicate-SD energy terms added
- `logistic_regression_including_energy_lasso/`
  parallel low-vs-high classifier with the same structural + MM/GBSA inputs but `L1` logistic regularization
- `binding_energy/`
  final MM/GBSA component summaries and fold-change comparison plots

Each analysis folder follows the same container layout:

`results/analysis/<analysis_name>/{plots,tables,config}/`


## What each folder is for

| Folder | Feeds |
| --- | --- |
| `susceptibility/` | Table 1; `tables/dor_susceptibility_values.csv` is the canonical fold-change table, also read by Supp. Table 3 and the FEP panel |
| `binding_energy/` | Table 2 and Supp. Table 3, via `tables/ddg_full.csv` |
| `fep_pmx/` | Table 2, Figure 2 and Supp. Figure 2; `legs/**/analysis/integ_*.dat` are the per-switch work values every ΔΔG is derived from |
| `mechanisms/` | Figure 3 and Table 3 |
| `pocket_volume/` | Table 3's pocket volume |
| `md_convergence/` | Supp. Figure 1 |
| `classification_performance/` | The counts quoted in Results |

Folders hold only what a manuscript artifact reads or what the documented
pipeline writes. Exploratory output is not kept here.
