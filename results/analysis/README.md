# Analysis Outputs

Active analysis folders:

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

Exploratory folders that were superseded during this cleanup were removed.
