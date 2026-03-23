# Logistic Regression

This folder contains the consolidated binary logistic classifier for DOR resistance.

## Main result

- Task: `low` (`fold < 10`) vs `high` (`fold >= 10`)
- Cross-validation: `5` stratified folds
- Performance: `accuracy = 0.789`, `balanced_accuracy = 0.800`
- Confusion matrix: `9/9` low classified correctly, `6/10` high classified correctly

The main outputs are:

- `plots/confusion_matrix.png`
- `plots/cv_probability_ranked.png`
- `plots/cv_probability_vs_log10_fold.png`
- `plots/feature_contributions.png`
- `plots/full_model_feature_coefficients.png`

Key tables:

- `tables/cv_summary.csv`
- `tables/cv_predictions.csv`
- `tables/false_negative_cases.csv`
- `tables/full_model_feature_coefficients.csv`
- `tables/selected_feature_matrix.csv`

Current high cases predicted as low:

- `G190E`
- `K103N+M230L`
- `A98G+F227C`
- `V106I+F227C`

The strongest full-model coefficients are currently:

- `ligand_palm_distance_angstrom_repstd`
- `residue_min_distance_LYS101_angstrom_repstd`
- `ligand_pose_rmsd_angstrom_mean`
- `residue_min_distance_PRO95_angstrom_repstd`
- `residue_min_distance_PRO236_angstrom_mean`

## Feature Screening

The `feature_screening/` subfolder keeps only the exploratory artifacts still used for interpretation:

- `tables/mutation_feature_matrix.csv`
- `tables/feature_target_associations.csv`
- `plots/feature_target_associations.png`
- `plots/feature_target_scatter_grid.png`
- `plots/feature_target_scatter/`

## Feature Triplets

The `feature_triplets/` subfolder contains the triplet-style histogram + KDE plots for the main classifier features:

- `ligand_palm_distance`
- `ligand_pose_rmsd`
- `residue_min_distance_LYS101`
- `residue_min_distance_PRO236`
- `residue_min_distance_VAL108`
