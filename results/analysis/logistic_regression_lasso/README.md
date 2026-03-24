# Logistic Regression Lasso

This folder contains the same low-vs-high resistance classifier as [logistic_regression](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/analysis/logistic_regression), but with `L1` logistic regularization instead of `L2`.

## Main result

- Task: `low` (`fold < 10`) vs `high` (`fold >= 10`)
- Cross-validation: `5` stratified folds
- Performance: `accuracy = 0.579`, `balanced_accuracy = 0.583`

This sparse `L1` variant performs substantially worse than the main structural `L2` model, so the existing [logistic_regression](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/analysis/logistic_regression) folder should remain the preferred classifier.

The full-data fit retains only a small number of nonzero coefficients, centered on the same leading structural features:

- `ligand_palm_distance_angstrom_repstd`
- `residue_min_distance_LYS101_angstrom_repstd`
- `ligand_pose_rmsd_angstrom_mean`
- `residue_min_distance_PRO95_angstrom_repstd`
- `residue_min_distance_TYR188_angstrom_repstd`
- `residue_min_distance_PRO236_angstrom_mean`
