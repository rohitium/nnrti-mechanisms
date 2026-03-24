# Logistic Regression Including Energy

This folder is the parallel binary logistic classifier built from the same mutation-level structural features as [logistic_regression](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/analysis/logistic_regression), with MM/GBSA component means and between-replicate SD terms added to the feature matrix.

## Main result

- Task: `low` (`fold < 10`) vs `high` (`fold >= 10`)
- Cross-validation: `5` stratified folds
- Performance: `accuracy = 0.632`, `balanced_accuracy = 0.644`

This energy-augmented variant performs worse than the main structural replicate-SD model, so it should currently be treated as a comparison analysis rather than the preferred classifier.

Main outputs:

- `plots/confusion_matrix.png`
- `plots/cv_probability_ranked.png`
- `plots/cv_probability_vs_log10_fold.png`
- `plots/cv_probability_vs_fold.png`
- `plots/feature_contributions.png`
- `plots/full_model_feature_coefficients.png`

Key tables:

- `tables/cv_summary.csv`
- `tables/cv_predictions.csv`
- `tables/false_negative_cases.csv`
- `tables/full_model_feature_coefficients.csv`
- `tables/selected_feature_matrix.csv`

Among the MM/GBSA terms, `binding_dg_electrostatic_mean` is the strongest surviving full-model coefficient, but the classifier still remains dominated by the structural features.

## Feature Screening

The `feature_screening/` subfolder stores the exploratory mutation-level matrix and univariate feature-target plots for the combined structural + MM/GBSA feature set.
