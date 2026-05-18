# Contact Occupancy Permutation Null

This analysis tests whether WT-referenced DOR contact occupancy shifts are larger than expected from trajectory-to-trajectory variability alone.

Replicate trajectories were represented by their contact-occupancy profiles across the 27 WT-contacted residues displayed in the all-mutation WT-referenced occupancy heatmap. Missing residue/trajectory rows in the replicate contact table were treated as zero-contact observations with the trajectory's frame count, matching the heatmap's zero-fill behavior for absent contacts.

For each permutation, mutation labels were shuffled across the 57 replicate trajectories while preserving the original label counts. Pooled contact occupancies were recomputed for WT and each mutant, and the maximum absolute WT-referenced occupancy shift across all screened residues and mutant comparisons was recorded. This global maximum provides a multiple-comparison-aware empirical null for the heatmap screen.

Generated outputs:

- `plots/mutation_max_shift_vs_permutation_null.png`: strongest residue-level occupancy shift per mutant compared with the global 90th percentile null threshold.
- `tables/mutation_max_shift_permutation_summary.csv`: per-mutation maximum observed shift, top residue, null thresholds, and empirical p-values.
- `tables/candidate_reporter_coordinates.csv`: residue/mutation coordinates exceeding the global 90th percentile null threshold.
- `tables/residue_shift_permutation_results.csv`: all mutation-residue WT-referenced shifts and empirical global p-values.
- `tables/permutation_null_distribution.csv`: retained null maxima for reproducibility.
