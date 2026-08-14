# Modern NNIBP MD descriptor suite

Regen metrics:
```bash
python -m src.analysis.cli.compute_modern_md_suite --workers 4
```

Replot / reinterpret (no traj reload):
```bash
python -m src.analysis.cli.replot_modern_md_suite
```

## Plain-language guide

### DCCM / ΔDCCM
**Dynamic Cross-Correlation Matrix.** For every pair of NNIBP Cα atoms, how correlated is their motion over the trajectory?
- **+1** = move in the same direction together
- **−1** = move in opposite directions
- **0** = uncorrelated

The **Δ heatmap** (mutant − WT) is *not* a contact map. Red = that pair became *more* correlated in the mutant; blue = *less*. See `plots/dccm_explainer_V106A.png`.

### H-bonds (`dor_hbond_*`)
**No — F227C is not the biggest real H-bond loss.** The dark blue at residue **227** for F227C genotypes was almost entirely the mutated site itself (Phe227→Cys changes the partner chemistry). Lys103 occupancy stays high (~0.95–0.99) for F227C. See:
- `dor_hbond_delta_heatmap.png` (mutation sites excluded)
- `dor_hbond_K103_occupancy.png`
- `dor_hbond_res227_occupancy.png`
- `dor_hbond_sum_delta_vs_wt.png`

### Ligand RMSF
The old **C1x spike to ~60 Å was a bug** (`md.rmsf` after pocket fit). True C1x RMSF ≈ **0.7 Å**. Fixed plots: `ligand_rmsf_by_genotype.png`, `ligand_rmsf_mean_by_genotype.png`.

### Contact Δ
Old top-8 bars were unreadable. Use `nnibp_contact_delta_matrices.png` (15×15 Δ frequency vs WT).

### PCA
Old plot = one point per *replicate mean structure* → no dynamics. Replaced by:
- `dor_pose_pca_density_by_genotype.png` — **DOR conformation** density in PC1–PC2
- `nnibp_motion_pca_density_by_genotype.png` — **pocket Cα motion** density

Pose “clusters” were scipy kmeans on COM+flexibility fingerprints; denser maps are the better view.

### Pocket volume
`pocket_volume_vs_experiment.png` — same style as the FEP scatter (vs log10 fold, no fit line).
