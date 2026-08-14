# How to read these plots

- `dor_hbond_delta_heatmap.png` — mutation-site residues excluded. F227C’s blue square at 227 was an artifact of mutating the partner residue; see `dor_hbond_res227_occupancy.png` + `dor_hbond_K103_occupancy.png`.
- `ligand_rmsf_by_genotype.png` — **fixed**. Old C1x spike was an `md.rmsf` bug; true CF3/core RMSF is ~0.5–1.5 Å.
- `dor_com_r_alpha_by_genotype.png` — 2D hist of depth $r$ vs polar angle α from WT mean COM axis (better than 1D $r$ alone).
- `dor_com_r_beta_by_genotype.png` — same with azimuth β around that axis.
- `dor_com_hist_by_genotype.png` / `dor_com_radial_*` — ligand COM after Kabsch to WT NNIBP; Δ = lig COM − pocket COM.
- `dor_pose_pca_density_by_genotype.png` — frame-level DOR conformational map (what pose clusters were trying to say).
- `nnibp_motion_pca_density_by_genotype.png` — frame-level pocket Cα motion map (replaces the old per-rep mean-structure PCA scatter).
- `nnibp_contact_delta_matrices.png` — replaces the hard-to-read top-8 bar panel.
- `dccm_explainer_V106A.png` — what DCCM / ΔDCCM means.
- `pocket_volume_vs_experiment.png` — volume vs log10(fold), no fit line.
