# results/

Generated output. Everything here is rebuilt by `workflows/`; nothing here is a
hand-edited input.

| path | what | rebuilt by |
|---|---|---|
| `md_runs/` | Per genotype and replicate: the stripped analysis DCD, its topology PDB, the run JSON and the energy log. ~6 GB, the deposited trajectory set. | `workflows/02` (cluster) |
| `analysis/fep_pmx/` | The pmx non-equilibrium FEP tree: per-switch work values, per-leg BAR/CGI/Jarzynski estimates, and the ΔΔG panel. | `workflows/03` (cluster) |
| `analysis/binding_energy/` | MM/GBSA interface energies, per replicate. | `workflows/04` |
| `analysis/mechanisms/` | Interface geometry behind Table 3 and Figure 3. | `workflows/04` |
| `analysis/modern_md_suite/` | NNIBP pocket volume (the `V(NNIBP)` column of Table 3). | `workflows/04` |
| `analysis/md_convergence/` | DOR pose RMSD and DOR–RT centre-of-mass distance (Supp. Fig. 1). | `workflows/04` |
| `analysis/dor_susceptibility_bar_chart/` | The susceptibility panel behind Table 1. | `workflows/05` |
| `analysis/classification_performance/` | Sensitivity, specificity, MCC and ROC AUC quoted in the text. | `workflows/05` |
| `plots/` | `figure1B_dor_schematic.pdf` — the one figure written outside an analysis directory. | `workflows/05` |
| `.checkpoints/` | Cached intermediate tables so a rerun does not redo hours of trajectory work. Safe to delete; deleting one forces that step to recompute. |  |

## Why `analysis/` is still here

An earlier plan proposed flattening `results/analysis/fep_pmx` to `results/fep`
and so on. It was not done, deliberately: **912 files under `analysis/fep_pmx/`
and `analysis/binding_energy/` record their own paths** in run configs and
provenance JSONs written at execution time. Renaming the directory would either
falsify those records or leave them pointing at a path that no longer exists.
One redundant path component is a smaller cost than a provenance tree that lies
about where it ran.

`md_runs/` was left alone for the same reason: 60 manifest rows and 81 run JSONs
name it.

## No symlinks between result trees

Two run topology PDBs used to be symlinks into `results/visualization/`.
Archiving that directory silently broke WT replicate 1 and Y188L replicate 1 —
the trajectories loaded fine, the topologies did not. They are real files now. If
you find yourself linking one result tree into another, copy instead.
