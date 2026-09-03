# results/

Generated output. Everything here is produced by `workflows/`; nothing here is
edited by hand.

| Path | Contents | Produced by |
|---|---|---|
| `md_runs/` | Per genotype and replicate: the stripped analysis trajectory, its topology, the run record and the energy log | `workflows/02` (cluster) |
| `analysis/fep_pmx/` | Non-equilibrium alchemical free energy: per-switch work values, per-leg BAR/CGI/Jarzynski estimates, the ΔΔG panel | `workflows/03` (cluster) |
| `analysis/binding_energy/` | MM/GBSA interface energies per replicate | `workflows/04` |
| `analysis/mechanisms/` | Interface geometry for Table 3 and Figure 3 | `workflows/04` |
| `analysis/pocket_volume/` | NNIBP pocket volume | `workflows/04` |
| `analysis/md_convergence/` | DOR pose RMSD and DOR–RT centre-of-mass distance | `workflows/04` |
| `analysis/dor_susceptibility_bar_chart/` | Susceptibility panel for Table 1 | `workflows/05` |
| `analysis/classification_performance/` | Sensitivity, specificity, MCC, ROC AUC | `workflows/05` |
| `plots/` | `figure1B_dor_schematic.pdf` | `workflows/05` |
| `.checkpoints/` | Cached intermediates. Safe to delete; removing one forces that step to recompute | — |

## Notes

- Directory names under `analysis/` and `md_runs/` are referenced by the run
  manifest and by run configurations stored inside those trees, so renaming one
  invalidates those records.
- Result trees reference each other by copy, not by symlink.
- Trajectories are not tracked in git; `ops/sync/` moves them between the
  cluster and a workstation.
