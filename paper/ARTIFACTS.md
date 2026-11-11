# Artifact map — submitted 10 November 2026

The reference package is [revisions/submission-ready/](../revisions/submission-ready/). Current numbering is **6 main figures, 2 main tables, 7 supplementary figures, 7 supplementary tables**. Historical filenames below are retained to avoid breaking analysis paths.

The machine-readable [artifact registry](artifacts.json) contains all 22 final labels.
Run `../workflows/05_manuscript_artifacts.sh` from an activated environment (or the
root-relative path `./workflows/05_manuscript_artifacts.sh` from the repository root).
The resulting `generated/paper/README.md` links every submitted artifact and supported
rebuild; `generated/reproduction.json` and `generated/validation/` record provenance and tests.

## Exact submitted artwork and tables

```bash
python -m nnrti.cli.manuscript --check --export
```

Exports appear in `paper/generated/submitted-figures/` and `submitted-tables/`. `captions.json` records each embedded image's document, internal path, SHA-256 and final caption. Main tables are extracted from the clean manuscript; SI workbooks are copied intact, including formulas. These exports reproduce what is actually in the documents, independently of an editable deck's state.

## Rebuilds from derived numerical data

Use `python -m nnrti.cli.manuscript --rebuild TARGET` (one or more targets), or `./workflows/05_manuscript_artifacts.sh` for every supported rebuild. Outputs go to ignored `paper/generated/`. These commands do not overwrite submitted documents. A rebuild is not guaranteed to reproduce final typography or manually composed artwork pixel for pixel.

Paths in the input column are relative to `results/analysis/` unless stated otherwise.

| Submitted artifact | Target | Numerical input / source | Generated file beneath `paper/generated/` |
|---|---|---|---|
| Figure 1: structure, mutation sites, DOR | Exact export; manual composition | `paper/sources/figure1/`, `paper/sources/Figures.pptx`; `nnrti.cli.build_figure1_pymol` and `plot_dor_schematic` generate components (PyMOL required for structure renders) | `submitted-figures/Figure-1.png` |
| Figure 2: FEP protocol and phenotype | Exact export; manual composition | `nnrti.fep.plot_protocol_schematic`, `nnrti.cli.plot_panel_by_resistance_category`; original composite at `archive/original-submission/figures/Fig2.pdf` | `submitted-figures/Figure-2.png` |
| Figure 3: five RT–DOR metrics | `heatmaps` | `structural_metrics/tables/structural_metrics_per_replicate.csv` | `heatmaps/plots/figure3A_rt_dor_contacts.png` |
| Figure 4: Y188L, V106A, G190 mechanisms | `mechanisms` | `mechanisms/mechanism_coordinates.csv`, `pocket_volume/tables/holo_volume_per_frame.csv` | `Figure-4.png` |
| Figure 5: RT metrics separately in holo/apo | `heatmaps` | Structural per-replicate CSV above | `heatmaps/plots/figure3CD_rt_by_state.png` |
| Figure 6: p66 RMSF relative to WT | `protein` | `protein_level/tables/ca_rmsf_per_residue.csv` | `figure6_protein_level.png` |
| Table 1: phenotypes | Exact export | Authoritative `data/DRM-susceptibilities.csv.xlsx`; `nnrti.cli.plot_dor_susceptibility_bars` produces numerical summaries | `submitted-tables/Table-1.csv` |
| Table 2: energetics | `table2` | `binding_energy/tables/ddg_full.csv`, `fep_pmx/panel_ddg.csv`; row order read from clean manuscript | `tables/Table-2-energetics.csv` |
| SFigure 1: independent experimental RT–DOR structures | Exact export; PyMOL source | `data/structures/nnibp-reproducibility.pml`, `nnrti.cli.compare_nnrti_crystal_poses`, `nnrti_pose_comparison/` | `submitted-figures/Supplementary-Figure-1.png` |
| SFigure 2: perturbed starting poses | `pose` | `pose_coordinates/tables/wt/sampled_frames.csv`, `pose_coordinates/tables/ladder_summary.csv` | `Supplementary-Figure-2.png` |
| SFigure 3: six-panel convergence | `convergence` | `md_convergence/tables/frame_traces.csv` and deposited convergence panel images | `convergence/supplementary_figure_1_convergence.png` |
| SFigure 4: switching work distributions | `work` | `fep_pmx/legs/*/*/rep_*/neq/analysis/{analysis.json,integ_fwd.dat,integ_rev.dat}` | `Supplementary-Figure-4.png` |
| SFigure 5: eight RT–DOR distributions | `dor-violins` | `mechanisms/mechanism_coordinates.csv`, `mechanisms/dor_moiety_contacts_per_frame.csv` | `supplementary_figure_dor_contact_violins.png` |
| SFigure 6: ten RT distributions, both states | `rt-violins` | `chi1_apo_holo/tables/chi1_per_frame.csv`, `pocket_volume/`, `protein_level/`, `protein_level_apo/`, structural per-replicate CSV | `rt-violins/supplementary_figure_rt_metric_violins.png` |
| SFigure 7: binding-associated geometric changes | `heatmaps` | Structural per-replicate CSV above | `heatmaps/plots/figure3B_rt_induced_fit.png` |

## Supplementary workbooks

The seven submitted workbooks can always be recovered intact with `--export`. The following distinguishes numerical regeneration from curated inputs.

| Workbook | Contents and provenance |
|---|---|
| STable 1 | Experimental susceptibility data; input workbook `data/DRM-susceptibilities.csv.xlsx`. Curated submitted layout is preserved. |
| STable 2 | Sequence/preparation records. Curated submitted workbook is the reference; preparation inputs are under `data/` and `manifests/`. |
| STable 3 | Per-replicate MM/GBSA and FEP; rebuild target `table3`. |
| STable 4 | Both-state structural replicate values, summaries, within-state and interaction tests; rebuild target `table4`. All derived statistics are live formulas. |
| STable 5 | Ligand charge sensitivity; source `results/analysis/relaxation_sensitivity/tables/charge_model_comparison.csv` and `ops/docx/build_si_workbooks.py`. |
| STable 6 | Snapshot relaxation sensitivity; source `results/analysis/relaxation_sensitivity/` and `ops/docx/build_si_workbooks.py`. |
| STable 7 | Experimental RT structure comparisons, including merged NNRTI tyrosine/atom sheets; source `crystal_pose_comparison/`, `nnrti_pose_comparison/`, `compare_dor_crystal_poses` and `compare_nnrti_crystal_poses`. Curated workbook layout is preserved. |

`ops/docx/build_si_workbooks.py` is an optional sensitivity-table builder, not part of the default artifact build. It requires checkpoint data excluded from Git. Existing submitted workbooks and derived source CSVs preserve the results when checkpoints are unavailable. Do not run historical document rewrite scripts to regenerate the submission.

## Statistical interpretation

Structural comparisons use independent replicate means, propagated SEMs and Student t intervals (df=4 within a state; df=8 for the four-group interaction). Side-chain χ1 means and differences use circular statistics. STable 4's flag reports CI exclusion of zero; figure highlighting also requires the Methods' minimum interpretable difference. The G190E apo intact-pocket filter is applied upstream and must not be dropped when regenerating source CSVs.

For raw trajectories, recovery limits and safe timing/PBC checks, see [REPRODUCE.md](../REPRODUCE.md) and [MD run data](../docs/methods/md-run-data.md).
