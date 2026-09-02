# `src/nnrti/` layout

The importable package. `pip install -e .`, or `export PYTHONPATH=src`, then
every entry point is `python -m nnrti.<subpackage>.<module>`.

Six subpackages:

## 1) `src/nnrti/structure_prep/`

Structure and mutation preparation logic.

- `config.py` — structure/run specs (e.g., DOR 4NCG)
- `preparation.py` — WT/mutant prep + manifest creation
- `mutation/` — mutation parsing, application, numbering helpers

## 2) `src/nnrti/md/`

MD execution runtime and cluster-facing utilities.

- `manifest.py` — task schema + CSV IO
- `worker.py` — one-task MD execution entry point
- `openmm/` — OpenMM simulation engine helpers:
  - `md_protocol.py` — staged heating + NPT production; `_StrippedDCDReporter`
    (see DCD timestamp note below)
  - `mmgbsa.py` — MM/GBSA single-point energy evaluation with H-relax
  - `platform.py` — GPU/CPU platform selection
- `sherlock/run_md_job.py` — SLURM job entrypoint (called by submission scripts)
- `cli/test_md_single.py` — local one-task smoke test (manual use only)

## 3) `src/nnrti/analysis/` — analysis library

Post-MD analysis and plotting.

- `susceptibility.py` — phenotype input loading (fold-resistance values)
- `metrics.py` — trajectory-derived structural metrics
- `result_collector.py` — MM/GBSA pipeline + profile collection + merge/correlation
- `plotting.py` — figure builders (convergence, scatter grids, pocket volume)

## 5) `src/nnrti/cli/` — manuscript artifact scripts

Only scripts that produce a numbered figure or table live here; there were 87
modules before the 2026-09-01 refactor and there are 16 now. The rest are in the
external archive (see the repository README). `paper/ARTIFACTS.md` maps
each one to its artifact.

## 6) `src/nnrti/fep/` — pmx non-equilibrium FEP

The alchemical pipeline, moved here from `ops/slurm/fep/`. Slurm submission
scripts and `.mdp` files stay under `ops/slurm/fep/`, since they are
cluster-side operations rather than importable code.

### command-line entrypoints

#### Data collection / computation

| Script | Description |
|---|---|
| `analyze_incremental.py` | Checkpointed pipeline runner (`--step collect\|metrics\|plots\|all`) |
| `compute_mmgbsa_safe.py` | Per-replicate MM/GBSA with H-relax, parallelized |
| `compute_structural_metrics_parallel.py` | Contacts, H-bonds, pocket volume (parallel) |
| `compute_nnbp_tunnel_dynamics.py` | NNBP gate distances over time (K101↔Y188, V106↔Y181, etc.) |
| `compute_dccm.py` | Dynamic cross-correlation matrix; NNBP↔domain allosteric coupling |
| `fix_pbc_trajectories.py` | Batch PBC correction for `*_analysis.dcd` files |
| `audit_pbc_trajectories.py` | Thresholded QC audit for residual PBC artifacts after correction |
| `align_trajectories_to_4ncg.py` | Batch-align corrected trajectories to `4NCG.cif` using protein C-alpha atoms |

#### Plotting

| Script | Description |
|---|---|
| `plot_all_mutation_drm_distances.py` | Per-mutation DRM sidechain↔DOR distance traces |
| `plot_all_mutation_dor_key_contacts.py` | Per-mutation crystal-derived DOR contact distances |
| `plot_interesting_drm_distance_traces.py` | Curated high-interest DRM traces (panel figure) |
| `plot_pocket_volume_distributions.py` | NNBP pocket volume distributions |
| `plot_resistance_heatmap.py` | Mutation × metric heatmap |
| `plot_mmgbsa_tables.py` | MM/GBSA component bar charts |
| `plot_mmgbsa_component_signatures.py` | ΔΔG component breakdown by mutation |
| `plot_key_contact_occupancy_heatmap.py` | Crystal-contact occupancy matrix |
| `plot_crystal_dor_key_contacts_2d.py` | 2D contact map from crystal structure (manuscript figure) |

#### Utilities

| Script | Description |
|---|---|
| `curate_interesting_drm_traces.py` | Score and rank DRM traces for PyMOL follow-up |
| `trim_for_pymol.py` | Topology/trajectory trimming helper for visualization |

## 4) `src/nnrti/utils/`

Shared cross-cutting helpers: paths, CIF parsing, mutation token utilities.

---

## DCD timestamp note

> See `ops/README.md` → "DCD trajectory time metadata" for full details.

**Never use `ts.time`** on DCDs produced by this pipeline. OpenMM's `DCDFile` writer
wrote a corrupt DELTA field (near-zero in AKMA units) because `interval` was passed
alongside `dt`. MDAnalysis fell back to `dt = 1.0 ps`, making all `ts.time` values
equal to the bare frame index.

The fix is in `src/nnrti/md/openmm/md_protocol.py` (`_StrippedDCDReporter`) and
`src/nnrti/analysis/result_collector.py` (all three profile workers). Any new code that
reads trajectory timestamps must use:

```python
time_ps = (frame_idx + 1) * production_ps / n_total_frames
# where production_ps = md_production_steps_completed * 2.0 / 1000.0
```
