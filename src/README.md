# src/ layout (consolidated)

Top-level packages are now intentionally limited to:

## 1) `src/structure_prep/`
Structure and mutation preparation logic.
- `config.py`: structure/run specs (e.g., DOR 4NCG)
- `preparation.py`: WT/mutant prep + manifest creation
- `mutation/`: mutation parsing/application/numbering helpers

## 2) `src/md/`
MD execution runtime and cluster-facing utilities.
- `manifest.py`: task schema + CSV IO
- `worker.py`: one-task MD execution
- `openmm/`: OpenMM simulation engine helpers
- `sherlock/`: Sherlock-specific runners/helpers
- `cli/test_md_single.py`: local one-task smoke test

## 3) `src/analysis/`
Post-MD analysis and plotting.
- `susceptibility.py`: phenotype input loading
- `metrics.py`: trajectory-derived structural metrics
- `result_collector.py`: MM/GBSA + merge + correlation pipeline
- `plotting.py`: figure builders
- `cli/`: analysis/plot command entrypoints
- `cli/trim_for_pymol.py`: topology/trajectory trimming helper for visualization
- `cli/fix_pbc_trajectories.py`: batch PBC correction for `*_analysis.dcd`

## 4) `src/utils/`
Shared cross-cutting helpers (paths, CIF parsing, mutation token utils).
