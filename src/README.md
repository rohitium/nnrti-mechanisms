# Scientific package

Install with `pip install -e .`, or set `PYTHONPATH=src`. CLI entry points use `python -m nnrti.cli.MODULE`; the manuscript entry point is `nnrti.cli.manuscript`.

| Package | Purpose |
|---|---|
| `structure_prep/` | Structures, mutation parsing, residue numbering and preparation manifests |
| `md/` | OpenMM preparation, production, minimization and MM/GBSA |
| `fep/` | pmx/GROMACS switching, estimators, corrections and panel assembly |
| `analysis/` | Trajectory metrics, PBC, susceptibility, panel definitions and workbook formulas |
| `cli/` | Analysis, plotting and table commands, including exploratory and legacy utilities |
| `utils/` | Shared structure and path helpers |

`nnrti.paths` defines common directories. [The artifact map](../paper/ARTIFACTS.md) identifies current publication commands; not every CLI produces a current figure.

DCD header times and run JSON status can be misleading. Some legacy analysis routines use timing fingerprints and interpolation; inspect them before analysing new or extended runs. Verify physical time from the step records and restart history rather than assuming frame count gives duration. See [MD run data](../docs/methods/md-run-data.md).
