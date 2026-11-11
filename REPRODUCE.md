# Reproducing the November 2026 revision

## Reproduce the paper

After activating the environment and installing the package, run:

```bash
./workflows/05_manuscript_artifacts.sh
```

This exports the submitted package, runs all 11 supported numerical rebuilds,
collects all 22 submitted figures and tables under their final labels, and runs
the tests. Start browsing at `paper/generated/paper/README.md`. The exact submitted
files are in `main/` and `supplementary/`; numerical candidates are in `rebuilt/`.
`paper/generated/reproduction.json` records input/output hashes and package versions.
`paper/generated/validation/` contains the test log, JUnit report and success record.
A nonzero command exit indicates a failed build or verification; inspect the log.
Tests that need an external spreadsheet engine report a skip if none is configured.

For exact export without rebuilding, use `python -m nnrti.cli.manuscript --check --export`.

## Setup and first check

Use `environment.yml` for the local `nnrti-prep` scientific environment, then `python -m pip install --no-deps -e .`. Set `PYTHONNOUSERSITE=1` to isolate it from user-installed packages. Sherlock simulations use site modules rather than Conda; the optional legacy pmx tools use a separate local environment. See [environment guide](docs/ENVIRONMENTS.md). Plotting and workbook generation require the Python dependencies in `pyproject.toml`; structure/trajectory operations additionally require the conda environment. GROMACS switching simulations used the cluster GPU build of GROMACS 2023.1. The environment specifies minor versions; it is not a complete lockfile for every transitive dependency.

```bash
python -m nnrti.cli.manuscript --check --export
python -m nnrti.cli.manuscript --list
./workflows/05_manuscript_artifacts.sh
```

`--check` verifies the frozen submission. `--export` extracts the exact 13 embedded figure images and their captions, exports the main tables as CSV, and copies all seven submitted SI workbooks into `paper/generated/`. This preserves final artwork even where figure composition was manual.

`--rebuild all` regenerates the supported analysis plots and Tables 2, S3 and S4 from derived data. Choose individual targets, for example:

```bash
python -m nnrti.cli.manuscript --rebuild heatmaps mechanisms protein table4
python -m nnrti.cli.manuscript --rebuild dor-violins rt-violins
```

See [paper/ARTIFACTS.md](paper/ARTIFACTS.md) for current numbering. Rebuilt files are candidates for review, not replacements for the submitted package. Fonts, plotting-library versions and manual composition can affect appearance. The build does not rerun MD or FEP, download missing data, or modify the submitted documents. Run from the repository root if invoking an individual legacy CLI directly.

## Statistical provenance

- Phenotypes: `data/DRM-susceptibilities.csv.xlsx` is authoritative; do not substitute manifest annotations.
- FEP: `results/analysis/fep_pmx/panel_ddg.csv`, with per-leg/per-replicate analysis beneath `legs/`. Recombining a subset of targets overwrites a subset panel: always supply the complete intended target list. Re-analysis with `--force` requires raw switching data.
- MM/GBSA: `results/analysis/binding_energy/tables/ddg_full.csv`, WT referenced across three replicates. The manuscript protocol uses 100 snapshots per replicate over the final 75 ns, with 2,000 unrestrained minimization steps. `mutation_component_summary.csv` is not the authoritative current table.
- Structural observables: `results/analysis/structural_metrics/tables/structural_metrics_per_replicate.csv` contains full-precision replicate means for both states. Frames are not independent replicates. Figure code and S4 share metric definitions and circular-angle handling.
- `build_si_structural_tables` writes live Excel formulas for means, SEMs, t statistics, p values and confidence intervals. Open the generated workbook in Excel or recalculate with LibreOffice before using cached values in a viewer. `openpyxl` writes formulas but does not evaluate them. The frozen submitted workbook already contains recalculated caches.
- The S4 flag means **95% CI excludes zero**. Figure highlighting additionally applies the manuscript's minimum interpretable differences; these are distinct criteria.

## Recomputing from trajectories

The numbered workflow scripts describe preparation, MD, FEP, trajectory analysis and artifact generation. Stages 1–4 are expensive and depend on raw files not distributed through Git. Stage 4 is a starting point, not a complete automatic rebuild of every revision-specific analysis; use the artifact map and individual CLI help for apo metrics, frame series and sensitivity calculations.

Do not rerun holo simulations as part of routine manuscript regeneration. Before analysing recovered trajectories, verify completed steps, physical times, PBC sidecar freshness and the matched apo/holo window. The G190E apo pocket filter must be retained: replicate 2 loses 158 of 200 frames to the 10 Å pocket-Cα radius-of-gyration ceiling, with an additional volume-specific filter. The relevant implementations are in `compute_pocket_volume`, `compute_block_metrics` and the structural table builder.

See [MD run data](docs/methods/md-run-data.md) for trajectory validation and [storage operations](ops/STORAGE.md) for transfers and archive coverage. Verify raw-data availability before starting a trajectory-level rebuild.

## Verification

```bash
pytest -q
python -m nnrti.cli.manuscript --check
```

Tests check panel definitions, energy-table consistency, current structural-table inputs and formula/statistical propagation. Submission checksums provide a separate guard against accidental editing of the submitted package.
