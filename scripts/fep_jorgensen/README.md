# Mutation-agnostic Jorgensen-style FEP

This workflow computes doravirine binding effects relative to WT using reusable
single-residue alchemical legs. It covers every mutant in `manifests/md_manifest.csv`.

Single mutants use one leg:

```text
WT -> V106A
```

Compound mutants use two sequential legs and reuse the first-leg result:

```text
WT -> V106A -> V106A+L234I
```

The target free energy and uncertainty are:

```text
ddG(WT -> target) = sum(ddG_leg)
sigma_target = sqrt(sum(sigma_leg**2))
```

## 1. Generate the panel plan locally

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.panel
```

This writes:

- `results/analysis/fep_jorgensen/prepare_all.sh`
- `results/analysis/fep_jorgensen/worker_manifest.csv`

The default panel contains 19 targets, 19 unique alchemical legs, and 418
OpenMM tasks (19 legs x 2 phases x 11 lambda states).

## 2. Prepare hybrid systems locally

Run the generated script in the full `nnrti-fep` environment. Perses,
OpenMMTools, OpenEye, and ligand parameterization packages are needed only for
this stage.

```bash
results/analysis/fep_jorgensen/prepare_all.sh
```

Each leg is serialized under `results/analysis/fep_jorgensen/legs/` as OpenMM
XML, PDB, and JSON inputs.

## 3. Run GPU sampling on Sherlock

Transfer the serialized leg directories, manifest, repository code, and an
OpenMM-only environment to Sherlock, then submit the array:

```bash
./scripts/sherlock/submit_fep_jorgensen_windows.sh
```

The Sherlock workers import only Python's standard library and OpenMM.
Concurrency can be limited with `SHERLOCK_MAX_CONCURRENT`.

## 4. Analyze locally

PyMBAR and NumPy are needed only on the local machine.

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --all-targets
```

Target summaries are written under `targets/`, and the combined table is
`manuscript_panel_summary.csv`. Positive values indicate weaker doravirine
binding than WT.

## Individual calculations

Prepare one arbitrary single-residue leg:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.prepare \
  --mutation Y181C --start-label WT --end-label Y181C \
  --input-complex-pdb results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb
```

Compound mutants must be expressed as sequential single-residue legs; a
compound label must never be passed as `--mutation`.
