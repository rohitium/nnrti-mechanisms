# NNRTI Resistance Mechanisms

This repository studies HIV-1 NNRTI resistance for doravirine (DOR) using:
- explicit-solvent MD simulations on Sherlock
- ensemble MM/GBSA-style energy decomposition from MD snapshots
- ensemble structural metrics (contacts, H-bonds, pocket volume)
- correlation of computed metrics with phenotypic fold-reduction data

## Protocol basis

This workflow is aligned to:
- Shao et al., PNAS (2009), DOI: `10.1073/pnas.0907304107`
- Supplementary Information for the same work:
  [https://www.pnas.org/doi/10.1073/pnas.0907304107#supplementary-materials](https://www.pnas.org/doi/10.1073/pnas.0907304107#supplementary-materials)

Key adopted ideas:
- staged minimization/heating/production MD
- ensemble averaging over snapshots for free-energy terms
- decomposition into van der Waals, electrostatic, polar solvation (GB), and nonpolar solvation (SA)
- mutation-level correlation to measured susceptibility shifts

## What is no longer used

The alchemical/FEP lambda protocol is deprecated in this repository and should not be used for new runs.

## Current end-to-end workflow

1. Local preparation:
- build WT + mutant structures
- minimize
- generate Sherlock-ready MD assets and manifest

2. Sherlock execution:
- run one explicit-MD job per manifest task (replicate)
- write trajectory (`.dcd`) and final structure

3. Local analysis:
- compute RMSD convergence profiles
- compute MM/GBSA-style components from ensemble snapshots
- compute WT-referenced ΔΔG
- compute ensemble contacts / H-bonds / pocket-volume proxy
- compute correlations vs phenotype from `data/DRM-susceptibilities.csv.xlsx`

## Quick start

### 1) Test Sherlock connectivity

```bash
export SHERLOCK_USER=<sunet_id>
./scripts/orchestrate.sh --test
```

### 2) Full run (prepare + sync + submit + wait + collect + analyze)

```bash
export SHERLOCK_USER=<sunet_id>
./scripts/orchestrate.sh
```

### 3) Collect/analyze only (no job submission)

```bash
export SHERLOCK_USER=<sunet_id>
./scripts/orchestrate.sh --collect-only
```

## Direct CLI usage

### Local prep only

```bash
python -m src.main \
  --prepare-local-openmm-only \
  --replicates 3 \
  --seed 42 \
  --mutation Y188L
```

### Generate SLURM script

```bash
python -m src.main --generate-slurm --use-openmm-module
```

### Collect and analyze

```bash
python -m src.main \
  --collect-results \
  --mmgbsa-snapshots 100 \
  --mmgbsa-discard-fraction 0.25
```

## Main analysis outputs

- `results/mmgbsa_replicate_metrics.csv`
- `results/ddg_full.csv`
- `results/ddg_summary.csv`
- `results/table1_like_energy_components.csv`
- `results/structural_metrics.csv`
- `results/rmsd_ca_profiles.csv`
- `results/rmsd_convergence_summary.csv`
- `results/boundness_qc.csv`
- `results/correlation_analysis.csv`
- `results/plots/all_metrics_vs_fold_reduction.png`
- `results/plots/fig_s1_like_mutation_landscape.png`
- `results/plots/fig_s2_like_ca_rmsd.png`
- `results/plots/boundness_qc_min_distance.png`

## Notes

- For historical compatibility, some filenames still contain `fep` (for example `results/fep_manifest.csv`, `results/fep_runs/`, `src/cluster/fep_worker.py`). These now carry MD tasks/results, not alchemical FEP jobs.
- Local prep is intended to run on CPU (`OPENMM_PLATFORM=CPU`) unless explicitly overridden.
- Structural metrics require MDAnalysis in the local analysis environment.
