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

## Workflow (Snakemake)

The pipeline is managed by [Snakemake](https://snakemake.github.io/) and runs natively on Sherlock. CPU rules (prep, analysis) execute on the login node; GPU rules (MD) are submitted to the `gpu` partition via the SLURM executor.

### Prerequisites (on Sherlock)

```bash
pip install --user snakemake snakemake-executor-plugin-slurm
```

### Quick start

```bash
# Dry run: see what would execute
snakemake -n

# Full production run on Sherlock
snakemake --profile workflow/profiles/sherlock

# Run with custom config overrides
snakemake --profile workflow/profiles/sherlock --config replicates=1 seed=123

# Visualize the DAG
snakemake --dag | dot -Tpng > dag.png
```

### Pipeline stages

1. **prep_wt_cif / prep_mutant_cif** -- Build WT and mutant structures from 4NCG
2. **prep_replicate** -- Minimize + solvate for each mutation/replicate (CPU)
3. **run_md** -- Heating (10-300K NVT) + production (300K NPT) on GPU via SLURM
4. **collect_and_analyze** -- MM/GBSA, structural metrics, RMSD/COM profiles, DDG
5. **generate_plots** -- All publication figures

### Configuration

All parameters are in `workflow/config.yaml`. Key settings:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `replicates` | 3 | Independent replicates per mutation |
| `seed` | 42 | Base seed for deterministic coordinate jitter |
| `md.production_ns` | 2.0 | Production MD length (ns) |
| `analysis.mmgbsa_snapshots` | 100 | Frames for MM/GBSA decomposition |
| `slurm.partition` | gpu | SLURM partition for MD jobs |

### Resume behavior

- Snakemake automatically skips rules whose output files already exist
- OpenMM checkpoint files (`.chk`) enable mid-simulation resume
- Re-running `snakemake` after a failure picks up where it left off

## Direct CLI usage (without Snakemake)

The `src.main` CLI is still available for development and debugging:

### Local prep only

```bash
python -m src.main \
  --prepare-local-openmm-only \
  --replicates 3 \
  --seed 42 \
  --mutation Y188L
```

### Collect and analyze

```bash
python -m src.main \
  --collect-results \
  --manifest results/md_manifest.csv \
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
- `results/plots/simulation_convergence.png`
- `results/plots/boundness_qc_min_distance.png`

## Notes

- The workflow uses `results/md_runs/` for per-mutation/replicate outputs.
- Local prep runs on CPU (`OPENMM_PLATFORM=CPU`) unless overridden.
- Structural metrics require MDAnalysis in the analysis environment.
- The legacy `scripts/orchestrate.sh` is deprecated; use Snakemake instead.
