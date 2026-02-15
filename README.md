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

## Workflow (Script-Based)

The production workflow is script-based:
- submit/extend MD runs on Sherlock with `scripts/sherlock/submit_md_batched.sh`
- sync `results/md_runs` + `results/md_manifest.csv` between Sherlock and local
- run checkpointed analysis locally

### Sherlock MD submission

```bash
# submit in batches (batch size=6, max queued jobs=12)
bash scripts/sherlock/submit_md_batched.sh 6 12

# extension-style rerun example (target 10 ns, skip tasks already at target)
MD_PRODUCTION_NS=10.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 SHERLOCK_TIME=12:00:00 \
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Monitor and completion checks:

```bash
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py --target-ns 10.0 --show-incomplete
```

### Sync results

```bash
# local -> Sherlock
SHERLOCK_USER=rsatija bash scripts/rsync_results.sh push

# Sherlock -> local (completed replicates only)
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=10.0 \
bash scripts/rsync_results.sh pull
```

## Direct analysis usage (without Snakemake)

For analysis from existing trajectories in `results/md_runs/`:

```bash
# end-to-end analysis (checkpointed)
./scripts/run_analysis.sh

# or run individual steps
python -m src.analysis.cli.analyze_incremental --step collect
python -m src.analysis.cli.compute_mmgbsa_safe --force --snapshots 100 --sample-window-ns 1.0 --workers 8
python -m src.analysis.cli.analyze_incremental --step metrics
python -m src.analysis.cli.analyze_incremental --step plots
```

## Main analysis outputs

- `results/mmgbsa_replicate_metrics.csv`
- `results/ddg_full.csv`
- `results/structural_metrics.csv`
- `results/rmsd_ca_profiles.csv`
- `results/com_distance_profiles.csv`
- `results/boundness_qc.csv`
- `results/correlation_analysis.csv`
- `results/plots/all_metrics_vs_fold_reduction.png`
- `results/plots/rmsd_convergence.png`
- `results/plots/com_distance_convergence.png`
- `results/plots/boundness_qc_min_distance.png`

## Notes

- The workflow uses `results/md_runs/` for per-mutation/replicate outputs.
- Local prep runs on CPU (`OPENMM_PLATFORM=CPU`) unless overridden.
- Structural metrics require MDAnalysis in the analysis environment.
