# MD run data

## Layout

```
results/md_runs/<genotype>/rep_NN/
  <genotype>_rep_NN.json                    run record: step counts and artifact paths
  <genotype>_repNN_analysis.dcd             stripped protein + ligand trajectory
  <genotype>_repNN_analysis_topology.pdb    matching topology
  <genotype>_repNN_md_state.csv             energy log
  <genotype>_minimized_repNN.pdb            post-minimisation structure
  assets/                                   prepared topology and serialised system
```

Apo (ligand-free) runs follow the same layout under
`results/md_runs/apo/<genotype>/rep_NN/`.

Production used a 2 fs timestep. Analysis frames are written at 500 ps
intervals, giving ~200 frames per 100 ns replicate.

## Ground truth

A run's progress is the last step recorded in `state.csv`, cross-checked against
the OpenMM checkpoint's `currentStep`. The JSON record is derived and can drift.

```bash
python ops/maintenance/md/audit_md_metadata.py [--check-checkpoints]   # read-only
python ops/maintenance/md/repair_md_metadata.py [--apply]              # dry-run by default
```

`repair_md_metadata.py` corrects paths and step counts from the trajectory and
energy log, moving anything it replaces aside and logging the move with a sha256
in `manifests/md_archive_log.csv`.

## Storage

Trajectories, checkpoints, per-switch `dgdl.xvg` and extracted frames are not
tracked in git and are not kept only on `/scratch`, which is purged after 90
days. They are held on durable storage with checksum manifests and moved with
`ops/sync/`. Prefer writing new files over overwriting existing ones, and verify
a checksum before trusting a copy.
