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

## Storage

Trajectories, checkpoints, per-switch `dgdl.xvg` and extracted frames are not
tracked in git and are not kept only on `/scratch`, which is purged after 90
days. They are held on durable storage with checksum manifests and moved with
`ops/sync/`. Prefer writing new files over overwriting existing ones, and verify
a checksum before trusting a copy.
