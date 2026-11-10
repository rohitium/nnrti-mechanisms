# MD run data

Storage inventory below was recorded in October–November 2026, not verified by the repository cleanup. Check current manifests and live job paths before assuming any backup is complete.

Where every class of simulation data lives, which copies are durable, and how
to tell a finished run from an interrupted one. `CLAUDE.md` carries the short
version; this is the reference.

## Layout

```
results/md_runs/<genotype>/rep_NN/
  <genotype>_rep_NN.json                    run record: step counts and artifact paths
  <genotype>_repNN_analysis.dcd             stripped protein + ligand trajectory
  <genotype>_repNN_analysis_topology.pdb    matching topology
  <genotype>_repNN_md_state.csv             energy log, and the authoritative step record
  <genotype>_repNN_md.state.xml/.json       portable resume state; absent once a run finishes
  <genotype>_repNN_md.chk                   binary checkpoint, GPU-specific
  <genotype>_minimized_repNN.pdb            post-minimisation structure
  assets/                                   prepared topology and serialised system
```

Apo (ligand-free) runs follow the same layout under
`results/md_runs/apo/<genotype>/rep_NN/`.

Production used a 2 fs timestep. The target is **100 ns per replicate =
50,000,000 steps**; there is no shorter acceptable length. Analysis frames are
written at 500 ps intervals, giving ~200 frames per 100 ns replicate — but see
the warning about frame counts below.

## Where the data actually lives

Four locations, two of which have confusingly similar names. **Sweep every path
under `/scratch/users/$USER/` before concluding that any data is missing.**

| path | role |
|---|---|
| `$SCRATCH/nnrti-mechanisms-git` | **the canonical Sherlock working tree.** 550 GB: all 60 apo replicates, `results/analysis/fep_pmx`, holo runs, pose ladders. Jobs are submitted against this. Its own git checkout is stale — ignore it |
| `$SCRATCH/nnrti-mechanisms` | historically a mostly purged shell as of 2026-10-31; some resumed jobs may still reference this path. Inspect squeue/logs before treating it as unused |
| `$GROUP_HOME/nnrti-mechanisms-repo` | the durable **code** checkout, kept in sync with GitHub |
| `$GROUP_HOME/nnrti-mechanisms` | durable **data**: `fep_raw/` archives and an `md_runs/` mirror |

On 2026-10-31 an agent searched only `$SCRATCH/nnrti-mechanisms`, found no apo
output, and reported that two days of apo simulation had "left no recoverable
output anywhere on the cluster." All of it was in `nnrti-mechanisms-git`, and
four replicates there were further along than the local copies — so four jobs
were also resumed from staler state than existed. Do not repeat this.

A git checkout in `$GROUP_HOME` protects **code, not trajectories**: git ignores
every `.dcd`. Simulation output needs an explicit copy.

## What is protected, and what is not

`$SCRATCH` is purged 90 days after last modification and is not backed up.
There is no Oak allocation. `$GROUP_HOME` has a 1 TB quota.

| data | size on scratch | durable copy |
|---|---|---|
| pmx FEP legs, tier-3 content | ~447 GB | **yes** — `fep_raw/2026-10-20`, 19 legs, 72 GB compressed, sha256 verified |
| FEP `equil.trr` endpoint trajectories | 230 GB | **no** — tier 3 excludes them by design (see below) |
| FEP `_archive/` superseded attempts | 93 GB | **no** — failed and superseded runs |
| apo MD runs | 9.8 GB | partial — mirrored to `md_runs/apo/` |
| pose_ladder / pose_test | 2.7 GB | mirrored to `md_runs/` |
| holo MD runs | ~1 GB | yes — `md_runs/<GENOTYPE>/` |

## FEP raw-data backup tiers

`ops/slurm/cluster/backup_fep_raw.sh` writes one `.tar.gz` per leg into
`$GROUP_HOME/nnrti-mechanisms/fep_raw/<date>/` with a sha256 manifest, and only
ever reads from scratch. It always runs as a batch job. Check it with
`ops/slurm/cluster/check_backup_status.sh`.

- **Tier 1** — `dgdl.xvg`, `equil.gro`, `system.{gro,top}`, itp, mdp, manifests,
  `analysis/`. Every free energy in the paper can be recomputed from this, and it
  compresses hard because `dgdl.xvg` is text.
- **Tier 2** — adds `equil.trr` and `equil.tpr`. This is most of the bulk, and it
  is what lets you re-extract snapshots and re-run switches without redoing
  equilibration.
- **Tier 3** — tier 1 plus `equil.xtc`/`.tpr`/`.cpt` and the switch seed
  snapshots, but **not** `equil.trr`. The `.trr` is full precision with
  velocities at ~2.0 GB per endpoint against the `.xtc`'s 313 MB, so dropping it
  removes ~120 GB while keeping every ensemble that is actually analysed.
  **This is the tier used for the JCIM revision archive.**

So the `.trr` files on scratch are deliberately unprotected. Losing them costs
the ability to re-run switches from the existing equilibration, not any number
in the paper.

Never archive a leg with jobs writing to it — tar would capture half-written
files. The script refuses to run if it sees queued pmx jobs unless `FORCE=1`.

## Backing up apo and pose data

Same hazard: a replicate with a job in flight must be skipped, not copied. The
backup job resolves live slots from `squeue` job names (`apo_<genotype>_<NN>`),
skips those, and refuses to overwrite a destination whose `md_state.csv` reaches
a higher step than the source. It writes a log and a sha256 manifest of the
state files into `$GROUP_HOME/nnrti-mechanisms/md_runs/`. Re-run it after jobs
finish to pick up the slots that were skipped.

## Telling a finished run from an interrupted one

The run JSON is metadata, not a measurement — `status` reads `"ok"` on truncated
runs, `md_production_steps_completed` is wrong in both directions, and the file
is written only on a clean exit, so a wall-clock kill leaves a full trajectory
and no JSON. Use the maximum `Step` in `*_md_state.csv`, cross-checked against
`*_md.state.json`. A run that finished cleanly has no `*_md.state.*` sidecar.

Resume **appends** to both the DCD and the state CSV, so frame counts exceed the
nominal 200 (up to 380 observed) with duplication at segment boundaries, and the
`Step` column is non-monotonic. Never derive a time axis from frame index or
frame count.

`results/md_runs/apo/README.md` has the full account with the specific
replicates that misled, and a one-liner for computing a replicate's real length.

## General

Prefer writing new files over overwriting existing ones, and verify a checksum
before trusting a copy. Transfers between the Mac and Sherlock go through
`ops/sync/`. Large copies within Sherlock run in a Slurm job, not on a login
node.
