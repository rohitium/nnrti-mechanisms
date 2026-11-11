# Simulation data storage

## Locate data

On Sherlock, determine active working directories from job records and logs.
Existing deployments may use either `$SCRATCH/nnrti-mechanisms-git` or
`$SCRATCH/nnrti-mechanisms`; inspect both before recovery or transfer.
Durable code and data locations are separate:

| Location | Contents |
|---|---|
| `$GROUP_HOME/nnrti-mechanisms-repo` | Code checkout |
| `$GROUP_HOME/nnrti-mechanisms/md_runs/` | MD mirrors |
| `$GROUP_HOME/nnrti-mechanisms/fep_raw/<date>/` | FEP archives and checksum manifests |

These are deployment paths, not an inventory or a guarantee of backup completeness.
Verify current files, manifests and checksums before relying on a copy. Scratch
storage is temporary; a code checkout does not preserve ignored simulation output.

## Transfer and verify

1. Identify source and destination paths and any jobs writing to them.
2. Copy completed, inactive runs; defer active runs until their jobs finish.
3. Compare integration-step records before replacing an existing MD copy.
4. Verify checksums and record the archive contents and coverage.

Use the helpers in `sync/` for workstation transfers. Run large cluster copies
and compression as batch jobs. Preserve source files until verification completes.

## FEP archives

`slurm/cluster/backup_fep_raw.sh` creates per-leg archives and a checksum manifest;
`slurm/cluster/check_backup_status.sh` checks backup status. Review the script's
selected legs and configuration before submission.

| Tier | Included data |
|---|---|
| 1 | Switching work, endpoint structures, topology, parameters, manifests and analysis |
| 2 | Tier 1 plus full-precision endpoint trajectories (`equil.trr`) and run input (`equil.tpr`) |
| 3 | Tier 1 plus compressed endpoint trajectories, run input, checkpoints and switch seed snapshots; excludes `equil.trr` |

Archive coverage determines which calculations can be repeated. Tier 3 omits
full-precision endpoint trajectories; repeating switches requires suitable
velocity-bearing starting states or renewed equilibration. Verify the archived
inputs for the intended calculation rather than assuming every intermediate is retained.

For trajectory format and completion checks, see [MD run data](../docs/methods/md-run-data.md).
