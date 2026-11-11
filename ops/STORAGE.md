# Simulation data storage

## Locate data

On Sherlock, determine active working directories from job records and logs.
Existing deployments may use either `$SCRATCH/nnrti-mechanisms-git` or
`$SCRATCH/nnrti-mechanisms`; inspect both before recovery or transfer.
Durable code and data locations are separate:

| Location | Contents |
|---|---|
| `$GROUP_HOME/nnrti-mechanisms/code` | Current code checkout |
| `$GROUP_HOME/nnrti-mechanisms/md_runs/` | MD mirrors |
| `$GROUP_HOME/nnrti-mechanisms/fep_raw/<date>/` | FEP archives and checksum manifests |

The project root contains a private `README.md` and `operations/` records for
current jobs and inventories. `archive/` holds preservation snapshots and retired
checkouts. Older checkout paths may be compatibility links; use `code/` for new work.

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

## Preserve a completed project

Run `ops/storage/preserve_sherlock.py --snapshot <unique-name>` inside a Slurm
allocation after checking that no other jobs are writing project data. The tool
inventories the reviewed project directories, refreshes the designated MD
mirror slots while preserving their prior copies, and writes separate archives
for FEP legs plus the remaining project and auxiliary files. Download and GPU
caches are excluded. The current script targets this project's established
Sherlock deployment; review its source list before reuse for another project.

Each archive includes a per-file JSONL inventory with sizes, modification times
and SHA-256 hashes. The archive checksum is verified by rereading the destination.
`*.verified.json` records completed archive verification; `COMPLETE.json` is written
only after every archive succeeds. `.part` files and an absent completion record
mean preservation is unfinished. After diagnosing an interruption, `--resume` rechecks finished archives against
their checksums and source metadata, then restarts unfinished units. Changes to
already preserved sources stop the resume. Allow sufficient batch time for the
full dataset and checksum rereads.
No source is deleted by this tool.

Use the archived member paths and inventories to restore only the inputs required
for a calculation into a fresh scratch workspace. Verify the archive checksum
before extraction. Preserve the existing directory relationships expected by the
workflow, and run analysis on scratch rather than against archived Group Home data.

For a local preservation snapshot, run
`python ops/storage/preserve_local.py --snapshot <unique-name>`. Completed earlier
preservation outputs are retained separately rather than nested into new archives.

After the scratch snapshot is complete and local transfers have finished, run
`ops/storage/finalize_sherlock.py` in Slurm with `--snapshot <name>` and
`--local-transfer <name>.incoming`. It independently verifies transferred archives,
checks scratch membership and modification records, and consolidates the working
copies under `$SCRATCH/nnrti-mechanisms-work/`. Current results remain at `results/`;
historical staging, code and auxiliary outputs move beneath `legacy/`. Links expose
the maintained Group Home code. The two old project-tree paths remain compatibility
links; no scientific file is deleted. Unexpected changes stop consolidation.
