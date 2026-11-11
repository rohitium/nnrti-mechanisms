# ops/

Cluster operations. Simulation and analysis code is the `nnrti` package under
`src/`; the entry points are in `workflows/`.

| Path | Contents |
|---|---|
| `slurm/cluster/` | Slurm submission and monitoring for the equilibrium MD arrays, with module-load helpers |
| `slurm/fep/` | Slurm submission for the alchemical free energy calculations, and the GROMACS `.mdp` files |
| `storage/` | Preservation, checksum verification and workspace consolidation |
| `sync/` | rsync helpers for moving trajectories between the cluster and a workstation |
| `maintenance/manuscript/` | Applies the ACS template to a Word draft |

The [environment guide](../docs/ENVIRONMENTS.md) defines the local/cluster split:
prepare and analyse locally in `nnrti-prep`; run Sherlock simulations with
OpenMM or GROMACS site modules. No Conda setup is required on Sherlock.

## On the cluster

```bash
cd "$GROUP_HOME/nnrti-mechanisms/code"
git pull --ff-only
cd "$SCRATCH/nnrti-mechanisms-work"
export PROJECT_ROOT="$PWD"
export PYTHONPATH="$GROUP_HOME/nnrti-mechanisms/code/src:${PYTHONPATH:-}"
```

Submitted jobs set `PYTHONPATH` themselves; interactive shells do not.

## Data management

See [storage operations](STORAGE.md) for deployment paths, transfers and archive
verification, and [MD run data](../docs/methods/md-run-data.md) for completion checks.
Stage only explicit code and provenance paths when working in a simulation checkout;
preserve ignored trajectories and checkpoints.

`docx/` contains reusable document inspection/editing utilities and the optional sensitivity-workbook builder. One-off scripts targeting retired revision drafts are preserved under `archive/legacy-document-edits/` and must not be used against the frozen submission. Manuscript export/rebuild uses `python -m nnrti.cli.manuscript`.
