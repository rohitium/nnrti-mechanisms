# ops/

Cluster-side operations. No science lives here — the analysis and simulation
code is the `nnrti` package under `src/`, and the entry points are in
`workflows/`.

| path | what |
|---|---|
| `slurm/cluster/` | Slurm submission and monitoring for the equilibrium MD arrays, plus module-load and environment helpers. |
| `slurm/fep/` | Slurm submission for the pmx non-equilibrium FEP campaign, the GROMACS `.mdp` files, and the operational runbooks. |
| `sync/` | rsync helpers moving heavy artifacts between the laptop and the cluster. git carries code and light provenance; rsync carries trajectories. |
| `maintenance/md/` | MD metadata audit and repair. `audit_md_metadata.py` is read-only; `repair_md_metadata.py` is dry-run by default and archives rather than deletes. |
| `maintenance/manuscript/` | Word tooling: applying the ACS template, applying a revision round. |

## Before running anything here

The package moved on 2026-09-01. On the cluster:

```bash
cd nnrti-mechanisms && git pull
export PYTHONPATH="$PWD/src"      # or: pip install -e .
```

Submitted jobs export `PYTHONPATH` themselves; the interactive shell does not.

## Ground rules

- **Never `git clean`, `git reset --hard` or `git add -A` on the cluster.** A
  `git clean -fd` there destroyed ~15 GB of raw FEP data. Stage explicit paths.
- Heavy data is not in git and not only on `/scratch`, which is purged after 90
  days. Trajectories and per-switch data go to durable storage with sha256
  manifests.
- How much a run actually did is its `state.csv` last step, never a JSON claim.
