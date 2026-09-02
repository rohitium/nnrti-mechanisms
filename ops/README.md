# ops/

Cluster operations. Simulation and analysis code is the `nnrti` package under
`src/`; the entry points are in `workflows/`.

| Path | Contents |
|---|---|
| `slurm/cluster/` | Slurm submission and monitoring for the equilibrium MD arrays, with module-load helpers |
| `slurm/fep/` | Slurm submission for the alchemical free energy campaign, GROMACS `.mdp` files, and operational runbooks |
| `sync/` | rsync helpers for moving trajectories between the cluster and a workstation |
| `maintenance/md/` | MD metadata audit and repair. `audit_md_metadata.py` is read-only; `repair_md_metadata.py` is dry-run by default |
| `maintenance/manuscript/` | Word tooling: applying the ACS template and revision rounds |

## On the cluster

```bash
cd nnrti-mechanisms && git pull
export PYTHONPATH="$PWD/src"      # or: pip install -e .
```

Submitted jobs set `PYTHONPATH` themselves; interactive shells do not.

## Notes

- Untracked trajectories and per-switch data sit alongside tracked provenance in
  the cluster checkout, so commits there stage explicit paths; `git clean`,
  `git reset --hard` and `git add -A` will take the untracked data with them.
- `/scratch` is purged after 90 days; durable copies carry a checksum manifest.
- A run's progress is the last step in its `state.csv`.
