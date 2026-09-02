# Repository structure

Layout and the conventions that keep it reproducible and provenance-tracked.
Paths are defined once in [`src/nnrti/paths.py`](src/nnrti/paths.py) — **import them, don't
hardcode directory strings.**

## Top-level layout

| Directory | Holds |
|---|---|
| `workflows/` | Five numbered entry points, in the order the study ran. |
| `src/nnrti/` | The package: `md/`, `fep/`, `analysis/`, `cli/`, `structure_prep/`, `utils/`, and the `paths.py` registry. Importable as `nnrti` after `pip install -e .` or `PYTHONPATH=src`. |
| `ops/` | Cluster operations, no science: `slurm/{cluster,fep}`, `sync/`, `maintenance/{md,manuscript}`. |
| `data/` | Deposited inputs: structures, ligands, prepared systems, the susceptibility spreadsheet. |
| `results/` | Generated output. See `results/README.md` for the map and for why `analysis/` is still a layer. |
| `paper/` | The deliverable: `submission/`, `tables/`, `sources/`, `ARTIFACTS.md`. |
| `manifests/` | Run manifests and archive inventories. |
| `docs/` | `methods/`, `runbooks/`, `decisions/`. |
| `env/` | Conda specs. The Python `.gitignore` template ignores `env/`, so the root directory is un-ignored explicitly. |
| `tests/`, `logs/` | Regression tests; Slurm stdout/stderr (git-ignored). |

## Results

- `results/md_runs/<genotype>/rep_NN/` — **holo** (DOR-bound) classical MD.
- `results/md_runs/apo/<genotype>/rep_NN/` — **apo** (ligand-free) classical MD.
- Earlier and stub run artifacts are never deleted. `ops/maintenance/md/repair_md_metadata.py`
  moves them out and logs every move in `manifests/md_archive_log.csv` with a sha256.
- `results/analysis/fep_pmx/` — pmx non-equilibrium alchemical FEP.

## Data-safety conventions (learned the hard way)

1. **Raw heavy data is not in git and not only on scratch.** Trajectories,
   checkpoints, per-switch `dgdl.xvg`, extracted frames → durable storage
   (`$GROUP_HOME`) with sha256 manifests. `/scratch` is purged after 90 days.
2. **git tracks light provenance + code only.** Heavy artifacts are ignored (see
   `.gitignore`); the pmx `legs/**` whitelist tracks only `analysis/**` +
   per-unit manifests. **Never run `git clean` over untracked heavy data.**
3. **Write-once / immutable.** Prefer new, segmented output files over appending
   to or overwriting authoritative trajectories and checkpoints. Archive before
   mutating; verify (checksums) before trusting a copy.
4. **Ground truth over metadata.** How much a run actually did = its `state.csv`
   last step (and the checkpoint `currentStep`), never a JSON claim.

## MD run metadata

- Each run dir carries a JSON with step counts + paths, plus `state.csv` (energy
  log) and `_md.chk` (OpenMM checkpoint). Timestep is 2 fs.
- **Audit:** `python3 ops/maintenance/md/audit_md_metadata.py [--check-checkpoints]` —
  read-only consistency check across every run (JSON vs `state.csv` vs `.chk`).
- **Repair:** `python3 ops/maintenance/md/repair_md_metadata.py [--apply]` — dry-run by
  default; corrects paths/step-counts from ground truth, archives (never deletes).

## Sync model

- **git/GitHub is the hub** for code + light results; Mac and Sherlock both stay
  in sync via `git pull`/`push`.
- **rsync only for the heavy MD artifacts** that are too large for git.
