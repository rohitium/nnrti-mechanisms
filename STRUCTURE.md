# Repository structure

Layout and the conventions that keep it reproducible and provenance-tracked.
Paths are defined once in [`src/nnrti/paths.py`](src/nnrti/paths.py) — **import them, don't
hardcode directory strings.**

## Top-level layout

| Directory | Holds |
|---|---|
| `workflows/` | The five numbered entry points, in the order the study ran. |
| `src/nnrti/` | The Python package: `md/` (simulation protocol + cluster workers), `fep/` (pmx non-equilibrium FEP), `analysis/` (analysis library), `cli/` (only the scripts that build a manuscript artifact), `structure_prep/`, `utils/`, and the `paths.py` registry. |
| `scripts/` | Cluster-side operations only: Slurm submission (`sherlock/`, `fep_pmx/`), rsync helpers, MD metadata audit/repair (`md/`), docx tooling (`manuscript/`). No analysis code. |
| `data/` | Inputs: `data/structures/` (experimental CIF/PDB), `data/ligands/` (SDF), `data/prepared/`, and the susceptibility spreadsheet. |
| `results/` | Manuscript-facing output only (see below). |
| `figures/` | Curated, manuscript-facing figures (distinct from per-run plots under `results/`). |
| `manifests/` | Provenance manifests & logs, including `archive_2026-09-01_manifest.csv`. |
| `manuscript/` | Manuscript drafts and assets; `ARTIFACTS.md` is the artifact authority. |
| `logs/` | Run logs (SLURM stdout/stderr). |
| `docs/`, `tests/`, `envs/` | Docs, tests, environment specs. |

The package is importable as `nnrti` once you `pip install -e .` or set
`PYTHONPATH=src`. Every entry point is `python -m nnrti.<subpackage>.<module>`;
there are no loose analysis scripts left to guess between.

## Results

- `results/md_runs/<genotype>/rep_NN/` — **holo** (DOR-bound) classical MD.
- `results/md_runs/apo/<genotype>/rep_NN/` — **apo** (ligand-free) classical MD.
- `results/md_runs/_archive/` — superseded/stub artifacts, never deleted; every
  move is logged in `manifests/md_archive_log.csv` with a sha256.
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
- **Audit:** `python3 scripts/md/audit_md_metadata.py [--check-checkpoints]` —
  read-only consistency check across every run (JSON vs `state.csv` vs `.chk`).
- **Repair:** `python3 scripts/md/repair_md_metadata.py [--apply]` — dry-run by
  default; corrects paths/step-counts from ground truth, archives (never deletes).

## Sync model

- **git/GitHub is the hub** for code + light results; Mac and Sherlock both stay
  in sync via `git pull`/`push`.
- **rsync only for the heavy MD artifacts** that are too large for git.
