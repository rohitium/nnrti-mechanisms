# Working notes for agents

Project-specific facts that are easy to get wrong. Read before analyzing MD
output or touching the cluster.

## Environment

- Code lives under `src/nnrti/`; `PYTHONPATH=src` is required. Commands are
  `python -m nnrti.cli.<name>`. The old `src.analysis.cli.*` entry points are
  dead.
- **Analysis runs locally**, in the `nnrti-prep` conda environment
  (`~/miniconda3/envs/nnrti-prep/bin/python`). Bare `python` lacks MDAnalysis
  and MDTraj. Sherlock is for MD simulation only.
- Stage 4 of `workflows/04_analysis.sh` is the analysis entry point.

## Holo vs apo

The **holo** simulations are the ones reported in the submitted manuscript. Do
not modify, re-image, or re-run them without being asked — a reviewer response
rests on them unchanged. The **apo** simulations are the revision work. When a
task concerns the apo legs, scope every command to `results/md_runs/apo/` and
not to `results/md_runs/`.

## Measuring MD run length

The run JSON is metadata, not a measurement. `status` reads `"ok"` on truncated
runs, `md_production_steps_completed` is wrong in both directions, and the file
is written only on a clean exit — so any wall-clock kill leaves a full
trajectory and no JSON. Use the maximum `Step` in `*_md_state.csv`, cross-checked
against `*_md.state.json`. 50,000,000 steps = 100 ns at the 2 fs timestep.
See `results/md_runs/apo/README.md` for the full account.

## Time axes

Never derive simulation time from a DCD header, from `MDAnalysis`'s `ts.time`,
or from frame index or frame count. Resumed runs append to the DCD, so frame
counts exceed the nominal 200 (up to 380) with duplication at segment
boundaries. Derive time from the step record.

## Apo analysis windows

`compare_apo_ensembles` compares genotypes on a matched time window set by
`--match-ns`. The module default was `10.0` long after every apo replicate
reached 100 ns, so running the script with no arguments silently produced a
10 ns analysis: WT apo volume 88.7 instead of 76.4, and every genotype off by
15–35 A^3 against the documents. It has caused a wrong recompute three times.
The default is now `100.0`, but check the window the run reports before
believing any apo number, and re-derive the documents' numbers rather than
assuming the defaults still match them.

## PBC-corrected sidecars

`nnrti.cli.fix_pbc_trajectories` skips existing `*_analysis_pbcfix.dcd` files,
so sidecars written before a run was extended stay stale and nothing warns you.
Compare their frame counts against the raw `*_analysis.dcd` before use. DCDs
are gitignored, so back originals up before passing `--force` — git cannot
recover them.

## Sherlock — where the data actually is

There are two directories under `$SCRATCH/` with confusingly similar names.
**Sweep both, and every other path under `/scratch/users/$USER/`, before
concluding that anything is missing.**

| path | what it is |
|---|---|
| `$SCRATCH/nnrti-mechanisms-git` | **the canonical working tree.** 549 GB: all 60 apo replicates, `results/analysis/fep_pmx` (536 GB of raw pmx FEP), holo runs, pose ladders |
| `$SCRATCH/nnrti-mechanisms` | mostly a purged shell, but **not empty and not disposable**: see the warning below |
| `$GROUP_HOME/nnrti-mechanisms-repo` | the durable **code** checkout, kept in sync with GitHub |
| `$GROUP_HOME/nnrti-mechanisms` | older data archive: `fep_raw` (77 GB) and an `md_runs` snapshot that stopped at 2026-10-25 |

On 2026-10-31 I searched only `$SCRATCH/nnrti-mechanisms` and reported that two
days of apo simulation "left no recoverable output." It was all in
`nnrti-mechanisms-git`, and four replicates there were further along than the
local copies. Never draw a conclusion about missing runs from one directory.

**The trap runs both ways.** On 2026-11-01, five apo jobs that were still
running turned out to be resuming from and writing into
`$SCRATCH/nnrti-mechanisms`, not the canonical tree — they were submitted
before the two trees were consolidated, and the resume path is baked into a
running job and cannot be changed. Read a job's own log line
(`Resumed production from portable state <path> at step N`) to learn which tree
it writes to; never assume it is the canonical one. Consequences:

- Judging those runs by the canonical tree's files makes a healthy job look
  dead: its files there sit untouched at the step it resumed from while the
  real progress accumulates in the staging tree.
- Do not delete `$SCRATCH/nnrti-mechanisms` while any apo job is in the queue.
  Check `squeue` and each running job's log first.
- When those jobs finish, merge their slots into the canonical tree, keeping
  whichever copy has the higher max `Step`.

Note also that `$GROUP_HOME/nnrti-mechanisms-repo` is a code checkout, so it
protects code, not simulation output — git ignores every `.dcd`. Trajectories
on `$SCRATCH` are purged 90 days after last modification and are not backed up.

## Sherlock — rules

- `/etc/agents/AGENTS.md` is binding. Read it at the start of every session and
  again after any context compaction.
- Code belongs in `$GROUP_HOME` (never purged), job I/O in `$SCRATCH` (purged
  after 90 days). The checkout at `$SCRATCH/nnrti-mechanisms` has already lost
  its `.git` and all its Python sources to the purge. The durable checkout is
  `$GROUP_HOME/nnrti-mechanisms-repo`.
- Never `git clean` or `git add -A` there; a `git clean -fd` once destroyed
  ~15 GB of raw FEP data.
- Omitting `--time` does not mean unlimited — the `gpu` partition default is
  1 hour and its maximum is 2 days.
- No unattended polling loops against the scheduler or a filesystem.
- Two machines push to `main`; fetch before pushing.

## Durability

`$SCRATCH` is purged 90 days after last modification and is not backed up, and
there is no Oak allocation. Durable copies live in `$GROUP_HOME/nnrti-mechanisms`:
`fep_raw/<date>/` holds per-leg `.tar.gz` archives with sha256 manifests, and
`md_runs/` mirrors the MD runs. Never archive or copy a replicate or FEP leg
that has a job writing to it — the copy captures half-written files. Resolve
live slots from `squeue` first.

`docs/methods/md-run-data.md` is the full reference: every data location, which
copies are durable, the FEP backup tiers and what tier 3 deliberately omits, and
how to tell a finished run from an interrupted one.

## Authoritative data sources

- DOR fold-resistance values: `data/DRM-susceptibilities.csv.xlsx`. The
  manifest has drifted from it more than once — re-verify against the
  spreadsheet, not the manifest.
- MM/GBSA components: `results/analysis/binding_energy/tables/ddg_full.csv`.
  `mutation_component_summary.csv` is stale and does not reconcile with the
  manuscript.
- Any re-run must reproduce the published manuscript numbers before its new
  results are reported. If it does not, find the methodological cause first.
