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

## PBC-corrected sidecars

`nnrti.cli.fix_pbc_trajectories` skips existing `*_analysis_pbcfix.dcd` files,
so sidecars written before a run was extended stay stale and nothing warns you.
Compare their frame counts against the raw `*_analysis.dcd` before use. DCDs
are gitignored, so back originals up before passing `--force` — git cannot
recover them.

## Sherlock

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

## Authoritative data sources

- DOR fold-resistance values: `data/DRM-susceptibilities.csv.xlsx`. The
  manifest has drifted from it more than once — re-verify against the
  spreadsheet, not the manifest.
- MM/GBSA components: `results/analysis/binding_energy/tables/ddg_full.csv`.
  `mutation_component_summary.csv` is stale and does not reconcile with the
  manuscript.
- Any re-run must reproduce the published manuscript numbers before its new
  results are reported. If it does not, find the methodological cause first.
