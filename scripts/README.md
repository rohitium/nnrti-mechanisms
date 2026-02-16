# scripts/ README

This directory contains shell entrypoints around the `src/` Python modules.

## Recommended Sherlock MD submission

The currently preferred launcher is:

```bash
chmod +x scripts/sherlock/submit_md_batched.sh
./scripts/sherlock/submit_md_batched.sh 6 12
```

This submits MD jobs in batches and waits for queue pressure to drop before submitting the next batch.
Current defaults in `submit_md_batched.sh` are set for extension runs (10 ns target, force rerun enabled, skip tasks already at target, 12h walltime).

For extension reruns (e.g., 2 ns to 10 ns) without Snakemake:

```bash
MD_PRODUCTION_NS=10.0 MD_FORCE_RERUN=1 SHERLOCK_TIME=12:00:00 ./scripts/sherlock/submit_md_batched.sh 6 12
```

## Analysis from existing trajectories

Run the local checkpointed analysis pipeline:

```bash
./scripts/run_analysis.sh
```

Force recomputation:

```bash
./scripts/run_analysis.sh --force
```

Pipeline steps executed by `run_analysis.sh`:
1. PBC correction for all `results/md_runs/*/*_analysis.dcd` trajectories
2. metadata collection
3. structural metrics
4. plot generation from available checkpoints
5. MM/GBSA (`src.analysis.cli.compute_mmgbsa_safe`, expensive step)
6. final plot regeneration (includes MM/GBSA/ddG outputs)

MM/GBSA tuning variables:

```bash
MMGBSA_SNAPSHOTS=100
MMGBSA_DISCARD_FRACTION=0.25
MMGBSA_WORKERS=8
```

Notes:
- `MMGBSA_SNAPSHOTS=100` is the default protocol.
- `sample_window_ns=1.0` is used internally to select snapshots from the last 1 ns.
- The discard fraction is only a fallback when timing metadata is unavailable.

## DCD timing metadata note (important)

Older stripped analysis DCDs can report inflated `dt` (for example `~50000 ps/frame`).
This came from interval double-counting in the DCD writer metadata path.

What is fixed now:
- `src/md/openmm/md_protocol.py` now writes stripped DCD timing metadata correctly.
- Analysis readers (`src/md/openmm/mmgbsa.py`, `src/analysis/metrics.py`) now normalize legacy inflated `dt` values and still respect last-window selection logic.

Interpretation guidance:
- For existing legacy DCDs already on disk, raw `u.trajectory.dt` may still look wrong.
- Current analysis code corrects this automatically before applying "last 1 ns" filtering.

## Sync helpers

- `scripts/rsync_results.sh`: rsync full `results/md_runs/` plus `results/md_manifest.csv` between local and Sherlock. Default is `push`; use `pull` to download. With `COMPLETE_ONLY=1` and `pull`, it transfers only replicate directories that reached the target production steps.
- `scripts/rsync_json_results.sh`: rsync only JSON files under `results/md_runs/` (plus `results/md_manifest.csv`). Default is `push`; use `pull` to download.

Both require `SHERLOCK_USER`.

## Sherlock submission scripts

- `scripts/sherlock/submit_md_batched.sh`: submit missing MD jobs in batches (preferred).
- `scripts/sherlock/submit_md_only.sh`: submit one job per prepared system immediately.
- `scripts/sherlock/submit_all_md.sh`: similar direct submit loop for all prepared systems.
- `scripts/sherlock/submit_all_tasks.sh`: SLURM array script that runs one manifest task via `src.md.worker`.
- `scripts/sherlock/submit_serial_tasks.sh`: submit selected manifest task IDs one-by-one, waiting for each to finish.
- `scripts/sherlock/report_md_progress.py`: summarize completion vs target steps, running tasks, and latest-log segfaults.
- `scripts/sherlock/test_one_job.sh`: short interactive sanity test on one prepared system.
- `scripts/sherlock/test_extension_resume.sh`: smoke-test checkpoint resume extension (e.g., 2.0 -> 2.01 ns) on one prepared system.
- `scripts/sherlock/run_md_only.sh`: thin wrapper around `submit_md_batched.sh`.

## Deprecated orchestration

- `scripts/orchestrate.sh` is deprecated.
- Use the Sherlock submit scripts above.
