# scripts/ README

This directory contains shell entrypoints around the `src/` Python modules.

## Recommended Sherlock MD submission

The currently preferred launcher is:

```bash
chmod +x scripts/sherlock/submit_md_batched.sh
./scripts/sherlock/submit_md_batched.sh 6 12
```

This submits MD jobs in batches and waits for queue pressure to drop before submitting the next batch.

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
1. metadata collection
2. MM/GBSA (`src.analysis.cli.compute_mmgbsa_safe`)
3. structural metrics
4. plot generation

## Sync helpers

- `scripts/rsync_results.sh`: rsync full `results/fep_runs/` and `logs/` from Sherlock with retries.
- `scripts/rsync_json_results.sh`: rsync only JSON files from `results/fep_runs/`.

Both require `SHERLOCK_USER`.

## Sherlock submission scripts

- `scripts/sherlock/submit_md_batched.sh`: submit missing MD jobs in batches (preferred).
- `scripts/sherlock/submit_md_only.sh`: submit one job per prepared system immediately.
- `scripts/sherlock/submit_all_md.sh`: similar direct submit loop for all prepared systems.
- `scripts/sherlock/submit_all_tasks.sh`: SLURM array script that runs one manifest task via `src.md.worker`.
- `scripts/sherlock/submit_serial_tasks.sh`: submit selected manifest task IDs one-by-one, waiting for each to finish.
- `scripts/sherlock/test_one_job.sh`: short interactive sanity test on one prepared system.
- `scripts/sherlock/run_md_only.sh`: Snakemake path for submitting only MD rules.

## Deprecated orchestration

- `scripts/orchestrate.sh` is deprecated.
- Use Snakemake (`workflow/profiles/sherlock`) or the Sherlock submit scripts above.
