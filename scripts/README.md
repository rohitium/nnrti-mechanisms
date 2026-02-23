# scripts/

Shell entrypoints for the two main workflows: **holo MD** and **apo MD**.

## Workflow overview

```
1. Prep structures locally
        ↓
2. Sync to Sherlock + submit MD jobs
        ↓
3. Sync results back
        ↓
4. Run analysis
```

---

## Holo workflow (DOR-bound, 19 mutations)

### Step 1 — prep (already done; only needed for new mutations)

```bash
# All mutations in susceptibility xlsx:
OPENMM_PLATFORM=CPU ~/miniconda3/envs/nnrti-prep/bin/python -m src.structure_prep.preparation

# One specific mutation (e.g. F227C):
OPENMM_PLATFORM=CPU ~/miniconda3/envs/nnrti-prep/bin/python -m src.structure_prep.preparation \
    --mutations F227C
```

### Step 2 — submit to Sherlock

```bash
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Extension reruns (e.g. after syncing back partial results):

```bash
MD_PRODUCTION_NS=10.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 SHERLOCK_TIME=12:00:00 \
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Monitor:

```bash
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py --target-ns 10.0 --show-incomplete
```

### Step 3 — sync results back

```bash
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=10.0 \
bash scripts/rsync_results.sh pull
```

### Step 4 — analysis

```bash
bash scripts/run_analysis.sh
```

---

## Apo workflow (ligand-free, 7 priority mutations)

### Step 1 — prep locally

```bash
OPENMM_PLATFORM=CPU python -m src.md.dor_md_pipeline_apo
```

This strips DOR from each holo minimized PDB and writes `results/apo_md_manifest.csv`.

### Step 2 — push apo assets to Sherlock + submit

```bash
# Push apo assets (parallel, one Duo auth)
SHERLOCK_USER=rsatija bash scripts/rsync_apo.sh push

# Then on Sherlock:
bash scripts/sherlock/submit_apo_md_batched.sh 6 12
```

### Step 3 — sync results back

```bash
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=10.0 \
bash scripts/rsync_apo.sh pull
```

### Step 4 — analysis

```bash
bash scripts/run_apo_analysis.sh
```

---

## File inventory

| File | Purpose |
|---|---|
| `rsync_and_analyze.sh` | **One-command**: pull completed results + run full analysis (holo + apo when available) |
| `run_analysis.sh` | Full analysis pipeline (holo core + apo/holo comparative analyses when `results/apo_md_manifest.csv` exists) |
| `run_apo_analysis.sh` | Apo analysis (PBC fix → tunnel dynamics → DCCM, apo vs holo comparison) |
| `rsync_results.sh` | Push/pull `results/md_runs/` to/from Sherlock (parallel push, one Duo auth) |
| `rsync_apo.sh` | Push/pull `results/apo_md_runs/` to/from Sherlock (same parallel logic) |
| `sherlock/submit_md_batched.sh` | Submit holo MD jobs in batches with queue monitoring |
| `sherlock/submit_apo_md_batched.sh` | Submit apo MD jobs (same logic, targets `results/apo_md_runs/`) |
| `sherlock/report_md_progress.py` | Summarize job completion vs target steps, flag errors |
| `sherlock/test_one_job.sh` | Single-job smoke test for debugging on Sherlock |
| `src/analysis/cli/validate_manuscript_citations.py` | Check that manuscript figure files exist |

---

## Environment note

Analysis scripts (`run_analysis.sh`, `run_apo_analysis.sh`) use the `nnrti-prep`
conda environment which contains MDAnalysis. Never use bare `python -m` for analysis.

## DCD timestamp note

Never trust `ts.time` from trajectories produced by this pipeline — the DCD
header DELTA field was historically corrupt. Always derive frame timestamps as:

```python
time_ps = (frame_idx + 1) * production_ps / n_total_frames
# where production_ps = md_production_steps_completed * 2.0 / 1000.0
```

See `README.md` → "DCD trajectory time metadata" for the full explanation.
