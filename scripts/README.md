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
MD_PRODUCTION_NS=100.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 SHERLOCK_TIME=12:00:00 \
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Monitor:

```bash
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py --target-ns 100.0 --show-incomplete
```

### Step 3 — sync results back

```bash
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=100.0 \
bash scripts/rsync_results.sh pull
```

### Step 4 — analysis

```bash
bash scripts/run_analysis.sh
```

---

## Apo workflow (ligand-free, defaults to all mutations)

### Step 1 — prep locally

```bash
OPENMM_PLATFORM=CPU python -m src.md.dor_md_pipeline_apo
```

This strips DOR from each holo minimized PDB and writes `manifests/apo_md_manifest.csv`.

### Step 2 — push apo assets to Sherlock + submit

```bash
# Push apo assets (parallel, one Duo auth)
SHERLOCK_USER=rsatija bash scripts/rsync_apo.sh push

# Then on Sherlock:
bash scripts/sherlock/submit_apo_md_batched.sh 6 12
```

### Step 3 — sync results back

```bash
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=100.0 \
bash scripts/rsync_apo.sh pull
```

### Step 4 — analysis

```bash
bash scripts/run_apo_analysis.sh
```

---

## Boltz-2 affinity workflow (WT RT + DOR, Sherlock)

### Step 1 — one-time environment install (login node)

```bash
bash scripts/sherlock/boltz/setup_boltz_env.sh
```

No conda env is created; this installs Boltz into your active Python/user site.

### Step 2 — request an interactive GPU

```bash
bash scripts/sherlock/boltz/salloc_boltz_gpu.sh
```

### Step 3 — run Boltz WT affinity

```bash
cd /scratch/users/$USER/nnrti-mechanisms
bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh
```

The run script auto-generates the Boltz input YAML from:
`data/prepared/dor_4ncg/wt_4ncg.cif`

Default output:
`/scratch/users/$USER/nnrti-mechanisms/results/boltz/wt_affinity`

If your Sherlock account needs explicit module loads for Python/CUDA:

```bash
export SHERLOCK_MODULES="python/3.12.1"
```

Optional speed/sampling override:

```bash
BOLTZ_EXTRA_ARGS="--sampling_steps 50 --diffusion_samples 1" \
bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh
```

Low-memory mode for smaller GPUs (for example 12 GB cards):

```bash
BOLTZ_LOW_MEM=1 bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh
```

### Step 4 — run control mutations (K103N, Y181C, V106A, Y318F)

Single mutation:

```bash
BOLTZ_LOW_MEM=1 BOLTZ_EXTRA_ARGS="--seed 1001" \
bash scripts/sherlock/boltz/run_boltz_mutation_affinity.sh K103N
```

Control panel with replicates:

```bash
BOLTZ_LOW_MEM=1 BOLTZ_REPLICATES=10 BOLTZ_SEED_START=1001 \
bash scripts/sherlock/boltz/run_boltz_control_panel.sh
```

Summarize replicate outputs:

```bash
python3 scripts/sherlock/boltz/summarize_boltz_panel.py \
  --glob "/scratch/users/$USER/nnrti-mechanisms/results/boltz/control_panel/*/replicates/affinity_seed*.json" \
  --out-csv "/scratch/users/$USER/nnrti-mechanisms/results/boltz/control_panel/summary.csv"
```

---

## File inventory

| File | Purpose |
|---|---|
| `rsync_and_analyze.sh` | **One-command**: pull completed results + run full analysis (holo + apo when available) |
| `run_analysis.sh` | Full analysis pipeline (holo core + apo/holo comparative analyses when `manifests/apo_md_manifest.csv` exists) |
| `run_apo_analysis.sh` | Apo analysis (PBC fix → tunnel dynamics → DCCM, apo vs holo comparison) |
| `rsync_results.sh` | Push/pull `results/md_runs/` to/from Sherlock (parallel push, one Duo auth) |
| `rsync_apo.sh` | Push/pull `results/md_runs/apo/` to/from Sherlock (same parallel logic) |
| `sherlock/submit_md_batched.sh` | Submit holo MD jobs in batches with queue monitoring |
| `sherlock/submit_apo_md_batched.sh` | Submit apo MD jobs (same logic, targets `results/md_runs/apo/`) |
| `sherlock/report_md_progress.py` | Summarize job completion vs target steps, flag errors |
| `sherlock/test_one_job.sh` | Single-job smoke test for debugging on Sherlock |
| `sherlock/boltz/setup_boltz_env.sh` | Install/update Boltz-2 on Sherlock without creating a conda env |
| `sherlock/boltz/salloc_boltz_gpu.sh` | Request an interactive GPU allocation sized for Boltz tests |
| `sherlock/boltz/make_boltz_wt_input.py` | Generate RT + DOR Boltz affinity YAML from prepared mmCIFs |
| `sherlock/boltz/run_boltz_wt_affinity.sh` | Run WT RT + DOR Boltz affinity prediction in an interactive GPU session |
| `sherlock/boltz/run_boltz_mutation_affinity.sh` | Run one mutation (or WT) Boltz affinity prediction from prepared mutant CIF |
| `sherlock/boltz/run_boltz_control_panel.sh` | Run the default control panel (`K103N Y181C V106A Y318F`) with N replicates |
| `sherlock/boltz/summarize_boltz_panel.py` | Aggregate replicate affinity JSONs into per-mutation mean/SD/95% CI |
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
