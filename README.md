# NNRTI Resistance Mechanisms

This repository studies HIV-1 NNRTI resistance for doravirine (DOR) using:
- explicit-solvent MD simulations on Sherlock (OpenMM, AMBER ff14SB + GAFF2)
- ensemble MM/GBSA binding free energy decomposition from MD snapshots
- ensemble structural metrics (Cα RMSD, COM distance, pocket volume, key contacts)
- correlation of computed metrics with phenotypic fold-resistance (FR) data

## Quick reference

| What | Command |
|---|---|
| Submit holo MD jobs | `bash scripts/sherlock/submit_md_batched.sh 6 12` |
| Sync results from Sherlock | `SHERLOCK_USER=rsatija COMPLETE_ONLY=1 bash scripts/rsync_results.sh pull` |
| Run full analysis | `bash scripts/run_analysis.sh` |
| Submit apo MD jobs | `bash scripts/sherlock/submit_apo_md_batched.sh 6 12` |
| Run apo analysis | `bash scripts/run_apo_analysis.sh` |
| Monitor queue | `python3 scripts/sherlock/report_md_progress.py --target-ns 100.0 --show-incomplete` |
| Interactive GPU for Boltz | `bash scripts/sherlock/boltz/salloc_boltz_gpu.sh` |
| Install Boltz (no conda env) | `bash scripts/sherlock/boltz/setup_boltz_env.sh` |
| Run WT RT+DOR Boltz affinity | `bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh` |
| Run one mutant RT+DOR Boltz affinity | `bash scripts/sherlock/boltz/run_boltz_mutation_affinity.sh K103N` |
| Run control panel (10x seeds each) | `bash scripts/sherlock/boltz/run_boltz_control_panel.sh` |

---

## Repository structure

Top-level layout:

| Path | Purpose |
|---|---|
| `data/` | Input data (including susceptibility sheet `DRM-susceptibilities.csv.xlsx`) |
| `src/` | Analysis and MD pipeline code (`python -m src...`) |
| `scripts/` | Sherlock submission/sync helpers and local runner scripts |
| `manifests/` | Holo/Apo manifest CSVs (`md_manifest*.csv`, `apo_md_manifest.csv`) |
| `results/md_runs/` | Holo MD trajectories + JSON/state outputs by mutation/replicate |
| `results/md_runs/apo/` | Apo MD trajectories + outputs |
| `results/analysis/` | Containerized analysis products (tables + plots grouped by analysis type) |
| `results/tables/` | Consolidated CSV outputs moved from `results/*.csv` |
| `results/plots/png/` | Consolidated PNG outputs moved from `results/plots/**/*.png` |
| `logs/` | Run logs and checkpoints |

Triplet contact-story output container:

| Path | Contents |
|---|---|
| `results/analysis/triplet_contact_story_100ns/plots/` | One figure per triplet (top: mean trace across replicates; bottom: pooled occupancy heatmap) |
| `results/analysis/triplet_contact_story_100ns/tables/selection_summary.csv` | Selected story residue + pooled occupancy scores per triplet |
| `results/analysis/triplet_contact_story_100ns/tables/mutation_occupancy.csv` | Mutation-level residue occupancy (unweighted mean + pooled) |
| `results/analysis/triplet_contact_story_100ns/tables/replicate_occupancy.csv` | Per-replicate residue occupancy values |
| `results/analysis/triplet_contact_story_100ns/tables/mean_traces.csv` | Mean/SEM distance traces used in plots |
| `results/analysis/triplet_contact_story_100ns/tables/timing_audit.csv` | Frame count + effective ns used per mutation/replicate |
| `results/analysis/triplet_contact_story_100ns/tables/fold_lookup.csv` | DOR fold-change lookup parsed from susceptibility sheet |
| `results/analysis/triplet_contact_story_100ns/config/triplets.txt` | Triplet definitions used for the run |

## Systems simulated

19 systems are currently in the manifest (3 replicates × 10 ns each).

Clinical categorization follows Stanford HIVDB (Shafer et al.):

| Category | Mutations | Notes |
|---|---|---|
| Wildtype | WT | Reference system |
| **Negative controls** — no meaningful DOR effect; common 1st-gen NNRTI (EFV/NVP) DRMs | K103N, Y181C, G190A | Expected to show no signal vs WT |
| **Canonical DOR DRMs** — single | V106A, V106I, Y188L, Y318F | Well-documented susceptibility shifts |
| **Canonical DOR DRMs** — combinations with large clinical effect | V106I+F227C, A98G+F227C, V106A+F227L, V106A+P225H, V106A+L234I | F227 mutations require V106A/I or A98G co-mutation (epistatic obligates) |
| **Highly suspicious** — limited clinical susceptibility data | G190E, K103N+M230L, K103N+P225H, L100I+K103N | Observed in DOR-treated patients; effect size not well-quantified |
| **Uncertain** — weak or ambiguous evidence | V106M, G190S | May have minor effects; data sparse |

## Protocol basis

Workflow aligned to Shao et al., PNAS (2009) DOI: `10.1073/pnas.0907304107`:
- staged minimization → heating (0→300 K) → NPT production (10 ns, 2 fs timestep)
- ensemble averaging over 100 snapshots from the last 1 ns for MM/GBSA
- decomposition into ΔG_VDW + ΔG_Elec + ΔG_GB + ΔG_SA

## Workflow

### 1. Sherlock MD submission

```bash
# submit in batches (batch size=6, max queued jobs=12)
bash scripts/sherlock/submit_md_batched.sh 6 12

# extension reruns (target 100 ns, skip already-complete)
MD_PRODUCTION_NS=100.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 SHERLOCK_TIME=12:00:00 \
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Monitor completion:
```bash
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py --target-ns 100.0 --show-incomplete
```

### 2. Sync results

```bash
# Sherlock → local (completed replicates only)
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=100.0 \
bash scripts/rsync_results.sh pull
```

### 3. Analysis

Run the full checkpointed analysis pipeline (use nnrti-prep env — MDAnalysis required):

```bash
bash scripts/run_analysis.sh
```

> **Important**: Do not use bare `python -m` — the base Python environment lacks MDAnalysis.
> `run_analysis.sh` handles the correct environment automatically.

## Boltz-2 affinity test on Sherlock (WT RT + DOR)

This is a separate workflow to test Boltz-2 affinity prediction on the WT RT/DOR
complex used in this project.

```bash
# 0) (on Sherlock login node) one-time install (no conda env creation)
bash scripts/sherlock/boltz/setup_boltz_env.sh

# 1) request an interactive GPU session
bash scripts/sherlock/boltz/salloc_boltz_gpu.sh

# 2) from inside the allocation, run the WT job
cd /scratch/users/$USER/nnrti-mechanisms
bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh
```

Notes:
- The WT Boltz input YAML is auto-generated from `data/prepared/dor_4ncg/wt_4ncg.cif`.
- Default output path is `/scratch/users/$USER/nnrti-mechanisms/results/boltz/wt_affinity`.
- If your Sherlock account requires module loads for Python/CUDA, set:
  `SHERLOCK_MODULES="python/3.12.1"` before setup/run scripts.
- On smaller GPUs (for example 12 GB TITAN Xp), run with low-memory mode:
  `BOLTZ_LOW_MEM=1 bash scripts/sherlock/boltz/run_boltz_wt_affinity.sh`.
- Override runtime knobs with `BOLTZ_EXTRA_ARGS`, e.g.
  `BOLTZ_EXTRA_ARGS="--sampling_steps 50 --diffusion_samples 1"`.

### Boltz controls panel (K103N, Y181C, V106A, Y318F)

```bash
# inside an interactive GPU allocation
cd /scratch/users/$USER/nnrti-mechanisms

# one mutation
BOLTZ_LOW_MEM=1 BOLTZ_EXTRA_ARGS="--seed 1001" \
bash scripts/sherlock/boltz/run_boltz_mutation_affinity.sh K103N

# full control panel with replicates
BOLTZ_LOW_MEM=1 BOLTZ_REPLICATES=10 BOLTZ_SEED_START=1001 \
bash scripts/sherlock/boltz/run_boltz_control_panel.sh

# summarize replicate outputs
python3 scripts/sherlock/boltz/summarize_boltz_panel.py \
  --glob "/scratch/users/$USER/nnrti-mechanisms/results/boltz/control_panel/*/replicates/affinity_seed*.json" \
  --out-csv "/scratch/users/$USER/nnrti-mechanisms/results/boltz/control_panel/summary.csv"
```

## Analysis outputs

### Primary result (contact occupancy triplets)

Current primary deliverable lives in:

- `results/analysis/triplet_contact_story_100ns/`

Re-generate it with:

```bash
MPLCONFIGDIR=/tmp/mplconfig PYTHONPATH=. \
/Users/rohitpro/miniconda3/envs/nnrti-prep/bin/python \
  -m src.analysis.cli.plot_triplet_contact_story \
  --manifest manifests/md_manifest.csv \
  --susceptibility-xlsx data/DRM-susceptibilities.csv.xlsx \
  --window-ns 100 \
  --contact-cutoff 4.0
```

Key files:

- `results/analysis/triplet_contact_story_100ns/plots/`
- `results/analysis/triplet_contact_story_100ns/tables/selection_summary.csv`
- `results/analysis/triplet_contact_story_100ns/tables/timing_audit.csv`

### Per-replicate data
| File | Contents |
|---|---|
| `results/tables/holo/mmgbsa_replicate_metrics.csv` | MM/GBSA ΔG components per replicate |
| `results/tables/holo/ddg_full.csv` | ΔΔG vs WT + structural metrics merged |
| `results/tables/holo/structural_metrics.csv` | ensemble-averaged contacts, H-bonds, pocket volume |
| `results/tables/holo/rmsd_ca_profiles.csv` | Cα RMSD over time (200-frame resolution) |
| `results/tables/holo/com_distance_profiles.csv` | DOR–pocket COM distance over time |
| `results/tables/holo/pocket_volume_profiles.csv` | NNBP pocket volume over time |
| `results/tables/holo/boundness_qc.csv` | DOR minimum distance QC per replicate |
| `results/tables/holo/drm_sidechain_distance_timeseries_all_mutations.csv` | DRM sidechain↔DOR distances |

### Plots
| Path | Description |
|---|---|
| `results/plots/png/all_metrics_vs_fold_reduction.png` | Scatter grid: all metrics vs FR |
| `results/plots/png/rmsd_convergence.png` | Cα RMSD convergence by mutation |
| `results/plots/png/com_distance_convergence.png` | COM distance convergence |
| `results/plots/png/boundness_qc_min_distance.png` | DOR boundness QC |
| `results/plots/png/drm_distances/` | Per-mutation DRM sidechain↔DOR distance traces |
| `results/plots/png/dor_key_contacts/` | Per-mutation crystal-derived DOR contact distances |
| `results/plots/png/pocket_volume_timeseries/` | Per-mutation NNBP pocket volume traces |
| `results/plots/png/resistance_heatmap.png` | Mutation × metric heatmap |
| `results/plots/png/manuscript_global_signatures.png` | Global signature figure |

## DCD trajectory time metadata — known issue and fix

> This is a subtle but important correctness issue that affects all time-resolved analysis.

### Root cause

`_StrippedDCDReporter` in `src/md/openmm/md_protocol.py` was calling:
```python
app.DCDFile(handle, topology, dt=timestep_ps, interval=25000)
```
OpenMM's `DCDFile` ignores the `interval` argument and writes `NSAVC=1` plus `DELTA = timestep_ps` (in AKMA units) into the binary DCD header. MDAnalysis reads back `DELTA ≈ 0` (corrupt/near-zero AKMA value) and falls back to `dt = 1.0 ps`. Every frame's timestamp becomes just its frame index in picoseconds — e.g. `ts.time = 0, 1, 2, ..., 199` ps instead of `0, 50, 100, ..., 9950` ps for a 10 ns / 200-frame DCD.

This meant convergence plots showed only 200 ps of data instead of 10 ns, and the "last 1 ns" sampling window in MM/GBSA was grabbing the first 1 frame.

### Fix (applied — both writer and readers)

**Writer** (`src/md/openmm/md_protocol.py`): now passes the per-frame time step directly and omits `interval`:
```python
dt_frame = timestep_ps * interval  # e.g. 0.002 ps × 25000 = 50 ps per frame
app.DCDFile(handle, topology, dt_frame)  # no interval kwarg
```

**Readers** (`src/analysis/result_collector.py`): all three profile workers no longer call `ts.time`. Instead they compute time from the JSON output file:
```python
production_ps = md_production_steps_completed * 2.0 / 1000.0  # timestep = 2 fs
time_ps = (frame_idx + 1) * production_ps / n_total_frames
```

**Rule**: never trust `ts.time` or any MDAnalysis-derived timestamp for DCDs produced by this pipeline. Always derive time from the JSON + the formula above.

## Apo simulation pipeline

Apo (ligand-free) simulations test two mechanistic hypotheses:

- **Hypothesis 1** (DCCM): F227C alone disrupts NNBP↔fingers allosteric coupling;
  epistatic co-mutations V106A/A98G restore it.
- **Hypothesis 2** (tunnel dynamics): V106A+P225H and K103N+M230L use a
  kinetic tunnel-opening mechanism visible in the gate distances without DOR.

By default, apo prep now uses WT + all mutations in `data/DRM-susceptibilities.csv.xlsx`.
Use `--mutations ...` to run only a subset (for example the prior 7-priority panel).

### 1. Apo system prep (local)

```bash
OPENMM_PLATFORM=CPU python -m src.md.dor_md_pipeline_apo \
    --holo-runs results/md_runs \
    --apo-runs  results/md_runs/apo \
    --manifest  manifests/apo_md_manifest.csv
```

This strips DOR (resname 2KW) from each holo minimized PDB, builds an amber-only
system, solvates, and writes `manifests/apo_md_manifest.csv`.

### 2. Apo MD on Sherlock

```bash
# Sync apo assets to Sherlock
SHERLOCK_USER=rsatija bash scripts/rsync_apo.sh push

# Submit (same batched workflow)
bash scripts/sherlock/submit_apo_md_batched.sh 6 12

# Monitor
squeue -u $USER
```

### 3. Apo analysis

```bash
# Sync completed apo trajectories back
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=100.0 \
bash scripts/rsync_apo.sh pull

# Run full apo analysis (PBC fix + tunnel dynamics + DCCM, apo and holo comparison)
bash scripts/run_apo_analysis.sh
```

Key output files:

| File | Description |
|---|---|
| `results/tables/apo/apo_nnbp_tunnel_summary.csv` | Gate mean±std per replicate (apo) |
| `results/tables/holo/holo_nnbp_tunnel_summary.csv` | Gate mean±std per replicate (holo, same mutations) |
| `results/tables/apo/apo_dccm_allosteric_coupling.csv` | NNBP↔fingers/palm/thumb coupling scalars (apo) |
| `results/tables/holo/holo_dccm_allosteric_coupling.csv` | Same, holo |
| `results/plots/png/nnbp_tunnel/apo/` | Per-gate distance timeseries + distributions (apo) |
| `results/plots/png/dccm/apo/` | DCCM heatmaps per replicate (apo) |

See `src/README.md` for full argument documentation.

## Notes

- Sherlock GPU allocation: use `salloc -p gpu --gres=gpu:1 --mem=32G` (not `sh_dev -g 1`, which lands on slow MIG-partitioned A30s).
- Local prep runs on CPU (`OPENMM_PLATFORM=CPU`) unless overridden.
- Structural metrics require MDAnalysis (`nnrti-prep` conda env).
