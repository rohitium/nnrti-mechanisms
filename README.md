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
| Monitor queue | `python3 scripts/sherlock/report_md_progress.py --target-ns 10.0 --show-incomplete` |

---

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

# extension reruns (target 10 ns, skip already-complete)
MD_PRODUCTION_NS=10.0 MD_FORCE_RERUN=1 SKIP_IF_AT_TARGET=1 SHERLOCK_TIME=12:00:00 \
bash scripts/sherlock/submit_md_batched.sh 6 12
```

Monitor completion:
```bash
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py --target-ns 10.0 --show-incomplete
```

### 2. Sync results

```bash
# Sherlock → local (completed replicates only)
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=10.0 \
bash scripts/rsync_results.sh pull
```

### 3. Analysis

Run the full checkpointed analysis pipeline (use nnrti-prep env — MDAnalysis required):

```bash
bash scripts/run_analysis.sh
```

> **Important**: Do not use bare `python -m` — the base Python environment lacks MDAnalysis.
> `run_analysis.sh` handles the correct environment automatically.

## Analysis outputs

### Per-replicate data
| File | Contents |
|---|---|
| `results/mmgbsa_replicate_metrics.csv` | MM/GBSA ΔG components per replicate |
| `results/ddg_full.csv` | ΔΔG vs WT + structural metrics merged |
| `results/structural_metrics.csv` | ensemble-averaged contacts, H-bonds, pocket volume |
| `results/rmsd_ca_profiles.csv` | Cα RMSD over time (200-frame resolution) |
| `results/com_distance_profiles.csv` | DOR–pocket COM distance over time |
| `results/pocket_volume_profiles.csv` | NNBP pocket volume over time |
| `results/boundness_qc.csv` | DOR minimum distance QC per replicate |
| `results/drm_sidechain_distance_timeseries_all_mutations.csv` | DRM sidechain↔DOR distances |

### Plots
| Path | Description |
|---|---|
| `results/plots/all_metrics_vs_fold_reduction.png` | Scatter grid: all metrics vs FR |
| `results/plots/rmsd_convergence.png` | Cα RMSD convergence by mutation |
| `results/plots/com_distance_convergence.png` | COM distance convergence |
| `results/plots/boundness_qc_min_distance.png` | DOR boundness QC |
| `results/plots/drm_distances/` | Per-mutation DRM sidechain↔DOR distance traces |
| `results/plots/dor_key_contacts/` | Per-mutation crystal-derived DOR contact distances |
| `results/plots/pocket_volume_timeseries/` | Per-mutation NNBP pocket volume traces |
| `results/plots/resistance_heatmap.png` | Mutation × metric heatmap |
| `results/plots/manuscript_global_signatures.png` | Global signature figure |

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

Priority apo systems: WT, F227C, V106A, V106A+P225H, K103N+M230L, A98G+F227C, V106I+F227C.

### 1. Apo system prep (local)

```bash
OPENMM_PLATFORM=CPU python -m src.dor_md_pipeline_apo \
    --mutations WT F227C V106A "V106A+P225H" "K103N+M230L" "A98G+F227C" "V106I+F227C" \
    --holo-runs results/md_runs \
    --apo-runs  results/apo_md_runs \
    --manifest  results/apo_md_manifest.csv
```

This strips DOR (resname 2KW) from each holo minimized PDB, builds an amber-only
system, solvates, and writes `results/apo_md_manifest.csv`.

### 2. Apo MD on Sherlock

```bash
# Sync apo assets to Sherlock
SHERLOCK_USER=rsatija bash scripts/rsync_results.sh push  # or rsync apo_md_runs directly

# Submit (same batched workflow)
bash scripts/sherlock/submit_apo_md_batched.sh 6 12

# Monitor
squeue -u $USER
python3 scripts/sherlock/report_md_progress.py \
    --manifest results/apo_md_manifest.csv --target-ns 10.0 --show-incomplete
```

### 3. Apo analysis

```bash
# Sync completed apo trajectories back
SHERLOCK_USER=rsatija COMPLETE_ONLY=1 MD_PRODUCTION_NS=10.0 \
bash scripts/rsync_results.sh pull  # sync apo_md_runs as well

# Run full apo analysis (PBC fix + tunnel dynamics + DCCM, apo and holo comparison)
bash scripts/run_apo_analysis.sh
```

Key output files:

| File | Description |
|---|---|
| `results/apo_nnbp_tunnel_summary.csv` | Gate mean±std per replicate (apo) |
| `results/holo_nnbp_tunnel_summary.csv` | Gate mean±std per replicate (holo, same mutations) |
| `results/apo_dccm_allosteric_coupling.csv` | NNBP↔fingers/palm/thumb coupling scalars (apo) |
| `results/holo_dccm_allosteric_coupling.csv` | Same, holo |
| `results/plots/apo_nnbp_tunnel/` | Per-gate distance timeseries + distributions (apo) |
| `results/plots/apo_dccm/` | DCCM heatmaps per replicate (apo) |

See `src/README.md` for full argument documentation.

## Notes

- Sherlock GPU allocation: use `salloc -p gpu --gres=gpu:1 --mem=32G` (not `sh_dev -g 1`, which lands on slow MIG-partitioned A30s).
- Local prep runs on CPU (`OPENMM_PLATFORM=CPU`) unless overridden.
- Structural metrics require MDAnalysis (`nnrti-prep` conda env).
