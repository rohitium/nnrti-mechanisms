# FEP on Sherlock: WT → V106A pilot and batch

Split compute between **local** (Perses hybrid prep + analysis) and **Sherlock GPU**
(lambda-window MD). The solvated Perses hybrid is ~200k atoms; do not run production
windows on a laptop CPU.

## Environments

| Location | Runtime | Used for |
| --- | --- | --- |
| Local Mac | `nnrti-prep` (conda) | `prepare --backend perses`, `panel`, `analyze` |
| Sherlock | `module load chemistry py-openmm/8.1.1_py312` | `worker.py` lambda windows (CUDA) |

Sherlock uses the same OpenMM module stack as `submit_md_batched.sh` — no conda env required.
The worker only needs OpenMM (stdlib otherwise); MBAR analysis runs locally.

One-time local setup:

```bash
bash scripts/fep_jorgensen/setup_perses_env.sh
```

---

## Phase 1 — Local: hybrid preparation

From the repo root on your Mac:

```bash
cd /path/to/nnrti-mechanisms
git checkout jorgensen-fep
git pull

# Perses hybrid for WT -> V106A (~15-20 min CPU)
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.prepare \
  --mutation V106A \
  --backend perses

# Worker manifest for this leg only (11 lambda states)
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.panel \
  --mutation V106A
```

Expected outputs:

```text
results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/hybrid_system.xml
results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/hybrid_topology.pdb
results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/schedule.json
results/analysis/fep_jorgensen/worker_manifest_v106a.csv
```

Verify prepare:

```bash
ls -lh results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/
head results/analysis/fep_jorgensen/worker_manifest_v106a.csv
```

---

## Phase 2 — Sync prepared leg to Sherlock

Replace `sherlock` with your SSH alias if different.

```bash
# From local repo root
rsync -av --progress \
  results/analysis/fep_jorgensen/legs/wt_to_V106A/ \
  sherlock:$SCRATCH/nnrti-mechanisms/results/analysis/fep_jorgensen/legs/wt_to_V106A/

rsync -av \
  results/analysis/fep_jorgensen/worker_manifest_v106a.csv \
  sherlock:$SCRATCH/nnrti-mechanisms/results/analysis/fep_jorgensen/
```

On Sherlock, confirm:

```bash
ssh sherlock
cd $SCRATCH/nnrti-mechanisms
ls -lh results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/
wc -l results/analysis/fep_jorgensen/worker_manifest_v106a.csv   # expect 12 lines (header + 11 states)
```

---

## Phase 3 — Interactive GPU pilot (recommended)

Request a GPU shell and run **one short λ window** before the batch array.

```bash
ssh sherlock
cd $SCRATCH/nnrti-mechanisms
git pull   # pick up scripts/sherlock/run_fep_jorgensen_pilot.sh if needed

bash scripts/sherlock/salloc_fep_jorgensen_gpu.sh
```

When the interactive shell starts:

```bash
cd $SCRATCH/nnrti-mechanisms-git   # or your repo clone
export PROJECT_ROOT=$PWD

# Optional: confirm CUDA + OpenMM (pilot script loads the module for you)
module load chemistry py-openmm/8.1.1_py312
nvidia-smi
python3 -c "from openmm import Platform; print(Platform.getPlatformByName('CUDA').getName())"

# One λ=0 window, ~50 ps production (pilot defaults)
bash scripts/sherlock/run_fep_jorgensen_pilot.sh
```

Pilot succeeded if you see:

```text
results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows/state_00_energies.csv
```

Optional: test another state or longer pilot:

```bash
FEP_STATE_INDEX=5 FEP_PROD_STEPS=100000 bash scripts/sherlock/run_fep_jorgensen_pilot.sh
```

Exit the allocation when done: `exit`

---

## Phase 4 — Batch: all 11 lambda windows

Still on Sherlock (login node is fine for `sbatch`):

```bash
cd $SCRATCH/nnrti-mechanisms
./scripts/sherlock/submit_fep_jorgensen_v106a.sh
```

Monitor:

```bash
squeue -u $USER
tail -f logs/fep_jorgensen.<JOBID>_0.out
```

Default production per window (from `FEPConfig`):

| Setting | Value |
| --- | --- |
| λ states | 11 |
| Equilibration | 250,000 steps (500 ps) |
| Production | 2,500,000 steps (5 ns) |
| Energy interval | 2,500 steps |
| Samples / window | 1,000 |
| Platform | CUDA |

Adjust walltime if needed:

```bash
SHERLOCK_TIME=48:00:00 SHERLOCK_MEM=32G \
  ./scripts/sherlock/submit_fep_jorgensen_v106a.sh
```

---

## Phase 5 — Sync results back and analyze locally

When all array tasks finish:

```bash
# On Sherlock — quick check
ls results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows/state_*_energies.csv | wc -l
# expect 11
```

```bash
# From local Mac
rsync -av --progress \
  sherlock:$SCRATCH/nnrti-mechanisms/results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows/ \
  results/analysis/fep_jorgensen/legs/wt_to_V106A/holo/windows/
```

Analyze:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze \
  --leg-dir results/analysis/fep_jorgensen/legs/wt_to_V106A

# or by target label
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --target V106A
```

Summary lands in:

```text
results/analysis/fep_jorgensen/legs/wt_to_V106A/summary.json
results/analysis/fep_jorgensen/targets/V106A/summary.json
```

---

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `Missing prepared holo artifact` on Sherlock | Re-run Phase 2 rsync |
| `CUDA` platform not found in pilot | Confirm `salloc --gres=gpu:1`; `nvidia-smi` |
| Pilot hangs / OOM | `SHERLOCK_MEM=64G`; reduce `FEP_PROD_STEPS` first |
| `Unsupported lambda parameter functions` | `schedule.json` must say `perses-default` (Perses prep) |
| Batch job 0 samples | Read `logs/fep_jorgensen.*.err`; often path or missing `module load py-openmm` |

---

## Full manuscript panel (later)

```bash
# Local: prepare all 19 legs (hours total)
bash results/analysis/fep_jorgensen/prepare_all.sh

# Local: full manifest
PYTHONPATH=. python -m scripts.fep_jorgensen.panel

# Rsync entire legs/ tree, then on Sherlock:
./scripts/sherlock/submit_fep_jorgensen_windows.sh
```
