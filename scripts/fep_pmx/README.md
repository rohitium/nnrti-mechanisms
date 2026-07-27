# DOR resistance FEP: pmx + GROMACS NEQ

Full-system alchemical free energy calculations for the **manuscript mutation panel** (19 targets,
19 unique legs).

**Method:** pmx hybrid topology + GROMACS non-equilibrium switching (Crooks/BAR) on rhombic
dodecahedron solvated HIV RT — holo (DOR-bound) and apo (ligand-free) per leg.

**Start here:** [`PLAN.md`](PLAN.md)

**Why not truncated / Perses:** [`APPROACHES.md`](APPROACHES.md)

Perses full-protein pilot (deprecated path): [`../fep_jorgensen/README.md`](../fep_jorgensen/README.md)

## Local pmx setup (Mac)

Hybrid prep and NEQ analysis run locally; GROMACS GPU jobs go to Sherlock.

```bash
bash scripts/fep_pmx/setup_pmx_env.sh
conda activate pmx
export GMXLIB=...   # printed by setup script
pmx -h
```

Do **not** use `pip install pmx` — install [deGrootLab/pmx](https://github.com/deGrootLab/pmx) `develop` from source.

## P0 workflow (local)

```bash
# 1) Inventory what MD assets exist vs missing
python scripts/fep_pmx/asset_manifest.py

# 2) Export DOR OpenFF → GROMACS (needs nnrti-prep env)
conda activate nnrti-prep
python scripts/fep_pmx/export_dor_itp.py

# 3) pmx hybrid structures for P0 legs (all reps × holo/apo)
conda activate pmx
bash scripts/fep_pmx/prepare_p0_hybrids.sh

# 4) GROMACS solvated hybrid systems (Sherlock login — needs gmx + GMXLIB)
source scripts/sherlock/load_gromacs_module.sh
conda activate pmx
REPLICATES=1 bash scripts/fep_pmx/build_p0_systems.sh   # smoke: rep01 only

# 5) Y188L apo MD (P0 blocker — assets exist, trajectories missing)
#    Smoke test on interactive GPU first, then batch:
bash scripts/fep_pmx/salloc_apo_gpu.sh
bash scripts/fep_pmx/test_y188l_apo_gpu.sh
bash scripts/fep_pmx/submit_y188l_apo_md.sh
```

## NEQ workflow (Sherlock GPU)

**Validate on an interactive GPU node first** (catches path/topology bugs before batch submit):

```bash
bash scripts/fep_pmx/salloc_neq_gpu.sh
# inside allocation:
bash scripts/fep_pmx/smoke_neq_em.sh          # V106A holo rep1 EM
LEG=wt_to_V106A PHASE=apo bash scripts/fep_pmx/smoke_neq_em.sh
```

After smoke passes, batch submit (one command — SLURM chains em → equil → extract → switch):

```bash
NEQ_SNAPSHOTS=100 REPLICATES=3 FORCE=1 bash scripts/fep_pmx/prepare_p0_neq.sh
bash scripts/fep_pmx/submit_p0_neq_pipeline.sh
```

Or stage-by-stage (manual dependency):

```bash
STAGE=em      bash scripts/fep_pmx/submit_p0_neq.sh   # normal/CPU (A-state min)
STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh   # gpu: per-λ min → C-rescale warmup → P-R production (PLAN.md §4.3)
STAGE=extract bash scripts/fep_pmx/submit_p0_neq.sh   # normal/CPU
STAGE=switch  bash scripts/fep_pmx/submit_p0_neq.sh   # gpu
# Y188L switches are 500 ps — use STAGE=switch SHERLOCK_TIME=03:00:00 if needed

# 3) Analysis (Mac or Sherlock login with pmx)
conda activate pmx

# 3a) per leg/phase/rep BAR+CGI+Jarzynski (writes analysis.json + work_dist.png)
python scripts/fep_pmx/analyze_neq.py --leg wt_to_V106A --phase holo --replicate 1

# 3b) ΔΔG_bind per genotype (= ΔG_holo − ΔG_apo, summed over legs), mean ± SEM
#     across replicates; correlates vs experiment; runs analyze_neq for any
#     missing leg automatically.
python scripts/fep_pmx/combine_neq.py --targets V106A Y188L --replicates 3

# 3c) QC: Crooks forward/reverse overlap, work outliers, BAR-vs-Jarzynski (§4.6)
python scripts/fep_pmx/qc_neq.py --replicates 3
```

Outputs land under `results/analysis/fep_pmx/`:
- `legs/{leg}/{holo,apo}/rep_*/neq/analysis/` — per-unit `analysis.json`, `results.txt`, `work_dist.png`, `integ_{fwd,rev}.dat`
- `targets/{genotype}/summary.json` — ΔΔG_bind ± SEM + per-replicate table
- `panel_ddg.csv`, `panel_ddg_vs_experiment.png` — ranking vs experiment (Spearman ρ once ≥3 genotypes)
- `panel_qc.csv`, `panel_crooks_overlap.png` — QC table + overlap histograms

**Note:** `pmx analyse` needs `numpy < 2.0` (newer numpy breaks pmx's estimators). The Sherlock `~/.venvs/pmx` env is fine; a local Mac `pmx` conda env may need `pip install 'numpy<2'`.
