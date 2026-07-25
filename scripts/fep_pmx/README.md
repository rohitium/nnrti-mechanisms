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
# or single leg:
python scripts/fep_pmx/prepare_hybrid.py --leg wt_to_V106A --phase holo --replicate 1

# 4) Y188L apo MD (P0 blocker — assets exist, trajectories missing)
#    Smoke test on interactive GPU first, then batch:
bash scripts/fep_pmx/salloc_apo_gpu.sh
bash scripts/fep_pmx/test_y188l_apo_gpu.sh
bash scripts/fep_pmx/submit_y188l_apo_md.sh
```

Outputs land under `results/analysis/fep_pmx/`.
