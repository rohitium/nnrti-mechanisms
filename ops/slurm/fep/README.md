# Alchemical free energy calculations

pmx hybrid topologies with GROMACS non-equilibrium switching, for the manuscript
genotype panel: 19 unique legs, each run in a DOR-bound (holo) and a ligand-free
(apo) box, three replicates per leg.

ΔΔG_bind = ΔG_holo − ΔG_apo. Multi-substitution genotypes are decomposed into
sequential single-residue legs, summed, with standard errors combined in
quadrature.

## Environment

Hybrid preparation and analysis run locally; GROMACS jobs run on the cluster.

```bash
bash ops/slurm/fep/setup_pmx_env.sh
conda activate pmx
export GMXLIB=...   # printed by the setup script
```

pmx is installed from [deGrootLab/pmx](https://github.com/deGrootLab/pmx); the
`pmx` package on PyPI is a different project. Its estimators require
`numpy < 2`.

## Preparation

```bash
python -m nnrti.fep.asset_manifest                    # what MD assets exist
conda activate nnrti-prep
python -m nnrti.fep.export_dor_itp                    # DOR OpenFF -> GROMACS
conda activate pmx
bash ops/slurm/fep/prepare_p0_hybrids.sh              # pmx hybrid structures
source ops/slurm/cluster/load_gromacs_module.sh
bash ops/slurm/fep/build_p0_systems.sh                # solvated hybrid systems
```

## Switching

SLURM chains the four stages: em → equil → extract → switch.

```bash
NEQ_SNAPSHOTS=100 REPLICATES=3 FORCE=1 bash ops/slurm/fep/prepare_p0_neq.sh
bash ops/slurm/fep/submit_p0_neq_pipeline.sh
```

Stages can also be submitted individually with
`STAGE={em,equil,extract,switch} bash ops/slurm/fep/submit_p0_neq.sh`.

Each end state is sampled for 5 ns; 100 evenly spaced frames per end state each
seed one switch, giving 100 forward and 100 reverse work values per phase per
replicate. Switches are 100 ps, or 500 ps for legs whose forward and reverse
work distributions are widely separated, which includes the charge-changing
WT→K103N leg inherited by every K103N-containing genotype.

## Analysis

```bash
conda activate pmx
python -m nnrti.fep.analyze_neq --leg wt_to_V106A --phase holo --replicate 1
python -m nnrti.fep.combine_neq --targets V106A Y188L --replicates 3
python -m nnrti.fep.qc_neq --replicates 3
```

`analyze_neq` gives BAR, Crooks Gaussian intersection and Jarzynski estimates per
leg, phase and replicate. `combine_neq` assembles ΔΔG_bind per genotype as a mean
± SEM across replicates, and runs `analyze_neq` for any leg that is missing.
`--replot-only` rebuilds the figures from the existing CSV.

Outputs land under `results/analysis/fep_pmx/`:

| Path | Content |
| --- | --- |
| `legs/{leg}/{holo,apo}/rep_*/neq/analysis/` | `analysis.json`, `results.txt`, `work_dist.png`, `integ_{fwd,rev}.dat` |
| `targets/{genotype}/summary.json` | ΔΔG_bind ± SEM and the per-replicate table |
| `panel_ddg.csv` | The panel, and the source for Table 2 |
| `panel_qc.csv`, `panel_crooks_overlap.png` | Forward/reverse overlap and estimator agreement |

## Syncing

```bash
SHERLOCK_USER=<user> bash ops/sync/rsync_fep_pmx.sh pull
```

Pulls analysis output without trajectories.
