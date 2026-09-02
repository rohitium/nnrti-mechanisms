# Molecular simulation of doravirine resistance in HIV-1 reverse transcriptase

Code and derived data for the manuscript *Molecular Simulation of Doravirine
Resistance in HIV-1 Reverse Transcriptase* (Satija, Tao & Shafer).

We simulated wild-type HIV-1 RT and 18 doravirine-associated genotypes in
explicit solvent, and asked whether computed changes in DOR binding track
measured phenotypic susceptibility. Two independent estimates were made for
each genotype: a non-equilibrium alchemical free energy (ΔΔG_bind, pmx) and a
wild-type-referenced MM/GBSA interface energy (ΔΔE_Total).

**If you are here to reproduce a number in the paper, read
[REPRODUCE.md](REPRODUCE.md).** It maps every table and figure to the command
that regenerates it, and says honestly which steps you can re-run and which need
a GPU cluster.

## Quickstart

```bash
git clone <repo> && cd nnrti-mechanisms
conda env create -f env/analysis.yml && conda activate nnrti-prep
pip install -e .                      # or: export PYTHONPATH=src

./workflows/05_manuscript_artifacts.sh   # regenerates every table and figure, ~5 min
```

## Layout

| Path | Holds |
|---|---|
| `workflows/` | The five numbered entry points, in the order the study ran. Start here. |
| `src/nnrti/` | The package. `md/` simulation protocol, `fep/` pmx non-equilibrium FEP, `analysis/` analysis library, `cli/` the scripts that build manuscript artifacts, `structure_prep/` system building. |
| `scripts/` | Cluster-side operations: Slurm submission (`sherlock/`, `fep_pmx/`), rsync helpers, MD metadata audit/repair, docx tooling. |
| `data/` | Deposited inputs: structures, the DOR ligand, and `DRM-susceptibilities.csv.xlsx` — the authoritative fold-change source. |
| `manifests/` | Run manifests and provenance logs, including the archive manifest with checksums. |
| `results/` | Manuscript-facing output only. Everything else was archived; see below. |
| `paper/` | Drafts, tables, and [`ARTIFACTS.md`](paper/ARTIFACTS.md) — the authority on what produced each figure and table. |
| `docs/` | Method notes and cluster runbooks. |

`src/nnrti/paths.py` defines the canonical directories. Import them; don't
hardcode directory strings.

## Workflows

| stage | script | where | time |
|---|---|---|---|
| 1 | `workflows/01_prepare_systems.sh` | laptop | minutes |
| 2 | `workflows/02_run_md.sh` | GPU cluster | ~2 weeks |
| 3 | `workflows/03_run_fep.sh` | GPU cluster | ~3 weeks |
| 4 | `workflows/04_analysis.sh` | laptop | ~2 hours |
| 5 | `workflows/05_manuscript_artifacts.sh` | laptop | ~5 minutes |

Stages 4 and 5 run from the deposited data. Stages 2 and 3 are documented but
need a cluster; `REPRODUCE.md` explains what is deposited so you can start at 4.

## What was archived, and where

The repository previously carried about three years of exploratory work
alongside the manuscript analyses, with no way to tell them apart. On
2026-09-01 everything no manuscript artifact depends on was moved to an external
archive: 1,069 result files (5.3 GB) and 76 code files. Nothing was deleted.

- The inventory, with a sha256 for every file, is committed at
  [`manifests/archive_2026-09-01_manifest.csv`](manifests/archive_2026-09-01_manifest.csv).
- The archive itself lives outside the repository and carries its own
  `MANIFEST.csv` and `README.md` explaining each group.
- Every archived path also remains in git history at tag
  `pre-refactor-2026-09-01`, so nothing is unrecoverable.

The selection was made from a static import graph rooted at the scripts that
produce each numbered figure and table, not by hand — which is how it caught
that `results/analysis/modern_md_suite` is load-bearing (it supplies the NNIBP
pocket volume column of Table 3) despite looking exploratory.

## Data-safety conventions

These were learned the hard way; see [docs/repository-layout.md](docs/repository-layout.md) for the full
set.

1. Raw heavy data is not in git and not only on `/scratch`, which is purged
   after 90 days. Trajectories, checkpoints and per-switch data go to durable
   storage with sha256 manifests.
2. git tracks light provenance and code. **Never run `git clean` over untracked
   heavy data** — doing so on the cluster in August 2026 destroyed ~15 GB of raw
   FEP data.
3. Ground truth beats metadata: how much a run actually did is its `state.csv`
   last step, never a JSON claim.
4. Fold-change values come from `data/DRM-susceptibilities.csv.xlsx`. Never
   hardcode them.

## License

See [LICENSE](LICENSE).
