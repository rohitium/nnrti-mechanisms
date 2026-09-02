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
| `ops/` | Cluster operations, no science: `slurm/` submission and runbooks, `sync/` rsync helpers, `maintenance/` MD metadata repair and Word tooling. |
| `data/` | Deposited inputs: structures, the DOR ligand, and `DRM-susceptibilities.csv.xlsx` — the authoritative fold-change source. |
| `results/` | Generated output, all of it rebuilt by `workflows/`. See `results/README.md`. |
| `paper/` | The deliverable: `submission/` is exactly what the journal receives, `tables/` the derived CSVs, `sources/` figure assembly, and `ARTIFACTS.md` the authority on what produced each figure and table. |
| `manifests/` | Run manifests and archive inventories with checksums. |
| `docs/` | `methods/`, `runbooks/`, and `decisions/` — the dated record of why things are as they are. |
| `env/` | Conda specs: `analysis.yml` and `fep.yml`. |
| `tests/` | Nine regression tests over the manuscript artifacts. |

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

The repository carried three years of exploratory work alongside the manuscript
analyses with nothing distinguishing them. Over 2026-09-01 and 09-02 everything
no manuscript artifact depends on was moved to an external archive: **2,379
files, about 27 GB.** Nothing was deleted.

- Inventories, with a sha256 for every file, are committed at
  `manifests/archive_2026-09-01_manifest.csv` and
  `manifests/archive_2026-09-02_manifest.csv`.
- The archive carries its own `MANIFEST.csv` and `README.md` per group.
- Every archived path also remains in git history, at tags
  `pre-refactor-2026-09-01` and `pre-structure-2026-09-02`.

What stayed was chosen mechanically: a static import graph rooted at the scripts
that produce each numbered figure and table, plus the data those scripts read.
That is how it caught `results/analysis/modern_md_suite` — which looks
exploratory and supplies the NNIBP pocket volume column of Table 3.

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
