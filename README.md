# Doravirine resistance in HIV-1 reverse transcriptase

Simulation code and derived data for *Molecular Simulation of Doravirine
Resistance in HIV-1 Reverse Transcriptase* (Satija, Tao & Shafer).

Wild-type HIV-1 RT and 18 doravirine-associated genotypes were simulated in
explicit solvent, three 100 ns replicates each. For every genotype the effect of
the substitution on doravirine binding was estimated two ways — a
non-equilibrium alchemical free energy (ΔΔG_bind, pmx) and a
wild-type-referenced MM/GBSA interface energy (ΔΔE_Total) — and compared with
measured phenotypic susceptibility.

## Install

```bash
conda env create -f environment.yml && conda activate nnrti-prep
pip install -e .          # or: export PYTHONPATH=src
pytest tests -q
```

## Quickstart

```bash
./workflows/05_manuscript_artifacts.sh    # every table and figure, ~5 min
```

To go further back, `workflows/04_analysis.sh` rebuilds the derived tables from
the trajectories (~2 h). Stages 1–3 build the systems and run the simulations
and need a GPU cluster; see [REPRODUCE.md](REPRODUCE.md), which maps every
numbered table and figure to the command that produces it.

## Layout

| Path | Contents |
|---|---|
| `workflows/` | Five numbered entry points, in the order the study ran |
| `src/nnrti/` | The package: `md/` simulation protocol, `fep/` alchemical free energy, `analysis/` analysis library, `cli/` artifact scripts, `structure_prep/` system building |
| `ops/` | Cluster operations: Slurm submission and runbooks, sync helpers, maintenance utilities |
| `data/` | Inputs: structures, ligand, prepared systems, susceptibility data |
| `results/` | Generated output — see [results/README.md](results/README.md) |
| `paper/` | The manuscript and its artifacts — see [paper/ARTIFACTS.md](paper/ARTIFACTS.md) |
| `manifests/` | Run manifests and archive inventories |
| `docs/` | Methods notes and the project's decision record |
| `environment.yml` | Conda specification (one file, all steps) |

Directory constants are in `src/nnrti/paths.py`.

## Data

Trajectories and per-switch free energy data are deposited separately;
[REPRODUCE.md](REPRODUCE.md) lists what each supports. Exploratory work not used
in the manuscript is held outside the repository, inventoried with checksums in
`manifests/archive_*.csv`.

## Citing

See [CITATION.cff](CITATION.cff). Licensed under the MIT License.
