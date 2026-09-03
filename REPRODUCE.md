# Reproducing the results in the manuscript

*Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase*

This file maps every table and figure in the paper to the command that produces
it.

## Setup

```bash
conda env create -f environment.yml && conda activate nnrti-prep
pip install -e .          # or: export PYTHONPATH=src
pytest tests -q
```

## What the repository carries

| In the repository | Supports |
| --- | --- |
| Per-switch work values, `results/analysis/fep_pmx/legs/**/analysis/integ_*.dat` | The free energy panel, Table 2, Figure 2, Supp. Fig. 2 |
| Susceptibility data, `data/DRM-susceptibilities.csv.xlsx` | Table 1 |
| Derived tables under `results/analysis/` | Tables 2 and 3, Supp. Tables 3 and 4 |

The free energies are recomputed rather than cached: the BAR, Crooks and
Jarzynski estimates are derived from the per-switch work values, which is the
step that determines every ΔΔG reported.

Production trajectories (13 GB stripped to protein plus ligand) are not in git.
They are needed only to rebuild the MM/GBSA side and the structural observables,
and are available from the corresponding author.

## Manuscript artifacts

Everything below is wrapped in `./workflows/05_manuscript_artifacts.sh`; the
trajectory analyses that feed it are in `./workflows/04_analysis.sh`. The table
says which command produces which number, which the wrapper cannot.

| # | Artifact | Command |
| --- | --- | --- |
| 1 | **Table 1** — DOR susceptibility panel | `python -m nnrti.cli.plot_dor_susceptibility_bars` |
| 2 | **Table 2** — ΔΔE and ΔΔG | `python -m nnrti.cli.build_table_2` |
| 3 | **Table 3** — RT–DOR interface observables | `python -m nnrti.cli.build_supplementary_table_4` |
| 4 | **Supp. Table 3** — per-replicate ΔΔE and ΔΔG | `python -m nnrti.cli.build_supplementary_table_3` |
| 5 | **Supp. Table 4** — per-replicate structural observables | `python -m nnrti.cli.build_supplementary_table_4` (same command as Table 3) |
| 6 | **Figure 1** — crystal structure | PyMOL session, not scripted; source `data/structures/` (PDB 4NCG) |
| 7 | **Figure 2** — FEP protocol (A–C) and ΔΔG vs fold (D) | `python -m nnrti.fep.plot_protocol_schematic` then `python -m nnrti.cli.plot_panel_by_resistance_category` |
| 8 | **Figure 3** — resistance mechanisms | `python -m nnrti.cli.compute_mechanism_coordinates` then `python -m nnrti.cli.plot_mechanism_panel` |
| 9 | **Supp. Figure 1** — MD convergence | `python -m nnrti.cli.compute_md_convergence` then `python -m nnrti.cli.plot_convergence_panel` |
| 10 | **Supp. Figure 2** — FEP work distributions | `python -m nnrti.cli.plot_fep_work_distributions` |
| 11 | Supp. Table 1 — isolate-level phenotypes | From the Stanford HIV Drug Resistance Database; not generated here |
| 12 | Supp. Table 2 — RT variant sequences | From the prepared systems in `data/prepared/dor_4ncg/` |

## Rebuilding Table 2 and Supp. Table 3

Both read `results/analysis/binding_energy/tables/ddg_full.csv` and
`results/analysis/fep_pmx/panel_ddg.csv`.

MM/GBSA, from the production trajectories (~25 h on 12 cores):

```bash
OPENMM_PLATFORM=CPU OPENMM_CPU_THREADS=1 python -m nnrti.cli.compute_mmgbsa_safe \
  --frame-sampling even --snapshots 100 --discard-fraction 0.25 \
  --snapshot-relaxation unrestrained --relaxation-iterations 2000 \
  --workers 12 --force \
  --output results/.checkpoints/.checkpoint_mmgbsa.csv

python -m nnrti.cli.rebuild_binding_energy_sources \
  --mmgbsa-csv results/.checkpoints/.checkpoint_mmgbsa.csv
```

Free energies, from the per-switch work values:

```bash
python -m nnrti.fep.combine_neq --targets \
  F227C G190A G190E G190S V106A V106I V106M Y181C Y188L Y318F \
  A98G+F227C V106A+F227L V106A+L234I V106A+P225H V106I+F227C \
  K103N K103N+M230L K103N+P225H L100I+K103N --replicates 3
```

`combine_neq --targets X` rewrites `panel_ddg.csv` with only those targets, so
pass the full list. `--force` re-analyses each leg from raw `dgdl.xvg`, which is
not distributed.

## Parameters that determine the numbers

- MM/GBSA scores 100 frames spaced evenly across the post-equilibration
  trajectory (final 75 ns of each 100 ns replicate).
- No frames are excluded. Close atomic contacts are resolved by minimization,
  since excluding configurations by energy biases an ensemble average.
- Minimization runs to convergence, 2000 iterations.
- The WT reference is the mean of the three WT replicates, not index-matched.
- Structural quantities are averaged within a replicate before being averaged
  across replicates, so uncertainties are between-replicate.

## Regenerating the simulations (GPU cluster, several weeks)

Not required to reproduce the published numbers.

**MD** — `ops/slurm/cluster/submit_md_batched.sh` (holo) and
`submit_apo_md_batched.sh` (apo). 19 genotypes × 3 replicates × 100 ns, roughly
38 GPU-hours per run. System preparation is
`src/nnrti/structure_prep/preparation.py`; it is idempotent, so `--replicates N`
adds replicates with fresh jitter seeds and leaves existing ones untouched.

**FEP** — pmx non-equilibrium switching; see `ops/slurm/fep/README.md`.

## Checksums

`manifests/` records what exists in durable storage. Verify a restored copy with
`sha256sum -c`.
