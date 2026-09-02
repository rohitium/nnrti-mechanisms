# Reproducing the results in the manuscript

*Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase*

This file maps every table and figure in the paper to the command that produces
it.

---

## Where to start

Simulation (stages 1–3) needs a GPU cluster and several weeks. Analysis (stages
4–5) runs on a workstation from deposited data, and regenerates every number in
the paper.

| Deposited | Location | Supports |
| --- | --- | --- |
| Stripped analysis trajectories, 94 runs, ~7 GB | Zenodo (§Data) | All structural analysis, MM/GBSA, Figures 1 and 3, Supp. Fig. 1 |
| Per-switch work values, `integ_*.dat`, 228 files | This repository, `results/analysis/fep_pmx/legs/**/analysis/` | The free energy panel, Table 2, Figure 2, Supp. Fig. 2 |
| Raw switch data, 5 legs | Durable storage, checksums in `manifests/fep_raw_backup_*.csv` | Re-integrating those legs from the trajectories |
| Susceptibility data | `data/DRM-susceptibilities.csv.xlsx` | Table 1 |

The free energies are recomputed, not cached: the BAR, Crooks and Jarzynski
estimates are derived from the deposited per-switch work values, which is the
step that determines every ΔΔG reported.

Regenerating the trajectories and the alchemical switching themselves is
documented under §Regenerating from scratch.

---

## Setup

```bash
conda env create -f env/analysis.yml && conda activate nnrti-prep
pip install -e .          # or: export PYTHONPATH=src
pytest tests -q           # 9 regression tests over the manuscript artifacts
```

The whole of §Artifacts below is wrapped in one script:

```bash
./workflows/05_manuscript_artifacts.sh
```

and the trajectory analyses that feed it are in `./workflows/04_analysis.sh`.
The per-artifact table is kept because it says which command produces which
number, which the wrapper cannot.

## Manuscript artifacts

| # | Artifact | Command | Cluster? |
| --- | --- | --- | :---: |
| 1 | **Table 1** — DOR susceptibility panel | `python -m nnrti.cli.plot_dor_susceptibility_bars` | no |
| 2 | **Table 2** — ΔΔE and ΔΔG | `python -m nnrti.cli.build_table_2 --docx paper/submission/DorDRM-MD-09-02-26-ACS.docx` | no |
| 3 | **Table 3** — RT–DOR interface observables | `python -m nnrti.cli.build_supplementary_table_4` | no |
| 4 | **Supp. Table 3** — per-replicate ΔΔE and ΔΔG | `python -m nnrti.cli.build_supplementary_table_3` | no |
| 5 | **Supp. Table 4** — per-replicate structural observables | `python -m nnrti.cli.build_supplementary_table_4` (same command as Table 3) | no |
| 6 | **Figure 1** — crystal structure | PyMOL session, not scripted; source `data/structures/` (PDB 4NCG) | no |
| 7 | **Figure 2** — FEP protocol + ΔΔG vs fold | `python -m nnrti.fep.plot_protocol_schematic` and `python -m nnrti.fep.combine_neq --replot-only` | no |
| 8 | **Figure 3** — resistance mechanisms | `python -m nnrti.cli.compute_mechanism_coordinates` then `python -m nnrti.cli.plot_mechanism_panel` | no |
| 9 | **Supp. Figure 1** — MD convergence | `python -m nnrti.cli.compute_md_convergence` then `python -m nnrti.cli.plot_convergence_panel` | no |
| 10 | **Supp. Figure 2** — FEP work distributions | `python -m nnrti.cli.plot_fep_work_distributions` | no |
| 11 | Supp. Table 1 — isolate-level phenotypes | From Stanford HIVDB; not generated here | no |
| 12 | Supp. Table 2 — RT variant sequences | From the prepared systems in `data/prepared/dor_4ncg/` | no |

### Upstream of Tables 2 and Supp. Table 3

Both read `results/analysis/binding_energy/tables/ddg_full.csv` and
`results/analysis/fep_pmx/panel_ddg.csv`. To rebuild the MM/GBSA side from the
deposited trajectories:

```bash
# 1. score 100 evenly spaced frames per replicate, no frame filtering,
#    minimised to convergence. ~25 h on 12 cores.
OPENMM_PLATFORM=CPU OPENMM_CPU_THREADS=1 python -m nnrti.cli.compute_mmgbsa_safe \
  --frame-sampling even --snapshots 100 --discard-fraction 0.25 \
  --snapshot-relaxation unrestrained --relaxation-iterations 2000 \
  --workers 12 --force \
  --output results/.checkpoints/.checkpoint_mmgbsa.csv

# 2. WT-reference, unit-convert and promote to the canonical tables
python -m nnrti.cli.rebuild_binding_energy_sources \
  --mmgbsa-csv results/.checkpoints/.checkpoint_mmgbsa.csv
```

Quality gate (not required, but it is what justifies the sampling window —
it checks DOR stays in one binding mode across the frames scored):

```bash
```

The FEP side is rebuilt from the deposited work values:

```bash
python -m nnrti.fep.combine_neq --targets \
  F227C G190A G190E G190S V106A V106I V106M Y181C Y188L Y318F \
  A98G+F227C V106A+F227L V106A+L234I V106A+P225H V106I+F227C \
  K103N K103N+M230L K103N+P225H L100I+K103N --replicates 3
```

⚠️ `combine_neq --targets X` **rewrites `panel_ddg.csv` with only those
targets** — always pass the full list. Do **not** pass `--force`: it re-analyses
every leg from raw `dgdl.xvg`, which no longer exists for most legs.

---

## Method choices that matter for reproduction

These are decisions, not defaults, and changing them changes the numbers.

- **MM/GBSA samples 100 frames spread evenly across the post-equilibration
  trajectory**, not a terminal window. The terminal-20 protocol used earlier was
  the dominant source of the reported error bars: only ~4% of between-replicate
  variance was frame noise, and adjacent terminal frames are strongly
  autocorrelated. Panel mean SEM 1.07 → 0.59 kcal/mol.
- **No frames are discarded.** Frames with close atomic contacts are repaired by
  minimisation, not excluded — excluding configurations on the basis of their
  energy biases an ensemble average.
- **Minimisation runs to convergence (2000 iterations).** At the old 100-iteration
  cap, overlap frames were unconverged by ~12 kcal/mol; at 2000 they agree with
  clean frames to ~1.3 kcal/mol, which is the genuine difference in pose.
- **WT reference is the mean of the three WT replicates**, not index-matched.
- Full rationale and the measurements behind each:
  `results/analysis/binding_energy/MMGBSA_METHOD_AND_RECOMPUTE.md`.

---

## Regenerating from scratch (needs a GPU cluster)

Neither of these is required to reproduce the published numbers from the
deposited data; they document how that data was made.

**MD** — `ops/slurm/cluster/submit_md_batched.sh` (holo),
`ops/slurm/cluster/submit_apo_md_batched.sh` (apo). 19 genotypes × 3 replicates ×
100 ns, ~38 GPU-h per run. System preparation is
`src/nnrti/structure_prep/preparation.py`, which is idempotent: `--replicates N`
generates rep_04 onwards with fresh jitter seeds and leaves existing replicates
untouched.

**FEP** — pmx non-equilibrium alchemical switching. Prepare hybrids, build
solvated systems, then chain em → equil → extract → switch:
`ops/slurm/fep/`, with `OPERATIONS.md` for the pipeline and
`RUNBOOK_G190E_SEM.md` for a worked campaign including failure recovery.

---

## Data

Heavy artifacts are not in git. `manifests/` records what exists in durable
storage with sha256 checksums:

- `md_runs_manifest_2026-08-28.txt` — 828 MD files
- `fep_raw_backup_2026-08-28.csv` — raw FEP switch data, 4 legs

Verify a restored copy with `sha256sum -c`.
