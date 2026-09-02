# Repository refactor for manuscript submission

**Status: proposal, not executed.** Drafted 2026-08-28 while the local MM/GBSA
rescore and the Sherlock G190E FEP campaign were both running.

Goal: a reader of the JCIM paper should be able to clone this repo and follow a
single documented path from the deposited inputs to every number, table and
figure in the manuscript — without tripping over three years of exploratory work.
Nothing is deleted; superseded material moves to a dated archive with checksums.

---

## What is wrong now

**1. There is no path from data to manuscript.** `scripts/run_analysis.sh` is the
nominal entry point, but of the 82 CLIs in `src/nnrti/cli/`, **61 are
referenced by no shell script and no doc** — and that unreferenced set includes
every script that builds a manuscript artifact:

| artifact | built by | in `run_analysis.sh`? |
| --- | --- | :---: |
| Table 1 | `plot_dor_susceptibility_bars` | no |
| Table 2 | `build_table_2` | no |
| Figure 2 | `src/nnrti/fep/plot_protocol_schematic.py`, `combine_neq.py` | no |
| Figure 3 | `compute_mechanism_coordinates` → `plot_mechanism_panel` | no |
| Supp. Figure 1 | `compute_md_convergence` → `plot_convergence_panel` | no |
| Supp. Figure 2 | `plot_fep_work_distributions` | no |
| Supp. Table 3 | `build_supplementary_table_3` | no |
| MM/GBSA inputs | `compute_mmgbsa_safe` → `rebuild_binding_energy_sources` | partly |

`run_analysis.sh` instead runs contact-occupancy, pocket-volume, DRM-distance and
tunnel analyses that the manuscript never cites.

**2. Exploratory work is indistinguishable from manuscript work.** 949 of 2416
tracked files under `results/` (~5.2 GB on disk) belong to lines of inquiry that
appear nowhere in the paper. Verified by searching the 09-02 draft: `logistic`,
`tunnel`, `DCCM`, `occupancy`, `pocket volume`, `sidechain deletion`, `medoid`,
`regression` all return **zero** mentions; `Boltz` and `Jorgensen` appear **only
in the bibliography**, not as analyses.

| directory | tracked | size | manuscript? |
| --- | ---: | ---: | --- |
| `results/analysis/fep_pmx` | 1397 | — | **yes** — Table 2, Fig 2, Supp Fig 2, Supp Table 3 |
| `results/analysis/binding_energy` | 50 | — | **yes** — Table 2, Supp Table 3 |
| `results/analysis/md_convergence` | 13 | — | **yes** — Supp Fig 1 |
| `results/analysis/mechanisms` | 4 | — | **yes** — Fig 3 |
| `results/analysis/dor_susceptibility_bar_chart` | 3 | — | **yes** — Table 1 |
| `results/analysis/fep_jorgensen` | 119 | 5.0 G | no — superseded by pmx NEQ |
| `results/boltz` | 294 | 144 M | no — separate Boltz-2 test |
| `results/analysis/new_logistic_regression` | 259 | 27 M | no |
| `results/analysis/triplet_story_analyses` | 107 | 51 M | no |
| `results/analysis/modern_md_suite` | 59 | 14 M | no |
| `results/analysis/custom_mechanism_*` (4 dirs) | 87 | 4.7 M | no |
| `results/analysis/openmm_sidechain_deletion_*` (4) | 11 | 44 K | no |
| `results/analysis/occupancy_stats` | 6 | 324 K | no |
| `results/analysis/ligand_pocket_features` | 2 | 1.8 M | no |
| `results/analysis/wt_original_vs_sherlock_rerun` | 5 | 20 K | no |

**3. Scripts are split across two trees with no rule.** `scripts/` holds both
SLURM submission (`scripts/sherlock/`, `scripts/fep_pmx/`) and one-off diagnostics
(`diagnose_mmgbsa_*.py`, `audit_mmgbsa_protocol_first5.py`) at top level;
`src/nnrti/cli/` holds 82 more. Which tree a thing lives in is historical.

---

## Safety: what may move, and when

The binding constraint is that **the Sherlock checkout syncs by `git pull`, and
running array elements resolve their paths at runtime.** A pull that relocated a
manifest, a task-id file or `scripts/fep_pmx/*` under a live job would break it
the same way the task-id collision did on 2026-08-15.

Therefore the refactor is staged, and **nothing that a running job reads moves
until its queue is empty.**

| stage | when | touches | safe because |
| --- | --- | --- | --- |
| **1** | now | additive only: `REPRODUCE.md`, archive tooling, this plan | no file moves at all |
| **2** | after local MM/GBSA rescore finishes | archive the exploratory `results/` dirs above | the FEP pipeline reads none of them; `results/analysis/fep_pmx/**` is untouched |
| **3** | after Sherlock queues drain | move code (`scripts/` ↔ `src/`), rewrite entry points | nothing is executing against those paths |

Stage 2 note: a Sherlock `git pull` during stage 2 applies the renames, so the
files still exist at their new paths. It is still safest to defer the pull until
the campaign finishes.

---

## Target layout

```
├── README.md               # what this is; 10-line quickstart
├── REPRODUCE.md            # THE submission entry point: artifact -> command
├── docs/                   # method notes, runbooks, this plan
├── env/                    # conda specs (renamed from envs/)
├── data/                   # deposited inputs: structures, ligand, susceptibilities
├── manifests/              # run manifests + archive log
├── src/nnrti/              # importable package (was src/)
│   ├── md/                 # simulation protocol + Sherlock workers
│   ├── fep/                # pmx NEQ (was scripts/fep_pmx/*.py)
│   ├── analysis/           # analysis library
│   └── cli/                # ONLY CLIs that produce a manuscript artifact
├── workflows/              # shell entry points, one per stage
│   ├── 01_prepare_systems.sh
│   ├── 02_run_md.sh              (Sherlock)
│   ├── 03_run_fep.sh             (Sherlock)
│   ├── 04_analysis.sh            (local)
│   └── 05_manuscript_artifacts.sh
├── results/                # manuscript-facing outputs only
├── paper/
└── archive/2026-08-28_pre-submission/
    ├── MANIFEST.csv        # original path, archived path, sha256, bytes, reason
    ├── README.md           # what each subdirectory was, and why it was set aside
    ├── boltz/
    ├── fep_jorgensen/
    ├── logistic_models/
    ├── triplet_story/
    ├── modern_md_suite/
    ├── sidechain_deletion_fep/
    └── code/               # the 61 unreferenced CLIs, frozen
```

`archive/` at top level rather than inside `results/` because it will hold code as
well as results. The existing `results/archive/<date>_<reason>/` convention (7
entries, e.g. `2026-08-28_binding_energy_terminal20`) stays as-is for
within-analysis supersessions — it is working and is referenced from the method
notes.

### On moving code (stage 3)

`src/nnrti/analysis/...` → `src/nnrti/analysis/...` breaks every `python -m
nnrti.cli.X` invocation, including ones written into the FEP runbooks and
into `RUNBOOK_G190E_SEM.md`. Two options:

- **(a) Rename the package**, update every call site and runbook in the same
  commit. Cleaner for a submitted artifact; one atomic breaking change.
- **(b) Leave `src/` as the package root**, and only prune `src/nnrti/cli/`
  down to the manuscript scripts. Much lower risk, still fixes the "which script
  made Figure 3" problem.

**(b) is recommended** unless you specifically want the package importable as
`nnrti`. The reproducibility win comes from `REPRODUCE.md` and the pruning, not
from the package name.

---

## REPRODUCE.md contract

One row per manuscript artifact, each with the command that regenerates it and
the inputs it consumes. Every command must run from a clean clone plus the
deposited data. Where a step needs a GPU cluster (MD, FEP), the row says so and
points at the runbook, and a cached intermediate is committed so downstream
steps still run.

Draft rows:

| artifact | command | needs cluster |
| --- | --- | :---: |
| Table 1 | `python -m nnrti.cli.plot_dor_susceptibility_bars` | no |
| Table 2 | `python -m nnrti.cli.build_table_2 --docx <draft>` | no |
| Supp. Table 3 | `python -m nnrti.cli.build_supplementary_table_3` | no |
| Figure 2 | `python -m nnrti.fep.plot_protocol_schematic` + `combine_neq.py --replot-only` | no |
| Figure 3 | `python -m nnrti.cli.compute_mechanism_coordinates && ... plot_mechanism_panel` | no |
| Supp. Figure 1 | `python -m nnrti.cli.compute_md_convergence && ... plot_convergence_panel` | no |
| Supp. Figure 2 | `python -m nnrti.cli.plot_fep_work_distributions` | no |
| MM/GBSA table | `python -m nnrti.cli.compute_mmgbsa_safe --frame-sampling even --snapshots 100` | no (GPU helps) |
| MD trajectories | `workflows/02_run_md.sh` | **yes** |
| FEP ΔΔG panel | `workflows/03_run_fep.sh` | **yes** |

---

## Open questions for the author

1. **Package rename** — option (a) or (b) above?
2. **Deposit scope** — do the 100 ns trajectories go to Zenodo, or only the
   analysis-ready stripped DCDs? This decides whether `REPRODUCE.md` starts from
   trajectories or from derived tables.
3. **Archive in-repo or out?** Keeping `archive/` in git preserves provenance but
   keeps ~5 GB of superseded results in history. Moving it to `$GROUP_HOME` with
   a checksum manifest committed instead would shrink the clone substantially —
   the `.git` directory is already 1.8 GB.
