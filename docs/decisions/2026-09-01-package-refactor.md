# Refactor executed — 2026-09-01

`docs/REFACTOR_PLAN.md` was the proposal; this is what actually happened, what
differed from the plan, and what you need to do on the cluster before running
anything there again.

Pushed to `main` as four commits, `b8dc3c9..71de407`. The pre-refactor state is
tagged `pre-refactor-2026-09-01` (also pushed), so the whole thing is one
`git checkout` away from being undone.

---

## ⚠️ Do this before the next cluster job

The package moved. `python -m src.md...` no longer exists, and `python -m
nnrti.md...` needs `src` on the path.

```bash
ssh <cluster>
cd nnrti-mechanisms
git pull                       # applies the renames
export PYTHONPATH="$PWD/src"   # or: pip install -e .
python -m nnrti.md.sherlock.run_md_job --help    # should print usage
```

Every submit script under `ops/slurm/cluster/` and `ops/slurm/fep/` now exports
`PYTHONPATH` itself, so submitted jobs are fine once the pull is done. It is the
*interactive* shell that needs the export.

---

## What changed

### 1. Exploratory work left the working tree

1,069 result files (5.3 GB) and 145 code files moved to
`~/Career/00_Github/_nnrti_archive/2026-09-01_pre-submission/`.

**Nothing was deleted.** Every file was moved with `mv` on the same filesystem,
then verified byte-for-byte by sha256 against an inventory taken beforehand —
1,214 files, 0 missing, 0 corrupt. The inventory is committed at
`manifests/archive_2026-09-01_manifest.csv` and duplicated in the archive as
`MANIFEST.csv`; the archive also carries its own `README.md` explaining each
group. Every archived path also remains in git history at the tag.

`results/analysis` went from 5.3 GB to 158 MB and now contains exactly the seven
directories a manuscript artifact reads.

### 2. `src/` became the `nnrti` package

```
src/nnrti/
  analysis/   library + the new panel.py
  cli/        artifact-producing scripts only  (87 -> 16)
  fep/        pmx non-equilibrium FEP  (was ops/slurm/fep/*.py)
  md/         simulation protocol + cluster workers
  structure_prep/  utils/  paths.py
```

`ops/` keeps only cluster-side operations: Slurm submission, rsync, MD
metadata audit/repair, docx tooling. `pyproject.toml` makes `pip install -e .`
work.

### 3. `workflows/` — five numbered entry points

`01_prepare_systems` · `02_run_md` · `03_run_fep` · `04_analysis` ·
`05_manuscript_artifacts`. Stages 4 and 5 run on a laptop from deposited data;
2 and 3 point at the cluster runbooks.

### 4. Tests that test the paper

The only two tests in the repo covered the superseded Jorgensen pipeline, and
broke on import once it was archived. Replaced with nine regression tests, each
aimed at a defect this project actually shipped:

- Table 2's components must sum to its total — the half-regenerated ΔΔE_SA
  column went unnoticed for three weeks.
- Supplementary Table 3's FEP sheet must reproduce `panel_ddg.csv` — it had
  drifted in *sign* for G190E.
- Table 3 must equal Supplementary Table 4's summary.
- Figure 2D must not import an exploratory plotting CLI.

```bash
PYTHONPATH=src pytest tests -q     # 9 passed
```

---

## Where the plan was wrong

**`results/analysis/modern_md_suite` is load-bearing.** The plan listed it as
archivable. It produces `pocket_volume_per_rep.csv`, which is the V(NNIBP)
column of the new Table 3 — a dependency that only came into existence on
2026-09-01, after the plan was written. Kept, along with
`compute_modern_md_suite.py`.

**Figure 2D depended on two "exploratory" scripts.** It imported the
Susceptible/Resistant/Uncertain genotype sets from
`plot_wt_referenced_occupancy_tick_lines`, which imports a helper from
`plot_triplet_contact_story`. Archiving either would have broken Figure 2D
silently. Those definitions are now `nnrti.analysis.panel`, and a test enforces
that the figure never re-acquires such a dependency.

Both were found by building a static import graph from the artifact-producing
scripts rather than by reading directory names. Neither would have survived a
by-hand pass.

**Option (a), not (b).** The plan recommended keeping `src/` as the package root.
You chose the rename; it is done, and every call site, runbook, shell script and
doc was rewritten in the same commit.

---

## Verification

| check | result |
| --- | --- |
| archive integrity | 1,214 files verified by sha256; 0 missing, 0 corrupt |
| package imports | 83 modules import cleanly |
| tests | 9 passed |
| shell scripts | all parse (`bash -n`) |
| `workflows/05` end-to-end | ran clean, 11 steps, no tracebacks |
| Table 3 | **byte-identical** after regeneration |
| Supp. Tables 3 and 4 | sheet-for-sheet identical (only zip timestamps differ) |
| Table 2 | all numbers identical; see below |

**One intentional difference.** Regenerating Table 2 through `workflows/05`
changes genotype labels from `A98G + F227C` to `A98G+F227C`, because the
workflow does not pass `--docx`. Every number is unchanged. The unspaced form
matches Tables 1 and 3, so I left it — but if you want Table 2 to keep the
spaced spelling, run it with
`--docx paper/submission/DorDRM-MD-09-02-26.docx` as
`REPRODUCE.md` documents.

---

## What this did not do

- **`.git` is still 1.9 GB** (1.56 GiB packed). Moving files out of the working
  tree does not shrink history; only a history rewrite would, and I did not do
  that. The 40 GB working tree is the win: `results/` is now 35 GB, essentially
  all of it the `md_runs/` trajectories.
- **`results/md_runs` (34 GB), `results/plots`, `results/tables`,
  `results/visualization`, `results/archive` and the loose `results/*.csv` files
  were left alone.** Some are certainly stale, but they were ambiguous enough
  that archiving them without you awake was the wrong trade. A second pass could
  reclaim several GB.
- **Figure 1A** is still an unscripted molecular-graphics session, and
  Supplementary Tables 1 and 2 are still external exports. Both are noted as
  such in `REPRODUCE.md`.
