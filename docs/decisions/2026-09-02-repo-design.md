# Repository design — proposal, 2026-09-02

**Status: plan, not executed.** Written after yesterday's package refactor, which
fixed the code but left the repository's *shape* untouched.

The question I asked myself: if this repo were handed to a research software
engineer on day one, with the paper in hand and no history, what would they build?

---

## The problem, in numbers

| | now | why it is wrong |
| --- | ---: | --- |
| Root entries | 26 | Should be ~10. `tmp_docx_extracts/`, `inputs/`, `figures/` (one README), `.env`, `.mplconfig` are all visible on `ls`. |
| `manuscript/` files | 112 | 18 superseded drafts, 19 dead figure exports (23 MB) for Figures 4–8 that no longer exist, 15 per-genotype notes, 16 dated working notes. The live submission set is **8 files** — and it lives in a directory named after a colleague. |
| `results/` | 35 GB | 27 GB is superseded: `.bak` analysis DCDs (6.9 GB), OpenMM checkpoints (1.9 GB), system XMLs (6.3 GB), an old in-tree archive (2.6 GB), 253 MB of tables and 315 MB of plots from analyses that were archived yesterday. The manuscript needs **~7.5 GB**. |
| `docs/` | 8 files | 5 of the 8 document the superseded Jorgensen pipeline. |

Three of these are the same failure: **superseded work sits beside current work
with nothing distinguishing them.** Yesterday I fixed that for code and for
`results/analysis/`. It is still true everywhere else.

---

## Design principles

1. **Root is a table of contents.** A directory per concern, nothing else. If
   you can't tell what the project is from `ls`, the root is wrong.
2. **Inputs, code, outputs, paper — four separate trees, never mixed.** Inputs
   are precious and read-only; outputs are regenerable and disposable; the paper
   is a deliverable, not a scratchpad.
3. **One directory is the submission.** Everything a journal receives, together,
   under names a stranger can read. No dates, no colleague names, no `-rev2`.
4. **Working notes are not deliverables.** They are the decision record — worth
   keeping, worth reading, and not worth mixing with the files you upload.
5. **Track what you cannot regenerate.** A 300 MB PNG tree that one command
   rebuilds does not belong in git; a 20 KB derived table that feeds a published
   number does.
6. **Nothing is deleted.** Same procedure as yesterday: `mv` on one filesystem,
   sha256-verified at the destination *before* the deletion is staged, inventory
   committed, git tag as the undo point.

---

## Target tree

```
nnrti-mechanisms/
├── README.md                  ← what this is, 10-line quickstart
├── REPRODUCE.md               ← artifact → command, the submission entry point
├── CITATION.cff               ← new
├── LICENSE  ·  pyproject.toml  ·  .gitignore
│
├── env/                       ← was envs/
│   ├── analysis.yml               (was nnrti-prep.yml)
│   └── fep.yml                    (was nnrti-fep.yml)
│
├── data/                      ← deposited inputs, read-only
│   ├── structures/  ligands/  prepared/
│   └── susceptibility/            DRM-susceptibilities.csv.xlsx, DRMs.csv
│
├── src/nnrti/                 ← unchanged (yesterday's refactor)
├── workflows/                 ← unchanged: 01…05
├── tests/                     ← unchanged: 9 regression tests
├── manifests/                 ← run manifests + archive inventories
│
├── ops/                       ← was scripts/. Cluster operations, not science.
│   ├── slurm/                     submission + .mdp  (was sherlock/, fep_pmx/*.sh)
│   ├── sync/                      rsync_*.sh
│   └── maintenance/               MD metadata audit/repair, docx tooling
│
├── results/                   ← regenerable outputs only, ~7.5 GB
│   ├── md/                        stripped analysis DCDs + topology + JSON
│   ├── fep/                       was analysis/fep_pmx
│   ├── mmgbsa/                    was analysis/binding_energy
│   ├── structure/                 was analysis/{mechanisms,modern_md_suite,md_convergence}
│   └── panel/                     was analysis/{dor_susceptibility_bar_chart,classification_performance}
│
├── paper/                     ← was manuscript/. The deliverable, and only that.
│   ├── ARTIFACTS.md               the artifact authority
│   ├── submission/                exactly what JCIM receives
│   │   ├── manuscript.docx
│   │   ├── supplementary-text.docx
│   │   ├── supplementary-table-{1,2,3,4}.xlsx
│   │   └── figures/               figure-{1,2,3}.pdf, figure-s{1,2}.pdf
│   ├── tables/                    Table-2-energetics.csv, Table-3-structural.csv
│   └── sources/                   Figures.pptx, ACS .dotx template
│
└── docs/
    ├── methods/                   protocol notes
    ├── runbooks/                  cluster runbooks (was ops-side .md)
    └── decisions/                 handoffs, refactor records, the notes below
```

Root goes from 26 entries to 12. `paper/` goes from 112 files to 15.

---

## Phases

Ordered by value per unit of risk. Each is independently revertible, and each
ends with the same gate: **9 tests pass and `workflows/05` regenerates every
artifact unchanged.**

### Phase 1 — root and `docs/` (30 min, near-zero risk)

- `envs/` → `env/`, files renamed to `analysis.yml` / `fep.yml`.
- Delete `figures/` (one README describing a convention we no longer use;
  `results/` and `paper/submission/figures/` cover it).
- Archive `inputs/boltz/` (6 YAMLs for the archived Boltz test) and
  `tmp_docx_extracts/` (10 text dumps of superseded drafts, regenerable).
- Add `.DS_Store`, `.mplconfig`, `.pytest_cache` to `.gitignore` — currently
  ignored only by accident of pattern.
- `docs/` → `methods/`, `runbooks/`, `decisions/`. The five Jorgensen/perses
  documents move to `decisions/` (they explain a real fork in the project) or to
  the archive if you would rather they were gone.
- Add `CITATION.cff`.

**Nothing here is read by any script.** I verified: no code references `figures/`,
`inputs/`, or `tmp_docx_extracts/`.

### Phase 2 — `paper/` (1 hour, low risk, biggest legibility win)

`manuscript/` → `paper/`, and the live submission set is lifted out of
`post-feedback-from-atanu/` into `paper/submission/` under plain names.

Archived (85 files, ~60 MB):

| what | count | note |
| --- | ---: | --- |
| Superseded drafts, Apr–Aug | 18 + 4 | Every one is in git history too. |
| `Figure-1..8.{png,jpg,tiff}` | 19 | 23 MB. Figures 4–8 no longer exist; 1–3 are regenerated by `workflows/05`. |
| `*_resistance_mechanism.md` | 15 | Per-genotype literature notes; superseded by the Discussion. |
| Dated working notes | 16 | CHANGES/EDITS/ERRATA/paragraph drafts → `docs/decisions/`, not the archive: they are the record of why the numbers are what they are. |
| `appendix_a_*.csv`, `Supplementary Table 4.xlsx` (June orphan) | 3 | The orphan is unrelated to the real Supp. Table 4 and cited nowhere. |

Kept in `paper/`: the 8 submission files, `ARTIFACTS.md`, the two derived table
CSVs, `Figures.pptx`, the ACS template, and the two SVG figure sources.

### Phase 3 — `results/` (2 hours, medium risk, biggest disk win: 35 GB → ~7.5 GB)

The rule: keep what `REPRODUCE.md` says is deposited, archive the rest.

| archive | size | why it is safe |
| --- | ---: | --- |
| `md_runs/**/*.bak` | 6.9 GB | Superseded 10 ns analysis DCDs; the 100 ns `.dcd` supersedes them. |
| `md_runs/**/*.xml` | 6.3 GB | Prepared system XMLs — regenerated by `workflows/01`. |
| `md_runs/_archive/` | 2.6 GB | An in-tree archive of an in-tree archive. |
| `md_runs/**/*.chk` | 1.9 GB | OpenMM checkpoints; only needed to resume a run that finished 6 months ago. |
| `results/{plots,tables,visualization}` | 751 MB | Outputs of the analyses archived yesterday. `figure1B_dor_schematic.pdf` is the one exception and moves to `results/panel/`. |
| `results/archive/` | 270 MB | Pre-existing within-analysis archive; belongs with the rest. |
| loose `results/*.csv` | 131 MB | `dor_key_contacts_timeseries_all_mutations.csv` alone is 114 MB, from an archived analysis. |
| `static_mmgbsa_*`, `average_structures/` | 97 MB | Superseded static-structure MM/GBSA. |

Kept: 154 live `*_analysis.dcd` (6.2 GB), their topology PDBs, run JSONs,
`state.csv`, and the seven `results/analysis/` directories, renamed as above.

**This is the phase with real risk**, so it gets the strictest gate: full sha256
inventory before the move, verification after, and `workflows/04` re-run on three
genotypes to prove the trajectories still load.

### Phase 4 — `ops/` and `results/` renaming (1 hour, medium churn, optional)

Rename `scripts/` → `ops/{slurm,sync,maintenance}` and
`results/analysis/{fep_pmx,binding_energy,…}` → `results/{fep,mmgbsa,structure,panel}`.

This is polish, not correctness. It touches ~20 default paths in the CLIs and the
Slurm scripts, and it is the one phase I would happily drop if you want the
churn kept down. **Recommend doing it** — we proved yesterday that path rewrites
are safe here now that the tests exist — but it is the first thing to cut.

---

## What I would not do

- **Rewrite git history.** `.git` is 1.9 GB and would shrink to a few hundred MB
  with a `filter-repo` pass over the archived blobs. I am not proposing it: it
  invalidates every existing clone and every commit hash in the notes, on the eve
  of a submission, to save disk that is not scarce.
- **Delete anything.** Same as yesterday — move, verify, inventory, tag.
- **Touch `data/`.** 45 MB of deposited inputs. It is fine.
- **Split the repo.** A paper repo and a code repo would be cleaner in the
  abstract and worse in practice: the whole point of `REPRODUCE.md` is that the
  paper and the code that made it travel together.

---

## Decisions I need from you

1. **Phase 4 in or out?** Recommend in; it is the only one that is purely
   cosmetic.
2. **Do the Jorgensen/perses method notes stay in `docs/decisions/` or go to the
   archive?** They document a real fork — several weeks of work that the paper
   mentions only as "superseded". I lean toward keeping them: they are 5 files
   and they explain a gap a reader would otherwise notice.
3. **`paper/submission/` filenames.** I propose `manuscript.docx`,
   `supplementary-text.docx`, `supplementary-table-1.xlsx`… Journals usually want
   their own naming at upload time, so this is for your benefit, not theirs.
   Happy to keep `DorDRM-MD-09-02-26.docx` if you'd rather the date stayed
   visible.
4. **Are the 85 archived `paper/` files ever coming back?** If the 18 old drafts
   have no further use, they are also in git history at every commit they were
   touched — the archive copy is belt-and-braces. Say if you want them simply
   left where they are instead.
