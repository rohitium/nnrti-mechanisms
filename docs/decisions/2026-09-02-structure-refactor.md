# Structure refactor executed — 2026-09-02

`docs/decisions/2026-09-02-repo-design.md` was the proposal. This is what
happened, including the two places I changed the plan once the data disagreed
with it.

Tag `pre-structure-2026-09-02` is the undo point.

## Result

| | before | after |
| --- | ---: | ---: |
| root entries | 26 | 16 |
| `paper/` (was `manuscript/`) | 112 files | 18 |
| `results/` | 35 GB | 13 GB |
| `docs/` | 8 flat files | `methods/`, `runbooks/`, `decisions/` |
| `scripts/` | 50 files, mixed purpose | `ops/{slurm,sync,maintenance}` |

1,165 files (~22 GB) moved to `_nnrti_archive/2026-09-02_structure/`, each
sha256-verified at the destination before its deletion was staged. Combined with
2026-09-01: **2,379 files, ~27 GB archived, nothing deleted.**

## Where the plan was wrong

**Phase 4's `results/` rename was dropped.** The plan proposed flattening
`results/analysis/fep_pmx` → `results/fep` and so on. Then I checked what is
actually inside those trees: **912 files under `analysis/fep_pmx/` and
`analysis/binding_energy/` record their own paths** in run configs and
provenance JSONs written at execution time, and `md_runs/` is named by 60
manifest rows and 81 run JSONs. Renaming would either falsify those records or
leave them pointing at paths that no longer exist. Removing one redundant path
component is not worth a provenance tree that lies about where it ran. The
`scripts/` → `ops/` half of Phase 4 was done in full; `results/README.md`
documents the decision.

**The prepared system XMLs were archived, then restored.** 6.3 GB, the single
largest line item. The manifest and the FEP hybrid preparation reference them by
path. `workflows/01` regenerates them in principle, but that is not verified
bit-for-bit, and breaking the documented stage-2 path to reclaim disk on a
volume with 853 GB free is the wrong trade. Verifying prep reproducibility, then
archiving them, is a clean follow-up.

## What the gate caught

**Two run topology PDBs were symlinks into `results/visualization/`.** Archiving
that tree left WT replicate 1 and Y188L replicate 1 with dangling topologies —
the trajectories were untouched, so nothing looked wrong until the gate tried to
load all 60 runs. Both are now real files. A symlink between result trees turns
any directory move into silent data loss; `results/README.md` says so.

**`build_table_2`'s `--docx` default pointed at a June draft.** It had been
silently supplying the row order and genotype labels for Table 2 ever since;
archiving that draft is what surfaced it. Now points at the current submission
draft.

## Gate, run after every phase

| check | result |
| --- | --- |
| all 60 manifest runs resolve to a trajectory and topology | 60/60 |
| broken symlinks under `results/` | 0 |
| archived files verify by sha256 | 1,165/1,165 |
| package imports | 83 modules |
| tests | 9 passed |
| `workflows/05` end to end | clean, no tracebacks |
| shell scripts parse | all |

## Follow-ups

- Verify `workflows/01` reproduces the prepared system XMLs bit-for-bit, then
  archive them (6.3 GB).
- `.git` is still 1.9 GB. Only a history rewrite would shrink it, which
  invalidates every clone; not worth doing before submission.
- `results/md_runs/**/*.pdb` is 3.4 GB across `start`, `topology` and `final`
  PDBs. `final.pdb` is referenced by the run JSONs but read by nothing; worth a
  look once someone can confirm it is not needed for a restart.
