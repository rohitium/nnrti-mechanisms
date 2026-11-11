# Repository structure

## Start here

The root README gives setup and the one-command paper build. The submitted
reference is `revisions/submission-ready/`. Run `workflows/05_manuscript_artifacts.sh`
and open `paper/generated/paper/README.md` to find every current figure and table.
[ARTIFACTS.md](../paper/ARTIFACTS.md) explains inputs, methods and curated components;
`paper/artifacts.json` supplies the publication labels used by the collection.

## One role per directory

| Directory | Role |
|---|---|
| `src/nnrti/` | Scientific implementation and CLIs |
| `tests/` | Scientific consistency and reproduction checks |
| `workflows/` | Numbered preparation-to-publication entry points |
| `data/`, `manifests/` | Experimental inputs and simulation provenance |
| `results/analysis/` | Derived numerical sources; preserve stable analysis paths |
| `results/md_runs/` | Working simulation files; large files excluded from Git |
| `paper/sources/`, `paper/tables/` | Editable artwork and reference tables |
| `paper/generated/paper/` | Recreated collection under final publication labels |
| `paper/generated/validation/` | Recreated test logs and verification reports |
| `revisions/submission-ready/` | Frozen submitted documents and workbooks |
| `docs/`, `ops/` | Methods and operational procedures |
| `archive/` | Historical material and ignored local recovery copies |

Final publication numbering belongs in the artifact registry, not in numerical
source paths. Use one shared derived-data source for a figure and its related
table. New manuscript builds write only under `paper/generated/` and never modify
the submitted package. Historical scripts are retained as history, not production entry points.

## Sherlock

The durable project entry point is `$GROUP_HOME/nnrti-mechanisms/`:

```text
code/                     current source checkout
md_runs/                  accessible MD mirror
fep_raw/                  existing FEP archives
archive/                  verified preservation snapshots and retired checkouts
operations/               private inventories, job logs and preservation records
README.md                 entry point and recovery instructions
```

Group Home holds preserved data and code; scratch holds disposable computation
workspaces. Stage inputs to scratch before intensive analysis or simulation.
Scratch files expire 90 days after their last content modification. An archive
is usable only after verification; a directory name or job submission is not
proof of completeness. No scratch source is retired before its preservation
record is complete. Keep deployment inventories in `operations/`, not in public
method documentation. See [storage operations](../ops/STORAGE.md).

The scratch workspace is `$SCRATCH/nnrti-mechanisms-work/`: `results/` and `logs/`
are working outputs; `src/`, `ops/`, `workflows/`, `data/` and `manifests/` link to
the durable checkout. `legacy/` groups historical staging, old code and auxiliary
runs. Set `PROJECT_ROOT` to this scratch workspace before submitting simulations.
The older two project-tree names are compatibility links, not additional checkouts.
