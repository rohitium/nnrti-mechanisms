# Workflows

Five numbered entry points, in the order the study ran. Each is a thin, readable
shell script: the science lives in the `nnrti` package, and these say only what
runs, in what order, against which inputs.

| stage | script | where it runs | typical wall time |
| --- | --- | --- | --- |
| 1 | `01_prepare_systems.sh` | laptop | minutes |
| 2 | `02_run_md.sh` | GPU cluster (Slurm) | ~2 weeks, 60 jobs |
| 3 | `03_run_fep.sh` | GPU cluster (Slurm) | ~3 weeks, 19 legs x 2 phases x 3 replicates |
| 4 | `04_analysis.sh` | laptop | ~2 hours |
| 5 | `05_manuscript_artifacts.sh` | laptop | ~5 minutes |

**Stages 4 and 5 are the ones a reader can run.** They regenerate every number,
table and figure in the paper from the deposited trajectories and per-switch work
values. Stages 2 and 3 need a cluster and weeks of wall time; see `REPRODUCE.md`
for what is deposited so that you can start at stage 4.

All scripts assume the `nnrti` package is importable — either `pip install -e .`
or `export PYTHONPATH=src`. Set `PYTHON` to override the interpreter.
