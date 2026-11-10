# Workflows

Run from the repository root with the `nnrti-prep` environment active.

| Script | Purpose |
|---|---|
| `01_prepare_systems.sh` | Prepare structures and parameterized inputs |
| `02_run_md.sh` | Equilibrium production simulations |
| `03_run_fep.sh` | Alchemical FEP workflow |
| `04_analysis.sh` | Core trajectory analysis; revision-specific apo/frame analyses also have individual CLIs |
| `05_manuscript_artifacts.sh` | Export the frozen November submission and rebuild supported current plots/tables into `paper/generated/` |

Inspect actual filenames in this directory before invoking stages 1–4; those stages require raw data, simulation dependencies and substantial compute. They are not prerequisites for the derived-data manuscript build. [REPRODUCE.md](../REPRODUCE.md) gives the recovery and validation procedure.
