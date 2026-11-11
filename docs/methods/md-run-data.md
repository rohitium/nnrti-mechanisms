# MD run data

## File layout

Holo replicates use `results/md_runs/<genotype>/rep_NN/`; apo replicates use
`results/md_runs/apo/<genotype>/rep_NN/`.

| File | Purpose |
|---|---|
| `*_analysis.dcd` | Protein–ligand analysis trajectory (protein only for apo) |
| `*_analysis_topology.pdb` | Matching trajectory topology |
| `*_md_state.csv` | Recorded integration steps and energies |
| `*.json` run record | Run configuration and artifact paths |
| `*_md.state.xml`, `*_md.state.json` | Portable restart state, when retained |
| `*_md.chk` | Binary checkpoint; compatibility depends on the runtime |
| `*_minimized_repNN.pdb`, `assets/` | Minimized structure and prepared system |

## Duration and trajectory validation

Production uses a 2 fs timestep: 100 ns corresponds to 50,000,000 integration
steps. Analysis output is nominally sampled every 500 ps.

Verify duration using the maximum `Step` in `*_md_state.csv`, cross-checked
against available restart records. Run-summary JSON fields and the presence or
absence of restart files are insufficient on their own to establish completion.

Resumed runs can append duplicate frames and nonmonotonic step records. Derive
physical time from integration-step records, and reconcile resume boundaries
before selecting analysis windows. Frame count alone does not establish duration.
Before using PBC-corrected trajectories, verify that they correspond to the
current source trajectory and cover the intended analysis window.

## Data access and preservation

Raw trajectories and large simulation intermediates are excluded from Git.
The repository contains selected derived data and provenance; reproducing
trajectory-level analyses requires the corresponding raw files.

See [reproduction instructions](../../REPRODUCE.md) for analysis requirements
and [storage operations](../../ops/STORAGE.md) for transfers and archive verification.
