# Project operating notes

Start with README.md, REPRODUCE.md, docs/REPOSITORY.md and paper/ARTIFACTS.md. Run workflows/05_manuscript_artifacts.sh to build the numbered paper collection and test reports. The frozen reference is the local JCIM package submitted 10 November 2026 in revisions/submission-ready/. Never overwrite it during regeneration; use paper/generated/. Abstract and author-reviewed manuscript wording are not cleanup targets.

## Repository safety

- Never stage recursively, run git clean, reset --hard or discard unrelated work. Inspect status and stage only an explicit reviewed file list. Two machines may push main: fetch before pushing; never force push.
- Raw MD/FEP files and local backups are ignored, not disposable. GitHub protects code and selected derived data, not trajectories.
- One root archive/ holds historical documents, audits and local recovery copies. Do not create another archive directory. Historical notes/scripts describe earlier states and are not current editing instructions.
- The user is the sole commit author; do not add co-author attribution.

## Analysis

Use the local nnrti-prep conda environment with PYTHONNOUSERSITE=1 and src/nnrti, not the retired src.analysis namespace. Do not modify or rerun holo simulations during manuscript regeneration. Authoritative phenotype data: data/DRM-susceptibilities.csv.xlsx. Authoritative MM/GBSA input: results/analysis/binding_energy/tables/ddg_full.csv. Reproduce the current submitted results before reporting regenerated replacements.

MD run JSON status and DCD timestamps are unreliable. Cross-check maximum state.csv Step against restart records; at 2 fs, 50 million steps is 100 ns. Resumes can append duplicate frames and nonmonotonic step records. PBC-corrected sidecars may be stale after extension: check frame counts before reuse. Match apo/holo analysis windows (the old 10 ns apo default is inappropriate for 100 ns comparisons). Preserve current G190E apo pocket/volume filtering and circular χ1 handling.

## Environments

Read docs/ENVIRONMENTS.md. Preparation, analysis and paper builds run locally in nnrti-prep; Sherlock simulations use py-openmm/8.1.1_py312 or GROMACS site modules, not Conda. Legacy pmx is a separate optional local compatibility environment. Do not remove other environments without confirming their users and preserving provenance.

## Sherlock data

Read ops/STORAGE.md before recovery or transfer. Check both $SCRATCH/nnrti-mechanisms-git and $SCRATCH/nnrti-mechanisms; jobs may reference either. Determine live job paths from squeue/logs. Do not delete a supposedly obsolete tree or copy/archive a run while jobs are writing to it. Follow the cluster's /etc/agents/AGENTS.md when working there.

Durable code: $GROUP_HOME/nnrti-mechanisms/code (the old checkout path is a compatibility link). Durable data: $GROUP_HOME/nnrti-mechanisms, including dated fep_raw archives and MD mirrors. Scratch purges after 90 days; no Oak allocation. Backup completeness must be checked from current manifests, not assumed from these historical notes. Large cluster copies run as batch jobs.
