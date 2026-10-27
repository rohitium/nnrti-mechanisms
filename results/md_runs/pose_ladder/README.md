# Pose-convergence ladder

Reviewer 1, comment 2: *"test representative systems from independent starting
poses to determine whether the same DOR-bound state is recovered."*

## Why this replaced the earlier test

The first attempt placed DOR 5.5–8.2 Å from the crystallographic pose. Those
runs showed only that a grossly non-native pose is kinetically trapped: the
most displaced was extended to 64 ns and its median pose RMSD moved by 0.06 Å.
That is true but uninformative — nothing could have converged from there, so
the test had no power to detect a basin even if one exists.

The informative experiment is a **ladder**: starts at graded distance from the
crystallographic pose, read for where convergence stops. That radius is the
answer to the reviewer's question.

## Reporting convention

A pose is identified by its RMSD **as generated**, before minimization. What
matters is that a pose at a given distance converged, not whether the minimizer
or the dynamics did the work. (Some poses collapse into the basin during
minimization alone; that still counts as convergence from that start.)

## The wild-type basin

From 3 × 100 ns of crystal-seeded WT, contact-filtered (795 of 860 frames):

| coordinate | min | median | max |
|---|---|---|---|
| pose RMSD (Å) | 0.33 | 0.86 | **1.78** |
| centroid displacement (Å) | 0.04 | 0.28 | **1.75** |
| long-axis tilt (°) | 0.27 | 6.16 | **15.07** |

A run that ends inside that envelope has recovered the crystallographic state.

## The ladder

Twelve new systems at 2 ns, one replicate each, no jitter:

| pose | RMSD as generated (Å) |
|---|---|
| `ladr10_01` | 1.26 |
| `ladr10_02` | 1.32 |
| `ladr15_02` | 1.50 |
| `tilt_04` | 1.72 |
| `ladr15_01` | 1.93 |
| `ladr20_01` | 2.06 |
| `alt_06` | 2.23 |
| `tilt_01` | 2.44 |
| `ladr25_01` | 2.88 |
| `alt_10` | 3.26 |
| `alt_09` | 5.33 |
| `tilt_03` | 5.42 |

Four runs from the earlier attempt extend the ladder, reported the same way:
`alt_01` 2.46 Å (landed), `alt_05` 5.85, `tilt_02` 6.70, `alt_02` 6.81 (did
not). Sixteen points in total, 1.26 to 6.81 Å.

Extend whichever do not land; 2 ns is the triage tier.

## A result from the generation step

Beyond roughly 3 Å RMSD, sterically allowed rigid-body poses barely exist. The
3.0–3.5 and 3.5–4.0 Å bands returned **zero** accepted placements in 300,000
attempts each, even allowing a 4 Å centroid shift, under a 2.2 Å minimum
contact screen. The pocket is tight enough that most non-native placements
simply do not fit. That is worth a sentence in its own right.

## How to read the result

```bash
# per-pose coordinates, for every ladder run
PYTHONPATH=src python -m nnrti.cli.compute_pose_coordinates \
    --run-dir results/md_runs/pose_ladder/pose_<name> --label <name> \
    --output-dir results/analysis/pose_coordinates/tables/<name>

# the figure: WT density, 4NCG as a star, each run as start -> end
PYTHONPATH=src python -m nnrti.cli.plot_pose_convergence
```

`plot_pose_convergence` wants `ladder_summary.csv` with one row per run:
`pose, start_rmsd_angstrom, start_displacement_angstrom, start_tilt_deg,
end_displacement_angstrom, end_tilt_deg, end_rmsd_angstrom`. Start values come
from the pose diagnostics (as generated); end values from the final portion of
the run.

## Provenance

| what | where |
|---|---|
| generation | `src/nnrti/structure_prep/pose_perturbation.py` (`--min-rmsd`/`--max-rmsd` select a band) |
| new pose CIFs | `data/prepared/dor_poses_ladder/` |
| reused pose CIFs | `data/prepared/dor_poses/` |
| system prep | `src/nnrti/structure_prep/prepare_pose_systems.py` |
| manifests | `manifests/pose_ladder_manifest.csv`, `..._b.csv` |
| submission | `ops/slurm/cluster/submit_pose_ladder.sh` |
