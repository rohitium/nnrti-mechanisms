# Pose-independence test — prepared systems, ready to submit

Four wild-type RT–DOR systems that differ **only** in the DOR starting pose. Prepared
locally on 23 Oct 2026; **nothing has been submitted to Sherlock** (the queue was
saturated at the time — check before adding to it).

## Why these exist

Reviewer 1, comment 2: *"test representative systems from independent starting poses to
determine whether the same DOR-bound state is recovered."* The submitted replicates were
all seeded from the 4NCG pose with a 0.1 Å jitter, so they are not independent in the
sense asked. These four are.

## Provenance

| step | what | where |
|---|---|---|
| pose generation | random rigid placements of DOR in the 4NCG pocket, screened for sterics | `src/nnrti/structure_prep/pose_perturbation.py` |
| pose coordinates | pocket-frame displacement / tilt / RMSD | `src/nnrti/cli/compute_pose_coordinates.py` |
| selection | farthest-point sampling on pairwise heavy-atom RMSD | 5 picked, `alt_08` dropped as redundant with `tilt_02` |
| system prep | minimise → solvate → system XML | `src/nnrti/structure_prep/prepare_pose_systems.py` |

Pose CIFs: `data/prepared/dor_poses/*.cif`. Rotation-family poses were set aside in
`data/prepared/dor_poses/rotations/` — they are sterically equivalent to the alternative
placements but carry zero centroid displacement, so they add nothing.

## The four systems

Pose RMSD is DOR heavy-atom RMSD from the crystallographic pose, in the pocket frame
(NNIBP Cα superposition), with no ligand superposition.

| pose | RMSD as generated | RMSD after prep | tilt after prep | note |
|---|---|---|---|---|
| `alt_01` | 2.46 Å | **0.92 Å** | 6.8° | collapsed into the WT basin during minimisation — kept deliberately, the recovery is the result |
| `alt_02` | 6.81 Å | 8.22 Å | 165.4° | near end-for-end flip, the most independent start |
| `alt_05` | 5.85 Å | 5.50 Å | 99.3° | |
| `tilt_02` | 6.70 Å | 6.60 Å | 130.0° | low displacement (0.80 Å), high tilt |

**Reference for judging the outcome — the WT sampled basin** (3 × 100 ns, contact-filtered):
pose RMSD **0.34–1.78 Å**, centroid displacement **≤ 0.80 Å**, long-axis tilt **≤ 15.1°**.
A start that ends inside that envelope has recovered the crystallographic state.

## What to run

2 ns, one replicate each — a triage tier to see which poses survive, not a convergence
run. Manifest: `manifests/pose_test_manifest.csv` (4 tasks, all fields absolute).

Per task, following `submit_md_batched.sh`:

```bash
python -m nnrti.md.sherlock.run_md_job \
  --mutation pose_alt_02 --replicate 1 --production-ns 2 \
  --system-xml    results/md_runs/pose_test/pose_alt_02/rep_01/assets/pose_alt_02_md_rep01_system.xml \
  --topology-pdb  results/md_runs/pose_test/pose_alt_02/rep_01/assets/pose_alt_02_md_rep01_start.pdb \
  --minimized-pdb results/md_runs/pose_test/pose_alt_02/rep_01/pose_alt_02_minimized_rep01.pdb \
  --output-json   results/md_runs/pose_test/pose_alt_02/rep_01/pose_alt_02_rep01.json
```

Quote the sbatch heredoc and escape `\${PYTHONPATH:-}` — an unquoted heredoc expands it
on the login node and erases py-openmm's path, which killed 55 apo jobs.

## How to read the result

```bash
PYTHONPATH=src python -m nnrti.cli.compute_pose_coordinates \
    --run-dir results/md_runs/pose_test/pose_alt_02 --label alt_02
PYTHONPATH=src python -m nnrti.cli.plot_pose_landscape
```

Three outcomes, all publishable:

- **Converges into the WT basin** — the crystallographic pose is recovered from a start
  that did not contain it. This is the answer the reviewer asked for.
- **Stays put** — 2 ns cannot discriminate poses, which is the reviewer's criticism
  confirmed rather than deflected. Say so.
- **Goes somewhere else** — a second accessible minimum. Worth looking at closely before
  reporting; check it is not the contact artifact below.

2 ns tests whether a pose *survives*, not whether it *converges*. Any convergence claim
needs the longer runs behind it.

## Caveats

- **Jitter is off** (`jitter_angstrom=0.0`). The panel used 0.1 Å to decorrelate
  replicates seeded from one pose; here it would blur the coordinate under test.
- **Contact QC.** 31 of 60 panel replicates have minimum ligand–protein heavy-atom
  contacts below 2.2 Å and all are flagged `OK`, because `boundness_qc` only tests the
  6 Å upper bound. In WT rep_02 this affects 36% of frames and inflated the apparent WT
  basin from 1.78 Å to 2.60 Å. `compute_pose_coordinates` records
  `min_contact_angstrom` per frame and `plot_pose_landscape --min-contact` filters on it;
  apply the same filter to these runs. The wider question of which published numbers are
  affected is open.
- **The WT reference is three replicates** at unequal DCD stride (500 / 180 / 180
  frames); the plot subsamples to the shortest so one run does not dominate.
