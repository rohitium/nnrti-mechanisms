# Technical notes — pitfalls and data issues

Distilled from the session handoffs. These are the things that have cost time
or produced wrong numbers. Everything here is durable; the session narratives
they came from are in `_archive/`.

## Traps

- **`$SCRATCH` purge** has eaten Python sources, venv `bin/`, git objects, and
  `submit_apo_md_groupsrc.sh`. Use `$GROUP_HOME/nnrti-repo`.
- **Submitted figure artwork lives in the submitted `.docx`**, not in
  `results/`. The PNGs were regenerated post-submission and a panel title
  drifted from what the reviewer saw ("pose convergence" →
  "internal conformational stability"). Extract `word/media/` instead.
- **`SKIP_IF_AT_TARGET` cannot tell "finished" from "finished badly."** The one
  genotype most needing a rerun is the one it filters out.
- **Piping `module load` runs it in a subshell** — the environment change is
  discarded and `module list` then shows nothing.
- **openrsync** (macOS default) rejects `--info=progress2`. Use
  `rsync -rltvz --partial`.
- **Login-node `python3` is 3.6.8** — no `capture_output`, no `text=`.
- **Health-check calibration:** thresholds set from textbook expectations were
  useless. The population is T = 300.3 ± 0.3 and ρ = 1.016 ± 0.002; a window of
  0.95–1.08 hid a 6% error for two hours. Box-volume swing (>3% between polls)
  is the most sensitive barostat signal.

## Sherlock rules for automated tools — re-read after every compaction

Sherlock installs binding instructions for AI agents at
`/etc/agents/AGENTS.md`, with topic files beside it (`slurm.md`, `storage.md`,
`software.md`, `policy.md`). The SSH banner points at it on every login. It
applies whether the tool runs on a Sherlock node **or drives Sherlock remotely
over SSH**, so it governs this workflow from the Mac too.

**Rule 12 requires reading it at the start of each Sherlock session and again
after any summary or context compaction.** An earlier read does not count, and
neither does this note or a memory file. Keep this paragraph in any future
handoff: preserving the requirement is itself part of the rule.

The rules that bite this project in particular:

- **No unattended polling loops** (rule 10). Status checks at least 60 s apart;
  prefer `sbatch --dependency=afterok:<id>` or `sbatch --wait` over polling.
  This is why `watch_apo_runs.sh` must not be left running.
- **Never resubmit a failed job without finding the cause**, and never submit a
  duplicate because one is pending. `QOSMaxGRESPerUser` means wait, not retry.
- **`--exclude`/`-x` is unsupported**; select hardware with `-C` from
  `sh_node_feat`. (Pinning the GPU is still usually the wrong fix — §5.)
- **Nothing may exist only on `$SCRATCH`** (90-day purge), and artificially
  refreshing files to dodge the purge is an account-termination offence.
- **A faulty node is reported to the author**, who decides about contacting
  srcc-support@stanford.edu. Stop jobs with `scancel`, never by killing
  processes on compute nodes.
- **Scope `find`/`du` to our own directories**; never scan filesystem roots.

## Impossible ligand contacts in the stripped DCDs (2026-10-25)

**This is the most serious thing found in this session and it is unresolved.**

Measuring the minimum ligand-protein heavy-atom distance on *every* frame of
all 60 holo replicates -- not the 50-frame sample `boundness_qc` takes:

```
replicates with any frame < 2.2 A : 38 of 60
frames < 2.2 A                    : 1705 of 14769  (11.5%)
worst contact                     : 0.74 A  (G190S rep_01, 2KW F3x <-> VAL186 backbone C)
worst replicates                  : Y188L rep_01 67.8% of frames, G190E rep_02 61.7%,
                                    G190S rep_01 58.7%, V106I_F227C rep_02 52.0%
clean contrast                    : F227C rep_01, 0 frames, median min contact 2.87 A
```

Table: `results/analysis/contact_qc_per_frame.csv`.

**What is established.** The contacts are real in the stored coordinates.
Internal geometry of both molecules is correct (protein N-CA median 1.47 A,
ligand bonds 1.38 A), so there is no atom-order mismatch. Frames are correctly
paired: scanning ligand-against-protein offsets -3..+3, offset 0 is optimal
for both a clean and a pathological replicate. The PBC correction is not the
cause -- per-frame minima are identical to 1e-6 A with and without it. And the
MD's own potential energy is flat across the entire run (sd ~1700 kJ/mol, not
one excursion beyond 10 sd in 13,719 samples), which a 0.74 A F...C contact
cannot be: it is worth order 1e5-1e6 kJ/mol of Lennard-Jones energy. The
stored geometry is one the integrator appears never to have propagated.

**The mechanism is NOT identified.** Two hypotheses were tested and both fail;
do not re-propose either:

- *Periodic imaging.* For every frame checked the stored image is already the
  one placing the ligand nearest the pocket, and none of the 27 lattice shifts
  improves the contact.
- *Ligand-specific corruption, via `md_final.pdb`.* The contact does look
  better in `md_final.pdb` (2.95 vs 1.89 A for G190S rep_01), but superposing
  on the protein shows the **protein** differs by as much as the ligand
  (2.77 vs 3.07 A; Y188L rep_01 2.47 vs 1.69 A). The DCD's last frame is
  simply a different snapshot in those runs. The clean control matches
  `md_final.pdb` to 0.001 A on both, which is what makes the comparison
  interpretable -- and what shows it proves nothing here.

**Scope -- narrower than it first looks, because this was already handled
where it matters most.** MM/GBSA is *not* affected: `mmgbsa.py` defaults to
`snapshot_relaxation="unrestrained"` and runs a capped all-atom minimisation
on every snapshot before evaluating the energy, precisely to relieve steric
overlap (its docstring notes that the older hydrogen-only relaxation could
not). That is the pipeline that 1/r^12 repulsion would otherwise have
poisoned, and it is protected.

What reads the frames as stored, with no relaxation, is the geometric side in
`result_collector.py`: contacts, H-bonds, pocket volume, pose RMSD. Even there
do not overstate it -- a contact count with a 4 A cutoff still counts an
over-close atom once, so the effect is distorted geometry rather than a wrong
count, and whether any published number actually moves is something to
measure rather than assert.

This is almost certainly the same phenomenon as the single case already known
in V106I+F227C rep_01 (TYR185:O <-> F3x at 1.87 A) -- the same ligand atom,
F3x, recurs -- but at 38 replicates rather than one.

**Done:** `compute_boundness_qc` now has a lower bound and emits a `CLASH`
flag (commit 9e57d1f). It previously tested only the 6 A upper bound and
called every one of these `OK`. It still samples 50 frames at stride 20, so
it under-counts; the full measurement is the CSV.

**Needs the author's decision:** whether the geometric metrics are worth
recomputing with clashing frames excluded -- measure the difference first,
it may be negligible -- and whether the mechanism is worth chasing into the
reporter before the resubmission. The energies need nothing.

## compare_apo_ensembles was reporting nonsense (fixed, 2026-10-25)

The drift table read 826 A^3 for G190E against a ~100 A^3 scale. Three causes:

- `total_ns` came from the run JSON, which is only written at completion, so
  an in-flight replicate advertises its previous target. WT rep_02 and rep_03
  claimed 10 ns while holding 91 and 83 ns of trajectory.
- Uniform frame spacing was assumed (`dt = total_ns / n_frames`). An extension
  appends at the new run's stride, so a 380-frame replicate is 200 frames at
  0.05 ns followed by 180 at 0.5 ns, and its "10 ns matched window" was really
  1.9 ns. The two-stride model reproduces the elapsed time of all 44 apo
  replicates exactly.
- G190E rep_02, whose pocket has come apart (pocket Ca Rg 18.8 A, median
  volume 1766 A^3), was averaged into its genotype.

Frame times now follow the real schedule, the trajectory is read from the
`*_analysis_pbcfix.dcd` sidecar where one exists, and two classes of replicate
are excluded *and named in the output*: destroyed pockets, and runs shorter
than the matched window (a 3.6 ns run had been entering the table as a 10 ns
observation, because every frame it owns is inside the window). Drift fell to
the -48..+87 A^3 range.

**Still open:** the same uniform-spacing assumption lives in the three profile
workers in `src/nnrti/analysis/result_collector.py`
(`time_ps = frame_idx * production_ps / (n_total_frames - 1)`, three sites).
20 of 60 holo replicates have frame counts inconsistent with a single stride,
and two show a direct displacement changepoint (`L100I_K103N/rep_01` at frame
200, `V106A_P225H/rep_03` at frame 360). The state CSVs record no restarts, so
the true holo schedule is not recoverable from the files on disk. Means are
unaffected; time axes and any time-windowed statistic are, including the
MM/GBSA even-frame sampling, which is even in index and therefore weighted
toward the densely sampled early segment.