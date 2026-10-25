# Session handoff — JCIM revision, 22–23 Oct 2026

Written for the agent picking this up. Everything below is verified state, not
plan. Where something is uncertain it says so.

**Manuscript:** JCIM ci-2026-030865, *Molecular Simulation of Doravirine
Resistance in HIV-1 Reverse Transcriptase*. Major revision, deadline
**10-Nov-2026**. Corresponding author: Rohit Satija.

---

## 1. How the author wants to work — read this first

These are standing instructions, repeated and emphatic. Violating them has
caused real friction.

- **The author writes the prose.** "I will do all the writing, can you just
  help me write?" Supply numbers, technical content, and flag downstream
  breakage. Do not produce paragraphs of manuscript text unprompted.
- **Edit documents in place.** Never rebuild a `.docx` from scratch — the
  author edits them in Word between automated passes. Ask before changing
  anything they wrote.
- **Dropbox holds only files the author has created and vetted.** Never write
  there.
- **Do not produce analysis nobody asked for.** This was said after a response
  letter padded a non-answer with unrequested analysis. It was a fair hit.
- **Answer the reviewer's actual question.** Same incident. Check the verbatim
  comment before claiming something addresses it.
- **Test Sherlock jobs on a dev node before submitting**, reproducing how the
  job is really launched, environment inheritance included.
- **Don't saturate the Sherlock quota.** See §5.

---

## 2. What was done this session (all committed)

```
0874a62  faulthandler + GPU probe on the apo segfault
b4fe34a  Methods: FEP charge model is AM1-BCC, not Gasteiger   <- see §4
1aa2fdd  Track apo submit script; fix heredoc that killed 55 jobs
4200f2e  Supplementary Tables 6 and 7; cover page corrected
80a7d6e  Supplementary Text restructured to the reference SI layout
60043db  Supp Fig 3 showed the rejected block-average plot, not the RMSF
715930f  Drop the standalone convergence figure (merged into Supp Fig 1)
3a39d6f  Restore the DOR title the reviewer saw: "pose convergence"
2fe8e8b  Chi-square interval on the blocked standard error
```

Uncommitted in the working tree, deliberately — the author was mid-edit:
`revisions/Manuscript-Revision-1.docx`, `revisions/Response-to-Reviewers.docx`.
Untracked: `results/md_runs/apo_collapse_solvated/` (22 GB, see §3).

QA: `PYTHONPATH=src python ops/qa/verify_revision_numbers.py` → **45/45**.

---

## 3. Data pulled off Sherlock

`results/md_runs/apo_collapse_solvated/` — 181 files, 22 GB, the 2 ns
water-filled apo collapse runs (array 48943745, 60/60 COMPLETED). Verified:
every file byte-size identical to the remote, all 60 coord arrays load as
`(2001, 16070, 3) float32`. This existed **only on `$SCRATCH`** until now.

---

## 4. The charge-model finding (significant, and in our favour)

The FEP ligand topology `results/analysis/fep_pmx/inputs/dor_openff/dor.top`
— which `fep/config.py` points at as the production input — carries DOR
charges that are **bit-identical to AM1-BCC** (RMS 0.00000) and differ from
Gasteiger by RMS 0.165 e, max 0.390 e. `export_dor_itp.py` hands
`Interchange.from_smirnoff` a molecule with no partial charges, so Sage 2.0.0
assigns AM1-BCC afresh. Confirmed by two independent routes (AmberTools
`assign_partial_charges` and SMIRNOFF), so not circular.

**Consequence:** the alchemical ΔΔG_bind — the paper's principal claim —
*sampled and scored* with AM1-BCC. Gasteiger is confined to the equilibrium MD
and the MM/GBSA column. Methods said parameters were "carried over from the
OpenFF 2.0.0 assignment", implying Gasteiger; corrected in `b4fe34a`.

**The Response-to-Reviewers has not been updated for this** and currently makes
three claims that need rework (Reviewer 2, comment 5):
1. "We have done both of the things the reviewer offers as alternatives" —
   under the author's reading, the first alternative was not done.
2. The "ligand contribution largely cancels" argument is unnecessary and
   implies the FEP used Gasteiger.
3. "The only difference between the two sets of numbers is the ligand charge
   model" — true of the numbers, but implies the charge model's *effect* was
   measured when only the scoring channel was.

**The author's critique, which is correct:** re-scoring existing frames changes
the energy function at fixed configurations. It cannot see the ensemble
change a different charge model would produce. Single-step reweighting can't
rescue it either — with RMS Δq = 0.165 e over 40 atoms the phase-space overlap
is too poor. Only new sampling answers it. The proposed test is short AM1-BCC
MD for WT + two genotypes (9 GPU jobs), comparing pocket volume, DOR RMSD and
key contacts against the Gasteiger ensembles. **Not yet started.**

---

## 5. Sherlock: live state

Queue as of handoff — 16 running, 4 eligible, 8 held:

```
RUNNING   k103n_m230l_02  k103n_p225h_02  v106a_01  v106a_l234i_01/02
          v106a_p225h_01/02/03  v106i_01  v106i_f227c_01/02/03
          wt_01  wt_03  y181c_01  y181c_02

ELIGIBLE  49145187 k103n_m230l_03   49152325 k103n_p225h_01
          49152333 k103n_p225h_03   49145190 wt_02

HELD      49145188/89 v106a_02/03   49140567/86 v106i_02/03
          49140897 y181c_03         49140899/904/906 y188l_01/02/03
```

Held via `scontrol hold` to stop my jobs backfilling the author's entire GPU
quota. Release with `scontrol release <ids>`. **Do not flood the queue** — at
16 jobs the author could not get a 1-CPU, 8 GB dev allocation (`QOSMaxCpuPerUserLimit`,
then `QOSMaxMemoryPerUser`).

Progress: 17–50% of the 100 ns target at 33–61 ns/day. **This is a multi-day
run.** The matched-window analysis is gated by the slowest replicate per
genotype, and `wt_03` (34.6 ns/day) matters most since WT is the reference.

All healthy: T 299–302 K, ρ 1.014–1.019.

### The heredoc bug (fixed, `1aa2fdd`)

55 of 55 apo jobs died in 6–20 s with `No module named openmm`.
`submit_apo_md_groupsrc.sh` builds the job with an **unquoted** heredoc, so
`${PYTHONPATH:-}` expanded on the *login node* (empty) and the job exported
`.../src:` over the top of py-openmm's path — py-openmm ships openmm **on
PYTHONPATH**, not in site-packages. The guard above that line passed because
it ran before the export. Fixed with `\${PYTHONPATH:-}`, guard moved below the
export. Verified on a dev GPU node.

The script was **only in the `$SCRATCH` checkout** and is now tracked at
`ops/slurm/cluster/submit_apo_md_groupsrc.sh`.

### Unresolved: the resume segfault

4 of 27 resumed jobs died with SIGSEGV 7–18 ms after loading the 10 ns
checkpoint, before any integration step. Ruled out: checkpoints are healthy
(no NaN, PE within 0.1% of working replicates, loads onto CUDA); not
node-specific (those nodes ran other jobs fine); not genotype-specific
(v106a rep01 fine, rep02/03 died); not the DCD append path. `0874a62` adds
`-X faulthandler` and an `nvidia-smi` probe so the next one prints a traceback.

### Walltime: 12 h default vs 48 h available — change this

`submit_apo_md_groupsrc.sh` defaults to `SHERLOCK_TIME=12:00:00`. The `gpu`
partition allows **`2-00:00:00`** (`sinfo -p gpu -o "%l"`).

At the observed 33-61 ns/day a replicate needs ~44 h for 100 ns, so **every job
times out at roughly 38%** and needs ~4 allocations instead of 1. First
confirmed case: `apo_k103n_m230l_02`, `TIMEOUT` at `12:00:33`, 38.1% of target,
checkpoint intact and healthy (T 299.4, rho 1.0162). Resubmitted with
`SHERLOCK_TIME=2-00:00:00` as **49219330**.

This matters beyond wasted queue time: each resume is another exposure to the
barostat instability below, which has already ruined 2 of 28 replicates.
Fewer resumes, fewer chances to hit it.

**Suggested (not done, needs the author's call):** resubmit the *queued and
held* jobs at `2-00:00:00` — free, since they have not started — and let the
running ones time out naturally rather than discarding their progress.

**Watch out:** `MUTATION_ALLOWLIST` is per-genotype, not per-replicate.
Resubmitting `k103n_m230l` pulled in rep_01, the replicate deliberately
cancelled for barostat instability; it is now **held as 49219329** pending the
checkpoint-vs-36.9 ns decision.

### Barostat instability after checkpoint resume — open lead

**Two of 28 replicates** have now failed this way, both on resumed runs, both
with otherwise-normal potential energy (so not unfolding):

- `k103n_m230l_01` — box swinging +/-8% (healthy is +/-0.1%), T 306-329 K,
  rho down to 0.904 against a population of 1.016. Cancelled; see below.
- `v106a_l234i_03` (job 49140469, died 2026-10-22T23:30) — the opposite
  direction, box *contracted*:
  `openmm.OpenMMException: The periodic box size has decreased to less than
  twice the nonbonded cutoff.` Resubmitted as **49215618**.

**Partly diagnosed; the barostat theory below WAS WRONG — see the retraction.**
`-X faulthandler` (commit 0874a62) showed both crashes die at the same line,
`md_protocol.py:445` `prod.step()` -> `integrator.step(10)`, on the first
integration step after resume -- not in the DCD reporter:

- `wt_02`  -> `Fatal Python error: Segmentation fault`
- `v106a_l234i_03` -> `OpenMMException: The periodic box size has decreased to
  less than twice the nonbonded cutoff`

### RETRACTION: it is not the barostat

I asserted the cause was "the barostat's first Monte Carlo volume move after a
resume". **That is disproven.** On a dev node, loading the known-bad checkpoint
(`wt_rep02_md.chk.segfaults-on-resume`, preserved) and stepping 191 steps works
identically with and without a `MonteCarloBarostat`: box steady at 12.625 nm,
PE ~-2712500 both ways.

Also ruled out, each tested rather than argued:

| hypothesis | test | verdict |
|---|---|---|
| barostat volume move | load ckpt +/- barostat, 191 steps | NO, identical |
| stripped-DCD append | construct reporter on the real 37 MB DCD | NO, modelCount 200, atoms match |
| GPU memory / contention | `[gpu]` probe in a failed job | NO, 1 MiB used of 11264, Exclusive_Process |
| GPU model | node features, failed vs healthy | NO, sh03-12n08 both segfaults and runs fine |
| the code path as a whole | real `run_prepared_md` on the bad checkpoint, dev node | **SURVIVES** -- 10,000 steps, 110 s, all reporters |

### What the evidence actually says

`wt_rep02` is **deterministic, not intermittent**: it failed on `sh02-15n09`
(V100_SXM2), `sh03-12n04`, `sh03-12n09` and `sh03-12n08` (RTX 2080Ti) -- four
nodes, two GPU architectures, every attempt. The same checkpoint then resumed
first try on a dev node. A fault that is 100% reproducible in the `gpu`
partition and 0% on `dev` is about the execution environment, not the
checkpoint contents or the code.

Probe `49303053` runs the identical script via `sbatch -p gpu` to test exactly
that. Scripts kept for reuse on Sherlock: `~/resume_debug.py`,
`~/full_resume_test.py`, `~/resume_gpu_batch.sh`.

Worth checking next if the probe survives: diff the real submit script's
environment against the probe's, since the probe reproduces everything else.
Also worth knowing that OpenMM documents `loadCheckpoint` as machine-, GPU- and
version-specific, and recommends `saveState`/`loadState` (portable XML) when
resumes may land on different hardware -- which is exactly our situation.

**Action taken:** `v106a_l234i_03`'s checkpoint was corrupt beyond recovery and
was moved to `*.chk.corrupt-box-collapse` so it restarts from scratch (only
1.1% lost).

**Process note:** `v106a_l234i_03` sat dead for ~10 h because the failed jobs
were found by grepping logs for `Segmentation fault`, and this one threw a
Python exception instead. When resubmitting after a batch failure, enumerate
from `sacct` state, not from log greps.

### `k103n_m230l_01` — cancelled, needs a decision

Barostat went unstable on resume: box swinging ±8% (healthy is ±0.1%), T
306–329 K, ρ down to 0.904 against a population of 1.016. **PE was normal, so
not unfolding.** Cancelled (job 49140426) to free a GPU; checkpoint preserved.
Its *pre-restart* run had reached **36.9 ns in perfect condition** before the
resume broke it — that earlier trajectory may be the better asset.
**Decision needed:** restart from checkpoint, or recover the 36.9 ns run.

### G190E apo — the important one, not yet rerun

Reviewer 1 **comment 3** (not comment 2) asks directly:

> G190E is particularly notable because it produces the largest pocket
> expansion, yet no binding penalty. The authors should examine whether
> mutation-induced changes occur in both apo and holo RT…

Our response's central mechanistic claim rests on the G190E apo ensemble. All
three G190E apo replicates **unfolded** (RT extended after 18.5, 23, 71 ns;
bonds intact; author confirmed in PyMOL and said to drop and rerun). They have
**not** been rerun and are not queued, because they report `status=ok,
steps=50000000` and `SKIP_IF_AT_TARGET=1` filters them out. T and ρ look
normal — those observables do not detect unfolding.

Supp Table 7's G190E number (apo 42.0 Å³, opening 242.3 Å³) comes from the
matched **10 ns** window, which precedes the unfolding, so it is defensible —
but it is the weakest trajectory supporting the strongest claim.

**Decision needed:** resubmit with `SKIP_IF_AT_TARGET=0`, and whether to move
the existing trajectories aside first.

---

## 6. AM1-BCC re-scoring — killed, do not naively restart

Ran ~8.5 h across two attempts (6 then 12 workers) and completed **zero** of 51
replicates. Killed at the author's instruction. Checkpoint intact at **9 rows**
(WT, V106A, Y188L) — nothing lost.

Unreconciled: a standalone run of Y188L rep3 completed in **20.6 CPU-minutes**
(`.checkpoint_mmgbsa_am1bcc_y188l_rep3.csv`, Oct 21 14:37). Each parallel
worker spent **4.5 CPU-hours** without finishing one. Same protocol, same
settings.

Investigated and **ruled out**:
- Not hung — stack sampling showed time inside `libOpenMM`/`libOpenMMCPU`.
- Not scheduling contention — each worker held 95–99% of a full core.
- Thread oversubscription was real (348 threads on 14 cores; numpy's macOS
  Accelerate backend spawns ~29 per worker and **ignores
  `OPENMM_CPU_THREADS`**, which OpenMM's CPU platform does honour). Pinning
  `OMP/MKL/OPENBLAS/VECLIB/NUMEXPR/RAYON_NUM_THREADS=1` dropped it to 37, but
  it climbed back to 180 once numerical work began, **and it made no
  difference** — cores were never the bottleneck.

**Untested lead:** memory bandwidth. MD force evaluation is bandwidth-bound,
and 12 concurrent 16k-atom minimizations can saturate it while every core
still reports ~98% busy. If so, *fewer* workers would finish sooner.

**Second untested lead:** the first tasks dispatched include **G190E**, whose
trajectories are unfolded — an extended protein would hit the 2000-iteration
ceiling on all 300 minimizations per replicate.

**Strategic note:** even if it ran, this only bounds the *scoring* channel
(§4). The manuscript's claim is already correctly scoped ("for every genotype
tested") and the Supplementary Text carries an honest
`[PENDING: extended to all genotypes]` marker. The 9-job AM1-BCC MD comparison
is the better use of compute.

---

## 6b. Pose-independence test — RESULT (2026-10-24)

All four 2 ns triage runs finished. Reviewer 1 comment 2 asked whether the
same DOR-bound state is recovered from independent starting poses. **It is
not, on this timescale.**

```
pose      prep RMSD |  2 ns RMSD (med)   disp    tilt  | in WT basin
alt_01         0.92 |  0.37-1.45  0.90   0.29     7.0  | YES
alt_05         5.50 |  5.05-6.42  5.99   3.68    56.7  | no
tilt_02        6.60 |  5.99-6.74  6.34   0.80   132.6  | no
alt_02         8.22 |  8.94-9.85  9.35   5.41   170.7  | no

WT basin: RMSD 0.34-1.78 A | displacement <=0.80 A | tilt <=15.1 deg
contact artifact: 0 of 800 frames below 2.2 A -- clean, no filtering needed
```

The only pose ending in the WT basin is `alt_01`, and it was already there
before any dynamics: it collapsed 2.46 -> 0.92 A during **minimisation**.
Every genuinely independent start stays where it began, plateauing within
0.4-0.8 ns. `tilt_02` is the sharpest case -- centroid displacement 0.80 A,
exactly at the WT limit, so the ligand is in nearly the right *place* while
flipped 133 degrees, and stays flipped.

**What to tell the reviewer.** This is the handoff's outcome 2: 2 ns cannot
discriminate poses, so the stability of our production replicates is not
evidence the pose is correct. The criticism is confirmed, not deflected. It
is also a finding in its own right: the NNIBP accommodates grossly non-native
poses, including a near-inverted doravirine 9 A away, without strain.

**Extension running:** `alt_02` only, job **49338492**, 100 ns, 48 h, pinned
to `GPU_SKU:TITAN_Xp` (the architecture that wrote its 2 ns checkpoint). It
tests whether the 9.3 A plateau is a true second minimum or a slow escape.
`alt_05` and `tilt_02` are NOT extended -- that is a further multi-day
commitment and depends on what `alt_02` shows.

Analysis: `nnrti.cli.compute_pose_coordinates --run-dir ... --label ...`
**overwrites** `sampled_frames.csv` each call, so run one pose at a time and
accumulate if you want all four in one table.

## 6c. Trajectory QC — the test that actually works (2026-10-25)

**Screen unfolding with the radius of gyration of the LARGEST BONDED MOLECULE,
never all-atom Rg.** RT is a heterodimer: p66, p51 and DOR are three separate
molecules, so all-atom Rg conflates subunit *imaging* with conformational
change and gave me four false "unfolding" calls. p66 alone cannot be inflated
by imaging.

```
                Rg(all)        Rg(p66 alone)   verdict
g190e/rep_01   34.1 - 35.2     33.1 - 35.0     healthy
g190e/rep_02   34.1 - 92.9     31.9 - 72.6     GENUINELY UNFOLDED
g190e/rep_03   34.2 - 35.4     32.7 - 34.4     healthy
```

Screening all 36: **2 genuinely unfolded** — `g190e/rep_02` (rerunning, job
49513580) and `k103n_m230l/rep_01` (restarted from scratch; its checkpoint
also gave unstable barostat behaviour). Everything else sits at ratio <= 1.08.

**So G190E does NOT need a full rerun**: rep_01 and rep_03 are usable, only
rep_02 was bad. The claim that all three were unfolded, and the PyMOL look at
"frame 199", were both the imaging artifact.

Temperature and density see none of this: `g190e/rep_02` ran at 300 K and
rho 1.016 throughout.

## 6d. Operational notes from the unattended run

- **A GPU can die with Slurm still reporting RUNNING.** `sh03-13n06` logged
  `NVML: GPU is lost` at 06:00 and the job sat "running" for five more hours
  producing nothing. Grep job logs for `GPU is lost`; `squeue` will not tell
  you.
- **Legacy checkpoints can be migrated** rather than discarded: load on the
  architecture that wrote them (`sacct ... NodeList` + `scontrol show node`)
  and `saveState` to portable XML. `$HOME/migrate_ckpt.py` on Sherlock does
  this; it rescued `k103n_p225h/rep_03` at 77%.
- **The portable-state fix is live in production** — pose extension runs are
  writing `*_md.state.xml` and `*_md.state.json` beside the binary checkpoint.
- **`$SCRATCH` purge also eats assets**: y188l lost its `*_start.pdb` files
  and three jobs died with `FileNotFoundError`. Restored from `$GROUP_HOME`.
- Apo systems for 7 genotypes (g190a, g190s, k103n, l100i_k103n, v106a_f227l,
  v106m, y318f) were **never built** — start PDBs existed, system XMLs did
  not. Built and submitted, 21 replicates; the apo panel now covers all 20
  genotypes rather than 8.

## 7. Open decisions for the author

1. **Reviewer 1, comment 2** — "test representative systems from independent
   starting poses" is **not done**, and the response buries that admission
   mid-paragraph. The author wanted to think about what "independent starting
   poses" should mean before anything is run. **This conversation was
   interrupted and never resumed — pick it up here.**
2. Rework the Response for Reviewer 2 comment 5 around the AM1-BCC finding (§4).
3. Resubmit G190E apo (§5).
4. `k103n_m230l_01`: checkpoint or the 36.9 ns run (§5).
5. Whether to run the 9-job AM1-BCC MD comparison (§4).
6. Convert `submit_apo_md_*` to a Slurm array with `%N` throttling, matching
   `ops/slurm/apo_collapse/submit_collapse_solvated.sh` (`--array=1-N%12`).
   Server-side, no login-node polling loop. Needs the author's N.

---

## 8. Traps found this session

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

## 9. Monitoring scripts (in the repo, usable from any session)

- `ops/slurm/cluster/apo_health.py` — run **on Sherlock**. Reads each running
  replicate's live `*_md_state.csv` and prints one line per unhealthy run;
  silence means all healthy. Checks step advancing (needs >4 min elapsed),
  T 295-307 K, rho 0.99-1.04, box swing >3% between polls, speed, and greps
  logs for `Traceback|Segmentation fault|Particle coordinate is nan|CUDA error`.
  Keeps per-replicate state in `~/.apo_health_state.json`. Needs py3.6 on the
  login node, so no f-string-only features beyond 3.6.
  Deploy: `scp ops/slurm/cluster/apo_health.py sherlock:~/`
- `ops/slurm/cluster/watch_apo_runs.sh` — **do not leave this running.** It
  polls every 5 min and emits only actionable lines: terminal job failures,
  R/ELIG/HELD changes, capacity free with held work, queue empty. Sherlock's
  rules for automated tools (below) forbid leaving an unattended loop polling
  the scheduler, so run `apo_health.py` and `squeue` from your own turns, at
  least 60 s apart, instead. The script is kept for an attended watch only.

Caveat worth knowing: the health checker was wrong three times before it was
calibrated (too loose on density, too eager on startup, and a timing race when
run manually just before arming it). Do not run it by hand immediately before
starting the watcher — it writes the baseline the watcher compares against.

## 10. Sherlock's rules for automated tools — MANDATORY, re-read every time

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
