# pmx NEQ FEP — operations & troubleshooting runbook

How the Sherlock pipeline is wired, how to monitor it, and how to recover from the
failures we have actually hit. Read this before touching a running panel.

- **How to *run* a panel:** [`README.md`](README.md)
- **Scientific plan / stage design:** [`PLAN.md`](PLAN.md), [`docs/pmx-neq-fep-plan.md`](../../docs/pmx-neq-fep-plan.md)
- **Current live state (what's running *right now*):** [`STATUS.md`](STATUS.md)

Claude/agents cannot reach Sherlock. The human runs commands there and pastes output.
Every command block below is for the human to run on a Sherlock **login node**.

---

## 0. The one gotcha that bites every session

The login-node default `python3` is **3.6** and dies on the repo's `from __future__ import
annotations` with `SyntaxError: future feature annotations is not defined`. **Always** activate
the pmx env first (gives Python 3.9 + `GMXLIB` pointing at pmx's mutant force field):

```bash
cd /scratch/users/rsatija/nnrti-mechanisms-git
source scripts/sherlock/activate_pmx_env.sh   # prints [pmx] python: .../3.9.0 and GMXLIB=...
```

If you see the `future feature annotations` SyntaxError, you forgot this step.

---

## 1. Mental model (four facts that explain everything)

1. **Stages run in a fixed chain:** `em → equil → extract → switch`.
   - `em`, `extract` → **normal (CPU)** partition. `equil`, `switch` → **gpu** partition.
   - `submit_p0_neq_pipeline.sh` submits all four with SLURM `afterok` dependencies in one shot.
   - `submit_p0_neq.sh` submits **one** stage; env vars `STAGE`, `MANIFEST`, `DEPENDENCY` drive it.

2. **Dependencies are whole-array `afterok`.** extract depends on the *entire* equil array,
   switch on the *entire* extract array. **Consequence:** if even one equil element FAILS or
   TIMES OUT, `afterok` can never be satisfied and **all** downstream elements hang forever as
   `DependencyNeverSatisfied`. One bad element stalls the whole batch. This is the single most
   important failure mode — see §4.

3. **Every stage is idempotent.** `run_neq_task.py` skips any unit that already has its outputs
   (`em.gro`, `equil.gro`+`equil.trr`, extract frames, switch `dgdl.xvg`). So **re-submitting a
   whole stage safely re-runs only the missing/failed units** — completed ones return in seconds.
   You almost never need to target individual array indices; just resubmit the stage.
   - Nuance: within `equil`, the per-λ `min` and `warmup` sub-steps resume if partially done, but
     the 5 ns production restarts from `warmup.gro` (no `-cpi` resume). A killed production redoes
     its 5 ns — fine for correctness.

4. **GPU QOS cap:** `MaxSubmitPU=100` array elements. A 3-leg batch = equil 36 + switch 36 = 72
   GPU elements ≤ 100 (extract is CPU, doesn't count). This is why panels are batched in 3-leg
   groups. Don't submit a 4th leg's GPU stages while a batch's equil+switch are live.

Array-index ↔ task mapping: the stage submitter writes all task-ids for a stage to
`results/analysis/fep_pmx/neq_<stage>_task_ids_chunk000.txt`; **SLURM array index `i` = line `i+1`**
of that file. That's how a failed `..._<JOBID>_<IDX>` maps back to a leg/phase/rep.

---

## 2. Monitoring

```bash
source scripts/sherlock/activate_pmx_env.sh
# summary counts per stage:
python -m nnrti.fep.audit_neq_panel --manifest <manifest.csv> | grep '==='
# which units are incomplete and why:
python -m nnrti.fep.audit_neq_panel --manifest <manifest.csv> | grep -E 'FAIL|RUNNING'
```

- `RUNNING` = a unit mid-flight (equil checkpoint present, no final `.gro` yet) — **not** a failure.
- `FAIL` = output genuinely missing. If the job also isn't in `squeue`, it died — go to §3/§4.
- Goal state for a batch is `SWITCH (N/N ok)`.

Queue + accounting:
```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID%20,State,ExitCode,Elapsed,Timelimit -X | grep -Ev 'COMPLETED'
# which node each array element ran on (spot a bad node):
sacct -j <JOBID> --format=JobID%16,State,AllocTRES%45,NodeList%12 -X
```

Job logs (this exact path — we wasted time guessing it once):
```
logs/pmx_neq/<stage>/pmx_neq_<stage>_c0.<JOBID>_<ARRAYIDX>.err   # and .out
```
```bash
for idx in <failed indices>; do echo "== $idx =="; tail -n 20 logs/pmx_neq/<stage>/pmx_neq_<stage>_c0.<JOBID>_${idx}.err; done
```

---

## 3. Failure catalogue (symptom → cause → fix)

| Symptom in log / sacct | Cause | Fix |
|---|---|---|
| `future feature annotations is not defined` | ran with login-node Python 3.6 | `source scripts/sherlock/activate_pmx_env.sh` (§0) |
| `Device ID 0 did not correspond to any of the 0 detected device(s)` | landed on a node whose GPU is not visible to CUDA (**bad node**, not our code) | resubmit; node usually self-invalidates. §4 |
| `NVML: ... GPU is lost` then `TIMEOUT` after 12 h + `Unkillable job step` | GPU died mid-run; step stuck in D-state until wall clock | resubmit the element; the node goes `inval` on its own. §4 |
| equil element runs ~12 h (others ~20–30 min) then FAILS/TIMEOUT | CPU fallback / sick GPU on that node | same as above — §4 |
| `OMP_NUM_THREADS` disagrees with `-ntomp`, mdrun aborts at min | inherited OMP env | already fixed (`run_neq_task.py` pins `OMP_NUM_THREADS=ntomp`). If it recurs, check that fix is present |
| extract error `../eq_lambdaN/equil.trr does not exist` | historical path bug | fixed (extract paths use `../../eq_lambda`). Should not recur |
| switch produced `switch.xvg` but audit wants `dgdl.xvg` | GROMACS `-deffnm switch` names it `switch.xvg` | `run_neq_task.py` self-heals (renames). If reprocessing old runs, copy `switch.xvg`→`dgdl.xvg` |
| equil job dies in ~10 s with `Missing equil trajectory ... need equil.trr` (an *extract* error) or `task_id NN not found in <manifest>` | **task-id file collision between concurrent batches** — see §4c | cancel + re-submit with a per-batch `TASK_ID_FILE` |
| `combine_neq`/`qc_neq` reuse stale numbers after a re-run | cached `analysis.json` | pass `--force` to `analyze_neq`/`combine_neq`/`qc_neq` |
| qc_neq reports ~0.00 overlap / huge dissipation | **artifact** — do not negate reverse work; pmx stores W_R in the forward frame | already fixed in `qc_neq.py`. Overlap uses W_f vs W_r directly |

### Diagnosing a bad GPU node
```bash
sacct -j <JOBID> --format=JobID%16,State,NodeList%12 -X   # do all failures share one NodeList? → bad node
sinfo -n <node> -o "%n %t %E %G"                          # STATE inval/drain/down = out of pool already
```
If all failures share one node and that node is `inval`/`drain`/`down`, it is **already unschedulable**
— new jobs cannot land there, so no exclusion is needed.

**Do not assume a bad node self-invalidates — check `sinfo`.** Known repeat offenders that stay
schedulable:

| node | observed | `sinfo` after the failures |
|---|---|---|
| `sh03-12n12` | hangs equil to the wall; **also** fails switch in ~15 min with the device error | stays schedulable |
| `sh03-12n13` | fails switch with the device error (9 s to 4.5 h) — killed 4 K103N elements 2026-08-16 | `mix none gpu:4` — still advertising 4 GPUs |

Exclude both explicitly on any resubmit.

**The "hang to the wall" and "Device ID 0 ... 0 detected device(s)" symptoms are the SAME node
fault.** When a node's GPUs vanish, a job *already running* blocks forever on a dead device and burns
to the wall clock, while a job *starting afterwards* aborts at once because CUDA reports zero devices.
The symptom is decided by timing, not by a different defect — so treat either as "this node's GPUs
are gone" and exclude it. Consequence for a bundled switch task: it can also die **mid-bundle** after
hours, having written valid `dgdl.xvg` for the switches it already finished.

Nodes with `gpu:4` can be *partially* broken — on 2026-08-16 `sh03-12n13` failed one element in 9 s
while two others ran on it for 4.5 h. So a node is not exonerated by other elements succeeding there.
SLURM does requeue elements off such a node on its own (G190E elements 2 and 3 moved and completed
elsewhere without intervention), which is why letting an array drain before topping up is usually
better than cancelling.

Note the elapsed time when diagnosing: a GPU-visibility error is an *mdrun startup* failure and would
be instant, so an element that dies after **hours** lost its GPU **mid-bundle**. Its already-written
`dgdl.xvg` files are valid and complete — the idempotent resubmit skips them and redoes only the
remainder, so the loss is much smaller than "N failed elements x bundle size" suggests.

**Node exclusion only works at SUBMIT time, via an EXPORTED `EXCLUDE_NODES`.**

```bash
export EXCLUDE_NODES=sh03-12n12          # must be exported, see below
STAGE=equil bash scripts/fep_pmx/submit_p0_neq.sh
scontrol show job <id> | grep -io 'ExcNodeList=[^ ]*'   # verify: NOT (null)
```

Three ways this silently fails to apply — all verified 2026-08-15, all leaving
`ExcNodeList=(null)`:

1. **`scontrol update JobId=... ExcNodeList=...` on a pending job is a silent
   no-op.** It **returns exit 0**, so a `&&` chain proceeds as if it worked, and
   the value is simply never set. There is no way to add an exclusion after submit
   — resubmit or accept the risk.
2. **`SBATCH_EXCLUDE`** is ignored (the script builds its own `SBATCH_ARGS`).
3. **Assigning without `export`.** `submit_p0_neq.sh` runs as a child process and
   reads `EXCLUDE_NODES` from the *environment*, so a plain assignment does not
   reach it — while `echo $EXCLUDE_NODES` in your shell still shows the value,
   which makes this look fine. Beware especially of a missing space, as in
   `exportMANIFEST=x EXCLUDE_NODES=y`: bash parses that as two ordinary
   assignments (junk var `exportMANIFEST`, plus a NON-exported `EXCLUDE_NODES`),
   the line exits 0, and the `&&` chain continues with no exclusion applied.

**Verify with `SubmitLine`, NOT `ExcNodeList`.** On this cluster `ExcNodeList` reads `(null)` even
when `--exclude` is demonstrably on the sbatch command line — it does not reflect `--exclude` at all,
so treating `(null)` as failure sends you chasing a non-problem (it did, for two days):

```bash
scontrol show job -d <id> | grep -o 'SubmitLine=.*' | head -1  # THE reliable check
sacct -j <id> --format=JobID%16,State,NodeList%14 -X           # ground truth: no element on a bad node
```

**Use `scontrol`, not `sacct`, for this.** `--exclude` is appended LAST to the sbatch
command line, and `sacct --format=SubmitLine%250` truncates before reaching it (the
tell is a trailing `+` on the printed line) -- so it prints nothing and the exclusion
looks missing when it is actually present. `sacct ... SubmitLine%600` also works if
you prefer sacct.

Confirmed 2026-08-17 on job 39630950: `SubmitLine` carried `--exclude=sh03-12n12,sh03-12n13`,
`ExcNodeList` still printed `(null)`, and none of the allocated elements landed on either node.

And weigh the fix: cancelling a pending
job to add an exclusion forfeits accrued queue priority, whereas a node hang is
recoverable via §4. When the queue is the binding constraint, letting it run is
usually the better trade.

An access scare is almost never real: GPU **access** failures stop you at *scheduling* (job rejected,
or PD with a QOS/Assoc reason). If jobs got `gres/gpu=1` allocated and some **ran**, access is fine —
the rest is node health. Account/QOS check: `sacctmgr -n -P show assoc user=$USER format=Account,QOS`
(currently `rshafer | long,normal`).

---

## 4. The standard recovery: one bad element stalled the batch

This is the recipe for the most common incident — a few equil elements died on a bad GPU node, so
extract/switch are stuck `DependencyNeverSatisfied`. Because stages are idempotent (§1.3), the fix is
to cancel the doomed downstream and **re-run + re-chain** the whole equil→extract→switch. The 33 good
equil units skip instantly; only the failures recompute (on healthy nodes).

```bash
cd /scratch/users/rsatija/nnrti-mechanisms-git
source scripts/sherlock/activate_pmx_env.sh
export MANIFEST=results/analysis/fep_pmx/neq_p1a_manifest.csv   # the batch's manifest
export EXCLUDE_NODES=sh03-12n12                                 # avoid the known-bad node (repeat offender)

# 1) cancel the downstream jobs whose afterok can never be met
scancel <extract_jobid> <switch_jobid>

# 2) re-run equil (idempotent: completed units skip; failures recompute on healthy nodes)
EQUIL=$(STAGE=equil bash scripts/fep_pmx/submit_p0_neq.sh | tail -1);        echo "equil=$EQUIL"

# 3) re-chain extract (CPU) then switch (GPU)
EXTRACT=$(STAGE=extract DEPENDENCY=afterok:$EQUIL   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "extract=$EXTRACT"
SWITCH=$(STAGE=switch   DEPENDENCY=afterok:$EXTRACT bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "switch=$SWITCH"

# 4) verify the chain wired up (dependencies, not the transient squeue REASON):
for j in $EQUIL $EXTRACT $SWITCH; do echo "== $j =="; scontrol show job $j | grep -ioE 'Dependency=[^ ]*'; done
squeue -u $USER
```

Verification notes:
- Right after submit, `squeue` may show a dependent job's REASON as `(None)` for one scheduling cycle
  before it flips to `(Dependency)`. Trust `scontrol show job ... | grep Dependency`, not the transient
  REASON column.
- A dead node may keep an old element in `CG` (completing) for a long time — cosmetic, ignore it; it
  clears when admins reset the node.
- If a *re-run* element lands on **another** bad node and fails, the downstream `afterok` re-blocks —
  just repeat §4. It cannot corrupt completed work.

---

## 4b. Pre-flighting a fresh batch without a GPU

Before spending GPU-hours on a new panel, validate the inputs on CPU. The `em` stage is the built-in
GPU-free check: it runs `grompp` (topology/mdp/structure validation) + a CPU minimization on the
`normal` partition. Submit it alone; if it reaches `EM (N/N ok)` the inputs are sound and the GPU
stages (which reuse the same topology) will run. Only then chain equil→extract→switch.

```bash
export MANIFEST=results/analysis/fep_pmx/neq_<batch>_manifest.csv
STAGE=em bash scripts/fep_pmx/submit_p0_neq.sh                 # CPU only
python -m nnrti.fep.audit_neq_panel --manifest $MANIFEST | grep '==='   # want EM (18/18 ok)
# then the GPU stages, chained on the validated em:
EQUIL=$(STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
EXTRACT=$(STAGE=extract DEPENDENCY=afterok:$EQUIL   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
SWITCH=$(STAGE=switch  DEPENDENCY=afterok:$EXTRACT bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
```

Even cheaper first pass (zero jobs): confirm the manifest shape and that each unit's `system.top` +
`mdp/` exist on disk — a quick `csv`/`os.path` scan of `neq_<batch>_manifest.csv`. `missing inputs: 0`
means prep is complete for every leg/phase/rep.

New neutral P1 batch prep (mirrors how P1a/P1b were built):
```bash
python -m nnrti.fep.prepare_neq --legs wt_to_<A> wt_to_<B> wt_to_<C> \
  --replicates 3 --n-snapshots 100 \
  --panel-manifest results/analysis/fep_pmx/neq_<batch>_manifest.csv
```
Note the flag is `--panel-manifest` for `prepare_neq.py`, but `MANIFEST=` (env var) for
`submit_p0_neq.sh`. Neutral single legs use 100 ps switches (not in `LONG_SWITCH_LEGS`).

---

## 4c. Running TWO batches at once — set `TASK_ID_FILE` or you WILL corrupt one

**This bit us on 2026-08-15 and cost a full G190E equil wave.**

`submit_p0_neq.sh` writes the stage's task ids to
`results/analysis/fep_pmx/neq_<stage>_task_ids.txt` (+ `_chunk000`), a path with
**no batch qualifier** (`submit_p0_neq.sh:92`). Array elements read that file **at
runtime, not at submit time**. So submitting the same stage for a second batch
*overwrites the file out from under the first batch's still-queued elements*, and
they then execute whatever task id sits on their line — from the other batch.

What it looked like: a G190E equil array whose elements 0–6 completed normally
(they ran before K103N was submitted) while 7–11 died in ~10 s with either
`Missing equil trajectory ... need equil.trr` — the *extract* stage's error, in an
equil job — or `task_id 46 not found in ... neq_g190e_manifest.csv`.

**The dangerous part is what did NOT error.** K103N-500 has 9 tasks/unit (54 ids)
versus G190E's 7 (42), so ids ≥ 42 fall out of range and fail loudly — but ids
0–41 resolve *silently against the wrong manifest* and run a valid-but-unintended
task. We were lucky the failures landed out of range. Equil output is per-unit and
idempotent so a wrong-but-valid equil task is harmless, but do not count on that
for other stages.

**Fix — give every batch its own task-id file:**

```bash
export MANIFEST=results/analysis/fep_pmx/neq_<batch>_manifest.csv
EQUIL=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_equil_<batch>_task_ids.txt \
        STAGE=equil bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
EXTRACT=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_extract_<batch>_task_ids.txt \
        STAGE=extract DEPENDENCY=afterok:$EQUIL bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
SWITCH=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_switch_<batch>_task_ids.txt \
        STAGE=switch DEPENDENCY=afterok:$EXTRACT bash scripts/fep_pmx/submit_p0_neq.sh | tail -1)
```

Verify from the submit output: it must print
`file=neq_<stage>_<batch>_task_ids_chunk000.txt`. If it prints the bare
`neq_<stage>_task_ids_chunk000.txt`, the env var did not take — **stop**, or you
re-corrupt the other batch.

Corollary: `MANIFEST` is an exported shell var reused by both batches, so a stale
export silently audits/submits the wrong one. Pass `--manifest` explicitly to
`audit_neq_panel.py`. The switch count distinguishes them (100 ps = 2 switch tasks
per unit, 500 ps = 4).

**Recovery** is §4's: `scancel` the affected array plus its dependants (`scancel`
touches no filesystem, so it is safe even during a scratch outage), then re-submit
with per-batch `TASK_ID_FILE`. Stages are idempotent, so completed units skip.

## 4d. `$SCRATCH` degraded / hardware issue

The login banner reports cluster status; `- $SCRATCH: Hardware issue` means the
filesystem holding the repo and all job output is sick. Symptoms: commands that
stat many files hang and are **unkillable** (Ctrl-C cannot interrupt a process in
uninterruptible I/O wait), and jobs may requeue with `NODE_FAIL`.

- **Do not submit** into a degraded filesystem. Partial writes produce truncated
  outputs and completion markers that lie — worse than a clean failure, because
  the idempotency checks then skip a unit that never really finished.
- `scancel` is safe (SLURM control plane only).
- Check https://status.sherlock.stanford.edu before resuming.
- **Afterwards, verify anything that was mid-write.** A truncated `dgdl.xvg`
  parses into a wrong work value rather than an error:
  ```bash
  find results/analysis/fep_pmx/legs/<leg> -name 'dgdl.xvg' -exec wc -l {} + | sort -n | head
  ```
  All switches of a given length should have the same line count; short files are
  suspect. Same for `equil.gro` presence vs `equil.cpt`-only.

## 5. When a batch reaches `SWITCH (N/N ok)` — analysis

**Run the first analysis pass ON SHERLOCK, not the Mac.** The rsync (`scripts/rsync_fep_pmx.sh`)
deliberately **excludes `dgdl.xvg`** (light provenance only). `analyze_neq`/`combine_neq` need either the
`dgdl.xvg` *or* a pre-computed `analysis.json`. So on a freshly-finished batch the Mac has neither and
`combine_neq` silently **skips** every leg (`skip <G>: No fwd dgdl.xvg files`). The fix: run
`combine_neq` on Sherlock first — it auto-runs `analyze_neq` (which reads the `dgdl.xvg` and writes
`analysis.json` + `integ_{fwd,rev}.dat`), *then* rsync those light outputs to the Mac.

```bash
# 1) ON SHERLOCK (dgdl live here): generate analysis.json + the ΔΔG table.
python -m nnrti.fep.combine_neq --targets <G1> <G2> ... --replicates 3
python -m nnrti.fep.qc_neq --legs wt_to_<G1> wt_to_<G2> ... --replicates 3
# 2) THEN from the Mac: pull the now-existing light outputs and inspect/re-run locally (no --force).
SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull
```

Interface quirks that waste time (verified):
- **`combine_neq` takes `--targets` (genotypes, `V106I`); `qc_neq` takes `--legs` (leg ids,
  `wt_to_V106I`).** Different flags for the same set — don't mix them up.
- **`qc_neq` with no `--legs` defaults to `P0_LEGS`** (prints V106A/Y188L) — if you see the P0 legs when
  you expected new ones, you forgot `--legs`.
- **First analysis run is ~1 min *per unit*, not seconds** (it parses 200 switch `dgdl.xvg` per unit —
  tens of thousands of dH/dλ points each). A 3-leg batch = 18 units ≈ 10–30 min. It writes one
  `analysis.json` per unit as it goes, so watch progress with
  `find results/analysis/fep_pmx/legs/wt_to_<G> -name analysis.json | wc -l` (out of `6*n_legs`).
  Subsequent runs read the cached `analysis.json` and are fast — **unless** you pass `--force`, which
  re-parses (and, on a Mac with numpy≥2, breaks pmx's estimators; keep `numpy<2` locally, or just don't
  `--force` on the Mac).

Then free GPUs are available for the next batch (and the apo 100 ns extension — see `STATUS.md`).

---

## 6. Building out the panel — preparing NEW legs (hybrids → systems → neq)

The full manuscript panel is `scripts/fep_jorgensen/mutations.py::MANUSCRIPT_PLANS` (19 genotypes, 19
unique legs). Experimental folds: `results/analysis/dor_susceptibility_bar_chart/tables/dor_susceptibility_values.csv`.
Both `prepare_p0_hybrids.sh` and `build_p0_systems.sh` take a `LEGS="..."` override — that's the
mechanism for any leg, P0/P1/P2. Full prep for a set of legs (Sherlock login, pmx env):

```bash
LEGS="wt_to_G190S V106A_to_V106A_F227L ..."                                  # leg_ids from MANUSCRIPT_PLANS
LEGS="$LEGS" bash scripts/fep_pmx/prepare_p0_hybrids.sh                       # pmx hybrids, reps 1-3
REPLICATES="1 2 3" LEGS="$LEGS" bash scripts/fep_pmx/build_p0_systems.sh      # gmx solvate+ionize
python -m nnrti.fep.prepare_neq --legs $LEGS --replicates 3 --n-snapshots 100 \
  --panel-manifest results/analysis/fep_pmx/neq_<batch>_manifest.csv
```
All three are idempotent (SKIP existing). Then run via §4b (em pre-flight) + the em→…→switch chain.

**Gotchas that will otherwise cost you hours (all hit while preparing the P2 compounds):**

1. **Compound legs seed from an *intermediate* genotype.** `V106A_to_V106A_F227L` needs the MD
   `_start.pdb` for BOTH V106A (source) and V106A+F227L (endpoint), holo+apo, all reps. Check first:
   resolve `leg.input_complex_pdb / input_apo_pdb / endpoint_*_pdb` and confirm they exist.

2. **`results/md_runs/**/_start.pdb` live on the Mac, not Sherlock** (never git-committed; single legs
   didn't need them). Push once, Mac→Sherlock (small, but ~1.9 GB for all — compresses ~5×):
   ```bash
   rsync -avz --prune-empty-dirs --include='*/' --include='*_start.pdb' --exclude='*' \
     results/md_runs/ rsatija@login.sherlock.stanford.edu:/scratch/users/rsatija/nnrti-mechanisms-git/results/md_runs/
   ```

3. **Legs without a `fep_jorgensen` backend map need openmm in the pmx venv.** `prepare_hybrid` resolves
   the mutation site from `results/analysis/fep_jorgensen/legs/<leg_id>/prepare_backend.json` if present
   (the openmm-free path — most singles + A98G+F227C have one); otherwise it falls through to
   `resolve_mutation_site`, which imports openmm to read the PDBs and **verify** the unique old→new
   change (do NOT hand-write backend maps to skip this — a wrong residue id silently mutates the wrong
   site). One-time install into the pmx venv:
   ```bash
   pip install --no-deps openmm           # openmm 8.1.1 wheel
   pip install 'numpy<2'                  # REQUIRED: openmm 8.1.1's xtc module won't import under numpy 2
   python3 -c "import openmm.app; print('ok')"
   ```
   Safe to do mid-panel: the SLURM tasks (`run_neq_task.py`) import no numpy/pmx-python — they shell out
   to `gmx` — so a venv numpy change can't affect running jobs. (numpy<2 is also what `pmx analyse`
   wants; the venv had drifted to 2.0.2.)

4. **Proline mutations — FIXED (commit c21415b).** Mutating a proline (e.g. P225H, Pro→His) used to
   crash pmx `mutate` (`_set_conformation` → `IndexError`) while copying the A-state coords: proline's
   `HG2/HG3, HD2/HD3` must be `HG1/HG2, HD1/HD2` for pmx's hybrid, but `normalize_openmm_for_pmx` only
   fixed `HB`/`HA`. It now also renames proline HG/HD, **PRO-scoped** (His ring `HD1/HD2`, Asn/Gln
   amide, Arg/Lys methylenes keep their names) and idempotent. Verified end-to-end — mutate builds the
   P2H hybrid and `pdb2gmx` accepts the rename for all prolines (not just the mutated one). If a future
   non-proline residue hits the same `old_res[name]` IndexError, it's the same class of OpenMM→pmx
   methylene-naming gap (`2/3`→`1/2`); extend the rename residue-scoped the same way.

## 7. Lesson: validate convergence, not just execution (charge legs)

A protocol that **runs** and is formally **correct** can still be statistically **non-convergent** — and
that gap cost us a full charge-leg run. The co-alchemical ion (for K103N/G190E) kept the box neutral
(correct) and passed every execution gate (grompp accepts it, equil integrates, switches finish), so it
*looked* validated. But it doesn't converge: decoupling a **whole ion** (charge + LJ, a tens-of-kcal
transformation) in a 100–500 ps switch dissipates **~20–26 kcal/mol**, versus **~1–3** for a neutral
mutation. That ~10× dissipation drives the forward/reverse work distributions ~20 kcal apart with
near-zero overlap, so BAR cannot converge (SEM ~1.4; BAR–Jarzynski disagreement up to ~3.7).

The tell that it was the *ion*, not the mutation: dissipation was ~constant (~20 kcal) across K103N and
G190E even though their ΔGs differ wildly (~9 vs ~36) — the one thing they share is the Cl⁻ decoupling.

**Rules going forward:**
- The co-alchemical ion **relocates** the charge perturbation (mutation → bulk ion); it does not remove
  it. Annihilating a full ion is intrinsically dissipative in fast NEQ switching. **Do not use it here.**
- For charge-changing legs use the **analytical net-charge correction** (Rocklin/Hünenberger): run the
  leg raw (non-neutral box), then add a closed-form finite-size term — **zero** added simulation
  perturbation, so switch dissipation stays neutral-like.
- Before committing GPU to any *new* alchemical protocol, do a 30-second estimate: how large is the
  alchemical change being added, and can the switch length drive it near-reversibly? Then read the
  **overlap / dissipation QC** (`qc_neq`, `integ_{fwd,rev}.dat`) — "it finished" is not "it converged."
