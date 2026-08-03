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
python3 scripts/fep_pmx/audit_neq_panel.py --manifest <manifest.csv> | grep '==='
# which units are incomplete and why:
python3 scripts/fep_pmx/audit_neq_panel.py --manifest <manifest.csv> | grep -E 'FAIL|RUNNING'
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
| `combine_neq`/`qc_neq` reuse stale numbers after a re-run | cached `analysis.json` | pass `--force` to `analyze_neq`/`combine_neq`/`qc_neq` |
| qc_neq reports ~0.00 overlap / huge dissipation | **artifact** — do not negate reverse work; pmx stores W_R in the forward frame | already fixed in `qc_neq.py`. Overlap uses W_f vs W_r directly |

### Diagnosing a bad GPU node
```bash
sacct -j <JOBID> --format=JobID%16,State,NodeList%12 -X   # do all failures share one NodeList? → bad node
sinfo -n <node> -o "%n %t %E %G"                          # STATE inval/drain/down = out of pool already
```
If all failures share one node and that node is `inval`/`drain`/`down`, it is **already unschedulable**
— new jobs cannot land there, so no exclusion is needed. (We confirmed `SBATCH_EXCLUDE` and
`scontrol update ExcNodeList=...` are both silently ignored by this submit path, so **don't rely on
them** — rely on the node self-invalidating, which it does after a lost-GPU/unkillable event.)

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

## 5. When a batch reaches `SWITCH (N/N ok)`

```bash
# pull light provenance (no trajectories) to the Mac for analysis:
SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull      # from the Mac
# ΔΔG_bind per genotype + Spearman vs experiment (auto-runs analyze_neq for missing units):
python3 scripts/fep_pmx/combine_neq.py --targets <G1> <G2> ... --replicates 3
python3 scripts/fep_pmx/qc_neq.py --replicates 3
```
Then free GPUs are available for the next batch (and the apo 100 ns extension — see `STATUS.md`).
