# G190E — closing the SEM (2026-08-28)

**Goal:** bring `wt_to_G190E` from **ΔΔG = +2.00 ± 1.77** to SEM < 1 kcal/mol, so
G190E can carry a main-text point estimate in Table 2 / Supplementary Figure 2
instead of being flagged as unresolved.

Run on Sherlock at `$PROJECT_ROOT`. Agents cannot
reach Sherlock — paste output back. Read [`OPERATIONS.md`](OPERATIONS.md) before
acting; this runbook only covers what is specific to G190E.

Prerequisite every session:

```bash
cd $PROJECT_ROOT && git pull
source ops/slurm/cluster/activate_pmx_env.sh
```

---

## Where the error actually is

Current leg (verified from `neq_prepare.json`, all six units):
**500 ps switches, 5 ns endpoint equilibration, 100 snapshots, 3 replicates.**

Per-replicate BAR (kcal/mol), from the cached `analysis.json`:

| rep | holo ΔG | apo ΔG | **ΔΔG** | holo hyst | apo hyst | holo−apo hyst |
|---|---:|---:|---:|---:|---:|---:|
| 1 | −60.32 | −59.61 | **−0.71** | 19.33 | 10.20 | **+9.13** |
| 2 | −49.65 | −51.04 | **+1.39** | 20.37 | 21.64 | −1.27 |
| 3 | −47.69 | −53.01 | **+5.32** | 16.62 | 21.80 | −5.18 |

σ_DDG = **3.06** → SEM = 3.06/√3 = 1.77.

Three things follow, and they rule out two of the four available levers:

1. **Not switch length.** 500 ps was already the fix that worked for `wt_to_K103N`
   (σ_DDG 3.80 → 0.40) and it did not generalize here. `sd(holo−apo hysteresis)`
   is **7.39** for G190E versus 0.62 for K103N-500 — the triage rule in
   [`STATUS.md`](STATUS.md) puts anything > 5 in the "phases decoupled, longer
   switches will not be enough" bin. G190E still dissipates 16–22 kcal/mol at
   500 ps because Gly→Glu grows a charged carboxylate out of dummy atoms.
2. **Not more snapshots.** Per-rep `bar_err_boot` is 0.41; the within-replicate
   error is already an order of magnitude below the between-replicate spread.
   Tightening each replicate just makes three replicates confidently disagree.
3. **The endpoint ensembles differ between replicates.** Absolute per-phase ΔG
   moves by 12 kcal/mol across reps (holo −60.3/−49.7/−47.7) and holo and apo
   fail to track each other in rep 1 and rep 3. That is under-equilibrated /
   under-sampled endpoints, which leaves exactly two levers:

   * **A — longer endpoint equilibration** (`NEQ_EQUIL_NS`, currently 5 ns).
     Attacks σ_DDG at the root. Untested for this leg. Uncertain payoff.
   * **B — more replicates.** SEM = σ_DDG/√n is arithmetic, so this always
     works: n = 6 → 1.25, n = 10 → 0.97. Currently **blocked** (see below).

**Run A and B concurrently.** They are independent, they use different queues at
different times, and A changes what B costs: if A cuts σ_DDG to ~1.7, then n = 6
lands at 0.69 and the campaign ends early.

---

## Why B is blocked, and how to unblock it

`seed_extra_replicates.py` builds a rep_04 seed from a run's `*_md_final.pdb`. A
leg needs four of them — (source, endpoint) × (holo, apo). For `wt_to_G190E`:

| structure | path | status |
|---|---|---|
| holo source | `results/md_runs/wt/rep_01/wt_rep01_md_final.pdb` | ok |
| holo endpoint | `results/md_runs/G190E/rep_01/G190E_rep01_md_final.pdb` | ok |
| apo source | `results/md_runs/apo/wt/rep_01/wt_rep01_md_final.pdb` | ok |
| **apo endpoint** | `results/md_runs/apo/g190e/rep_01/` | **only `assets/` — apo G190E MD never ran** |

Verified on the Mac (`seed_extra_replicates.py --legs wt_to_G190E` reports
`apo endpoint BLOCKED ... no *_md_final.pdb and no usable analysis trajectory`).
Check Sherlock too before spending GPU — the Mac copy is rsync-filtered:

```bash
ls -la results/md_runs/apo/g190e/rep_0*/
```

If Sherlock also has only `assets/`, run the apo MD. The prepared system XMLs
already exist, so this is a submission, not a preparation:

```bash
MUTATION_ALLOWLIST=g190e MD_PRODUCTION_NS=100.0 \
  bash ops/slurm/cluster/submit_apo_md_batched.sh 3 3
```

Monitor with `python3 ops/slurm/cluster/report_md_progress.py --target-ns 100.0 --show-incomplete`.
Three 100 ns apo runs, ~38 GPU-h each under 12 h walls with checkpoint resume.

---

## Lever A — 20 ns endpoint equilibration, reps 1–3

Switch length, snapshots, charge handling and seeds are all unchanged; only
`NEQ_EQUIL_NS` moves. `prepare_neq.py --force` re-renders mdps and manifests and
**never deletes** em/equil/extract/switch outputs — so the existing 5 ns
`equil.gro` files would be treated as complete and skipped. **Archive them first**
or the run silently reuses 5 ns equilibration.

```bash
# A1 — archive the 5 ns equilibration so it is not silently reused.
ARCH=results/analysis/fep_pmx/_archive/wt_to_G190E_equil5ns_2026-08-28
mkdir -p $ARCH
for ph in holo apo; do for r in 01 02 03; do for l in 0 1; do
  d=results/analysis/fep_pmx/legs/wt_to_G190E/$ph/rep_$r/neq/eq_lambda$l
  [ -d "$d" ] && mkdir -p $ARCH/$ph/rep_$r && mv "$d" $ARCH/$ph/rep_$r/
done; done; done
find $ARCH -name equil.gro | wc -l          # want 12

# A2 — 20 ns equil, 500 ps switches, switches re-bundled at 10/task.
# All three env vars are required on THIS command. They are NOT needed on the
# later stages: submit_p0_neq.sh selects rows straight out of the manifest and
# run_neq_task.py reads the rendered .mdp files off disk, so the 20 ns and the
# bundling are already baked into what prepare_neq wrote.
NEQ_EQUIL_NS=20 NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E NEQ_SNAPSHOTS_PER_TASK_LONG=10 \
python -m nnrti.fep.prepare_neq --legs wt_to_G190E \
  --replicates 3 --n-snapshots 100 --force \
  --panel-manifest results/analysis/fep_pmx/neq_g190e_equil20_manifest.csv

# A3 — GATE. Do not skip; a missing env var has silently produced a wrong config twice.
python3 -c "import json; p=json.load(open('results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/neq_prepare.json')); print('switch_ps',p['switch_ps'],'equil_ns',p['equil_ns'],'n_tasks',p['n_tasks'],'|','OK' if p['switch_ps']==500.0 and p['equil_ns']==20.0 and p['n_tasks']==25 else 'WRONG')"
grep -H nsteps results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/mdp/npt_eq_lambda0.mdp \
               results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/mdp/nonequil_fwd.mdp
```

Want `switch_ps 500.0 equil_ns 20.0 n_tasks 25 | OK`, then
`npt_eq_lambda0.mdp: nsteps = 10000000` (20 ns at 2 fs -- `2500000` means
`NEQ_EQUIL_NS` did not take) and `nonequil_fwd.mdp: nsteps = 250000` (500 ps).

`prepare_neq` prints the same single `Wrote panel NEQ manifest` line whether or
not `--force` actually caused it to redo the work -- the printed output does not
distinguish the two cases, only the A3 gate does. Do not read that line as
confirmation.

**`n_tasks` is 25 with the bundling above**: 1 em + 2 equil + 2 extract + 20
switch (100 snapshots / 10 per task x 2 lambda). At the *default* 50/task it is
9 -- do not accept that; see "Switch bundling" below. (An earlier draft of this
runbook said 19, copied from a stale `neq_prepare.json`. It is wrong.)

### Switch bundling -- set it, do not take the default

At `SWITCH_SNAPSHOTS_PER_TASK_LONG = 50` each switch element runs 50 x 500 ps
sequentially. `config.py` already records what that costs: *"on 2026-08-17
G190E's last 96 switches sat in 2 elements at ~21 h; re-bundled at 6/task they
spread over 16 GPUs (~3 h)"*. 21 h overruns a 12 h wall.

| snapshots/task | switch elements (6 units) | est. h/element |
| ---: | ---: | ---: |
| 50 (default) | 24 | ~21 |
| 20 | 60 | ~8.4 |
| **10 (used here)** | **120** | **~4.2** |
| 6 | 204 | ~2.5 |

Re-bundling is safe and idempotent: `prepare_neq --force` rewrites only mdps and
manifests, and `run_neq_task` skips per-*switch* on an existing `dgdl.xvg`, not
per-task, so completed work is never redone.

```bash
# A4 — CPU pre-flight, then chain the GPU stages with PER-BATCH task-id files.
export MANIFEST=results/analysis/fep_pmx/neq_g190e_equil20_manifest.csv
TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_em_g190e_e20_task_ids.txt \
STAGE=em bash ops/slurm/fep/submit_p0_neq.sh
python -m nnrti.fep.audit_neq_panel --manifest $MANIFEST | grep '==='   # want EM (6/6 ok)

export EXCLUDE_NODES=sh03-12n12,sh02-16n06,sh03-13n01
EQUIL=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_equil_g190e_e20_task_ids.txt \
        STAGE=equil   bash ops/slurm/fep/submit_p0_neq.sh | tail -1); echo "equil=$EQUIL"
EXTRACT=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_extract_g190e_e20_task_ids.txt \
        STAGE=extract DEPENDENCY=afterok:$EQUIL   bash ops/slurm/fep/submit_p0_neq.sh | tail -1); echo "extract=$EXTRACT"
SWITCH=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_switch_g190e_e20_task_ids.txt \
        STAGE=switch  DEPENDENCY=afterok:$EXTRACT bash ops/slurm/fep/submit_p0_neq.sh | tail -1); echo "switch=$SWITCH"
for j in $EQUIL $EXTRACT $SWITCH; do echo "== $j =="; scontrol show job $j | grep -ioE 'Dependency=[^ ]*'; done
```

Every submit must print `file=neq_<stage>_g190e_e20_task_ids_chunk000.txt`. A bare
`neq_<stage>_task_ids_chunk000.txt` means `TASK_ID_FILE` did not take — **stop and
fix**, or a concurrently-submitted batch will corrupt this one mid-flight
(this destroyed a full G190E equil wave on 2026-08-15; OPERATIONS §4c).

Set `SHERLOCK_TIME` to at least `12:00:00` for equil: 20 ns is ~7h20m/task
against ~1h50m at 5 ns, and 12 elements run concurrently. Switch is unchanged at
~4 h/element over 120 elements at 10 snapshots/task.

---

## Lever B — replicates 4–6

Only after apo G190E MD has written `*_md_final.pdb`. Seeds are per-replicate:
rep_04 from MD rep_01, rep_05 from rep_02, rep_06 from rep_03, so all three are
independent conformations rather than three re-runs of one.

```bash
# B1 — seed. Dry run first; it validates atom count, PBC-split and ligand contact.
for s in 1 2 3; do
  python -m nnrti.fep.seed_extra_replicates --legs wt_to_G190E \
    --source-rep $s --dest-rep $((s+3))
done
# then repeat with --apply once all three report ok
```

```bash
# B2 — hybrids + solvated systems for the new reps only.
LEGS="wt_to_G190E" REPLICATES="4-6" bash ops/slurm/fep/prepare_p0_hybrids.sh
source ops/slurm/cluster/load_gromacs_module.sh && source ops/slurm/cluster/activate_pmx_env.sh
echo "$GMXLIB"          # must end in .../mutff
LEGS="wt_to_G190E" REPLICATES="4-6" bash ops/slurm/fep/build_p0_systems.sh
```

`REPLICATES` in the shell scripts is a **list/range**; `--replicates` in
`prepare_neq.py` is a **count**. Same word, opposite meaning.

```bash
# B3 — NEQ prep for reps 4-6. Match whatever equilibration Lever A settled on.
NEQ_EQUIL_NS=20 NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E \
python -m nnrti.fep.prepare_neq --legs wt_to_G190E \
  --rep-start 4 --replicates 6 --n-snapshots 100 --force \
  --panel-manifest results/analysis/fep_pmx/neq_g190e_reps456_manifest.csv
```

Then repeat A3 (gate, on `rep_04`) and A4 (submit) against
`neq_g190e_reps456_manifest.csv` with its own `TASK_ID_FILE` names
(`..._g190e_r456_...`). 3 leg-reps at 500 ps = 12 equil + 24 switch elements.

---

## Analysis — on Sherlock, and never with `--force`

The rsync excludes `dgdl.xvg`, so the Mac has nothing to analyse until
`analysis.json` exists. `--force` re-analyses every genotype, and the other 17
lost their raw `dgdl.xvg` to the 2026-08-13 `git clean -fd` — it would fail on
missing inputs.

```bash
python -m nnrti.fep.combine_neq --targets G190E --replicates 6
python -m nnrti.fep.qc_neq --legs wt_to_G190E --replicates 6
```

⚠️ `combine_neq --targets X` **rewrites `panel_ddg.csv` with only those targets.**
Rebuild the full panel afterwards with the complete genotype list:

```bash
python -m nnrti.fep.combine_neq --targets \
  F227C G190A G190E G190S V106A V106I V106M Y181C Y188L Y318F \
  A98G+F227C V106A+F227L V106A+L234I V106A+P225H V106I+F227C \
  K103N K103N+M230L K103N+P225H L100I+K103N --replicates 3
```

Pull to the Mac with `SHERLOCK_USER=rsatija bash ops/sync/rsync_fep_pmx.sh pull`.

---

## Read-outs, in order of what they decide

1. **σ_DDG after Lever A** (spread of the three per-rep ΔΔG, not the SEM). 3.06 →
   under 1.7 means n = 3 alone gets SEM < 1 and Lever B is unnecessary.
2. **`sd(holo−apo hysteresis)`** after A. 7.39 → under ~2 is the mechanism
   working: residual dissipation bias became correlated between phases and
   cancels in the double difference. This is what changed for K103N-500 (3.00 →
   0.62) while overlap stayed at ~0, so do not judge by overlap.
3. **Point-estimate stability.** The ΔΔG must stay within error of the two
   independent existing estimates — +2.00 ± 1.77 (500 ps / 5 ns) and +2.50 ± 1.41
   (archived co-alchemical run). A large shift means something other than noise
   changed.
4. **BAR vs Jarzynski vs CGI** within ~1 kcal/mol, and per-rep `bar_err_boot` in
   single digits. The archived co-alchemical arm had apo `bar_err` 35.49 — that
   is the non-convergence signature.

## Manuscript hooks

- `results/analysis/fep_pmx/panel_ddg.csv` → Table 2 `∆∆Gbind` column, rebuilt by
  `src/nnrti/cli/build_table_2.py`.
- `results/analysis/fep_pmx/panel_discussion_tiers.csv` → the main-text /
  show / omit tiering. `OMIT_MAIN_TEXT` in `combine_neq.py` still hardcodes
  `{"K103N", "G190E"}` from when both had SEM > 2; K103N is now 0.23. **Revisit
  that set** once G190E lands.
- Supplementary Figure 2 work distributions:
  `src/nnrti/cli/plot_fep_work_distributions.py`.
- Results text currently reads "in case of G190E because of large replicate
  variance ∆∆Gbind = 2.00 ± 1.77 kcal/mol"; Discussion reads "in some cases
  (e.g. G190E), this was not enough to mitigate the errors". Both need rewriting
  if the SEM drops.

---

## Lever C — repair holo rep 2's lambda=1 ensemble (2026-08-30)

After the 20 ns rebuild G190E sits at **+0.99 +/- 1.63**, sigma_DDG 2.82 (was
3.06). The equilibration lever worked on its own terms -- sd(holo-apo
hysteresis) fell 7.39 -> 4.52 -- but one unit carries the whole variance.

### The diagnosis, from the work distributions

Per-replicate ddG: **-0.57 / +4.24 / -0.71**. Replicates 1 and 3 agree to
0.14 kcal/mol; replicate 2 is 4.9 away. Mean work by direction (kcal/mol):

| phase | direction | rep 1 | rep 2 | rep 3 | spread |
| --- | --- | ---: | ---: | ---: | ---: |
| holo | forward (from WT, lambda=0) | -43.87 | -39.64 | -35.73 | 8.14 |
| holo | **reverse (from mutant, lambda=1)** | **-64.07** | **-55.43** | **-62.56** | 8.63 |
| apo | forward | -37.69 | -40.79 | -39.76 | 3.10 |
| apo | reverse | -64.09 | -61.13 | -64.26 | 3.12 |

In the forward direction rep 2 sits *in the middle* of a smooth 8 kcal spread --
ordinary scatter in the WT ensemble. In the reverse direction reps 1 and 3 agree
within 1.5 kcal and **rep 2 is 7.1 kcal away**. The fault is therefore the
**G190E mutant (lambda=1) endpoint ensemble of holo rep 2**, not the WT side,
not the switching, and not dissipation.

Corroborating, and all pointing the same way:

- rep 2's hysteresis is the **smallest** of the six units (15.79 vs 20-27) and
  its work range the **narrowest** (16.8 vs 24.5 / 31.0) -- it is internally
  well converged, just converged somewhere else;
- it is the only unit where BAR and Jarzynski disagree (-1.20; every other unit
  within 0.05), which is one of this pipeline's stated confidence criteria;
- no bimodality in any unit (largest sorted-work gap 1.3-3.9 kcal in ranges of
  17-31), so this is not a second basin being sampled within a run;
- the skew in every unit is driven by ~5 switches and vanishes when they are
  dropped, so rep 2 is not uniquely tailed.

This is the lambda=1 noisiness `STATUS.md` has flagged since the P0 pilot,
caught in a specific unit.

### The repair

**Re-run only holo rep 2's lambda=1 equilibration and everything downstream of
it** -- `eq_lambda1`, `snapshots/lambda1`, and the 10 reverse switch bundles.
The lambda=0 side of that unit is sound and its forward switches are reusable,
so this is about half the elements of a full unit re-run.

`npt_warmup.mdp` sets `gen_vel = yes` with no `gen_seed`, so GROMACS defaults to
`-1` and the re-run draws **independent velocities**. It is a genuine second
sample of the lambda=1 ensemble, not a repeat.

**Both outcomes are publishable, which is the point of doing it:**

- lands near -63 (like reps 1 and 3) -> the original was a fluke, sigma_DDG
  collapses, and G190E reports near **-0.6 +/- 0.1** -- sign flipped, meaning
  binding does *not* explain its 18-fold resistance;
- lands near -55 again -> rep 2's ensemble is real, sigma_DDG is genuinely ~2.8,
  the leg needs n ~ 8, and the paper reports G190E as unresolved with a reason.

What is *not* publishable is the current state: one unit failing an independent
convergence check, unexplained.

### Failure-minimisation, given an empty queue

Every one of these cost real time earlier in this campaign.

1. **Do not chain stages with `afterok`.** One transient
   `cudaErrorLaunchFailure` on 2026-08-28 made the dependency unsatisfiable and
   stalled everything behind it for a day. With a free queue there is no
   scheduling reason to pre-chain -- submit equil, verify, submit extract,
   verify, submit switch.
2. **Archive the stale `analysis/` directory BEFORE running**, not after.
   `combine_neq` reuses a cached `analysis.json` and will silently return the
   old numbers otherwise -- it did exactly that on 2026-08-30.
3. **Exclude the nodes that have faulted**: `sh03-12n12,sh02-16n06,sh03-13n01,sh03-12n01`.
4. **Per-batch `TASK_ID_FILE` on every stage.** Array elements resolve it at
   runtime; an unqualified name let a concurrent batch corrupt a live array on
   2026-08-15.
5. **Generous wall time.** Equil at 20 ns is ~7h20m/task; request 24 h.
6. **Verify `equil.gro` line counts** after equil -- a truncated one is treated
   as complete and skipped.

### Apo G190E MD -- give it a wall it can finish in

The apo runs have now failed to complete **twice**, both by wall-clock
truncation (12 h wall against a ~38 GPU-h job), leaving no output JSON either
time. They are the only thing blocking Lever B (extra replicates), and extra
replicates need *both* phases: pairing is doing real work here (paired sigma_DDG
2.82 against 5.69 unpaired), so holo-only replicates would not help.

**48 h is the hard ceiling, and the QOS cannot lift it.** The `long` QOS allows
7 days, but the `gpu` *partition* has `MaxTime=2-00:00:00`, and the partition
limit wins -- `SHERLOCK_QOS=long SHERLOCK_TIME=72:00:00` is rejected with
"Requested time limit is invalid (missing or exceeds some limit)". Verify with
`scontrol show partition gpu | grep -oE "MaxTime=[^ ]*"`.

So run on the default QOS at the ceiling; `long` would only cost you its
4-concurrent-job cap for no benefit:

```bash
SHERLOCK_TIME=48:00:00 MUTATION_ALLOWLIST=g190e \
  MD_PRODUCTION_NS=100.0 bash ops/slurm/cluster/submit_apo_md_batched.sh 3 3
```

48 h is 4x the wall that truncated these twice, against a ~38 GPU-h job, and the
runs resume from checkpoint so the ~20 GPU-h already spent counts toward the
100 ns.
