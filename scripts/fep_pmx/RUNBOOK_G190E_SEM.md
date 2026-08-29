# G190E — closing the SEM (2026-08-28)

**Goal:** bring `wt_to_G190E` from **ΔΔG = +2.00 ± 1.77** to SEM < 1 kcal/mol, so
G190E can carry a main-text point estimate in Table 2 / Supplementary Figure 2
instead of being flagged as unresolved.

Run on Sherlock at `/scratch/users/rsatija/nnrti-mechanisms-git`. Agents cannot
reach Sherlock — paste output back. Read [`OPERATIONS.md`](OPERATIONS.md) before
acting; this runbook only covers what is specific to G190E.

Prerequisite every session:

```bash
cd /scratch/users/rsatija/nnrti-mechanisms-git && git pull
source scripts/sherlock/activate_pmx_env.sh
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
  bash scripts/sherlock/submit_apo_md_batched.sh 3 3
```

Monitor with `python3 scripts/sherlock/report_md_progress.py --target-ns 100.0 --show-incomplete`.
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
python3 scripts/fep_pmx/prepare_neq.py --legs wt_to_G190E \
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
STAGE=em bash scripts/fep_pmx/submit_p0_neq.sh
python3 scripts/fep_pmx/audit_neq_panel.py --manifest $MANIFEST | grep '==='   # want EM (6/6 ok)

export EXCLUDE_NODES=sh03-12n12,sh02-16n06,sh03-13n01
EQUIL=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_equil_g190e_e20_task_ids.txt \
        STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "equil=$EQUIL"
EXTRACT=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_extract_g190e_e20_task_ids.txt \
        STAGE=extract DEPENDENCY=afterok:$EQUIL   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "extract=$EXTRACT"
SWITCH=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_switch_g190e_e20_task_ids.txt \
        STAGE=switch  DEPENDENCY=afterok:$EXTRACT bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "switch=$SWITCH"
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
  python3 scripts/fep_pmx/seed_extra_replicates.py --legs wt_to_G190E \
    --source-rep $s --dest-rep $((s+3))
done
# then repeat with --apply once all three report ok
```

```bash
# B2 — hybrids + solvated systems for the new reps only.
LEGS="wt_to_G190E" REPLICATES="4-6" bash scripts/fep_pmx/prepare_p0_hybrids.sh
source scripts/sherlock/load_gromacs_module.sh && source scripts/sherlock/activate_pmx_env.sh
echo "$GMXLIB"          # must end in .../mutff
LEGS="wt_to_G190E" REPLICATES="4-6" bash scripts/fep_pmx/build_p0_systems.sh
```

`REPLICATES` in the shell scripts is a **list/range**; `--replicates` in
`prepare_neq.py` is a **count**. Same word, opposite meaning.

```bash
# B3 — NEQ prep for reps 4-6. Match whatever equilibration Lever A settled on.
NEQ_EQUIL_NS=20 NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E \
python3 scripts/fep_pmx/prepare_neq.py --legs wt_to_G190E \
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
python3 scripts/fep_pmx/combine_neq.py --targets G190E --replicates 6
python3 scripts/fep_pmx/qc_neq.py --legs wt_to_G190E --replicates 6
```

⚠️ `combine_neq --targets X` **rewrites `panel_ddg.csv` with only those targets.**
Rebuild the full panel afterwards with the complete genotype list:

```bash
python3 scripts/fep_pmx/combine_neq.py --targets \
  F227C G190A G190E G190S V106A V106I V106M Y181C Y188L Y318F \
  A98G+F227C V106A+F227L V106A+L234I V106A+P225H V106I+F227C \
  K103N K103N+M230L K103N+P225H L100I+K103N --replicates 3
```

Pull to the Mac with `SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull`.

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
  `src/analysis/cli/build_table_2.py`.
- `results/analysis/fep_pmx/panel_discussion_tiers.csv` → the main-text /
  show / omit tiering. `OMIT_MAIN_TEXT` in `combine_neq.py` still hardcodes
  `{"K103N", "G190E"}` from when both had SEM > 2; K103N is now 0.23. **Revisit
  that set** once G190E lands.
- Supplementary Figure 2 work distributions:
  `src/analysis/cli/plot_fep_work_distributions.py`.
- Results text currently reads "in case of G190E because of large replicate
  variance ∆∆Gbind = 2.00 ± 1.77 kcal/mol"; Discussion reads "in some cases
  (e.g. G190E), this was not enough to mitigate the errors". Both need rewriting
  if the SEM drops.
