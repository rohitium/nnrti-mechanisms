# G190E run — 500 ps switches (charge-leg protocol)

Completes the 18-genotype panel. **G190E, not G190S.**

Run on a Sherlock login node at `$PROJECT_ROOT`.
Agents cannot reach Sherlock — paste output back.

> **2026-08-15: switched from 100 ps to 500 ps.** The original run was pinned to
> 100 ps to match K103N, but K103N is itself now re-running at 500 ps (its 100 ps
> arm is archived under `_archive/wt_to_K103N_100ps_2026-08-15/`), so 500 ps *is*
> the matched charge-leg configuration. Rationale for the charge family generally:
> `wt_to_K103N` dissipates 14–19 kcal/mol with forward/reverse work distributions
> 4.2–5.1σ apart — the regime where BAR is biased, not merely noisy — unlike the
> neutral legs (1–4 kcal/mol) that established switch-length invariance.
> Equil and extract are **switch-length independent**, so any equil already done at
> the 100 ps setting is reused unchanged; only the nonequil mdps and the manifest
> differ.

## Protocol

| setting | value | source |
|---|---|---|
| charge handling | raw non-neutral box + Rocklin/Hunenberger analytical correction | `USE_COALCHEMICAL_ION = False` |
| switch length | **500 ps** | `NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E` (env, not committed) |
| snapshots / endpoint / rep | 100 | `NEQ_SNAPSHOTS_DEFAULT` |
| endpoint equilibration | 5 ns | `NEQ_EQUIL_NS` |
| replicates | 3 | |
| warmup | 500 ps C-rescale -> P-R | `NEQ_WARMUP_PS` |
| delta_q | -1 (Gly0 -> Glu-) | `CHARGE_LEG_DELTA_Q` |

Everything is idempotent (SKIP existing) — a re-run costs nothing.

**The committed default for `wt_to_G190E` is still 100 ps** (it is not in
`_BASE_LONG_SWITCH_LEGS`). The 500 ps run is driven by the env var above, which
must be set on `prepare_neq.py`. If K103N-500 is adopted panel-wide, move
`wt_to_G190E` and `wt_to_K103N` into `_BASE_LONG_SWITCH_LEGS` and drop the env var.

## Is a 100 ps arm needed too? Normally no.

The deliverable is a confident ΔΔG, and 500 ps gives that directly. A matched
100 ps arm (switch stage only, reusing the same snapshots — cheap, since
equilibration is the expensive part) is only worth running if **K103N-500 comes
back ambiguous**: hysteresis falls but σ_DDG does not, or the two move
inconsistently. Then the matched-snapshot pair is the only way to separate "switch
length didn't help" from "re-equilibration hurt."

Attribution does not otherwise need it: **dissipation is a property of the driving
speed, not of endpoint equilibration**, so a hysteresis drop between the archived
K103N-100 and the new K103N-500 is attributable to switch length even though those
two arms differ in equilibration as well.

## Step 0 — sync

```bash
cd $PROJECT_ROOT
git pull
```

## Step 1 — activate pmx env (the gotcha that bites every session)

Login-node `python3` is 3.6 and dies on `from __future__ import annotations`.

```bash
source ops/slurm/cluster/activate_pmx_env.sh
```

## Step 2 — confirm the source PDBs are on Sherlock

`wt_to_G190E` seeds from the **WT** MD start structures (holo + apo, reps 1-3).
All 12 exist on the Mac; they were never git-committed, so verify before prep.

**Case matters.** The holo WT directory is lowercase `results/md_runs/wt/`
(`safe_label("WT")` returns `"wt"`), while mutant directories keep their
uppercase genotype label (`G190E/`) and apo directories are all lowercase
(`apo/wt/`, `apo/g190e/`). macOS is case-insensitive so a wrong-case path passes
on the Mac and fails only on Sherlock — always verify with the exact case below.

```bash
ls results/md_runs/wt/rep_0*/assets/wt_md_rep0*_start.pdb \
   results/md_runs/apo/wt/rep_0*/assets/wt_apo_md_rep0*_start.pdb \
   results/md_runs/G190E/rep_0*/assets/G190E_md_rep0*_start.pdb \
   results/md_runs/apo/g190e/rep_0*/assets/g190e_apo_md_rep0*_start.pdb
```

Want 12 files. If any are missing, push from the Mac:

```bash
rsync -avz --prune-empty-dirs --include='*/' --include='*_start.pdb' --exclude='*' \
  results/md_runs/ <user>@<cluster>:$PROJECT_ROOT/results/md_runs/
```

## Step 3 — hybrids

```bash
LEGS="wt_to_G190E" REPLICATES="1-3" bash ops/slurm/fep/prepare_p0_hybrids.sh
```

`REPLICATES` here is a **list/range**, not a count — `REPLICATES=3` means rep 3
only. The script guards this and warns, but use `"1-3"`.

## Step 4 — solvated systems

`build_p0_systems.sh` needs both GROMACS and pmx; load in this order so `GMXLIB`
ends up pointing at pmx's `mutff`.

```bash
source ops/slurm/cluster/load_gromacs_module.sh && source ops/slurm/cluster/activate_pmx_env.sh
echo "$GMXLIB"    # must end in .../mutff
LEGS="wt_to_G190E" REPLICATES="1-3" bash ops/slurm/fep/build_p0_systems.sh
```

G190E is a **growth** mutation (Gly -> Glu): the B-state sidechain is a dummy in
A, so it is never touched by the global A-state `em`. It is relaxed by the per-λ
minimization inside the equil stage (`run_neq_task.py::_run_equil` step 1,
PLAN §4.3) — this is expected, not a missing step.

## Step 5 — NEQ prep

```bash
NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E \
python -m nnrti.fep.prepare_neq --legs wt_to_G190E \
  --replicates 3 --n-snapshots 100 --force \
  --panel-manifest results/analysis/fep_pmx/neq_g190e500_manifest.csv
```

Note the asymmetry: `--replicates 3` here IS a count (reps 1..3), unlike the
shell `REPLICATES`. Same word, opposite meaning.

`NEQ_EXTRA_LONG_SWITCH_LEGS` is what makes this a 500 ps run — the committed
default for this leg is 100 ps. **It must be set on this command**; forgetting it
silently renders 100 ps mdps, which step 5b catches.

**`--force` is required here, and it is safe.** `prepare_neq` skips any unit that
already has a `neq_manifest.csv`, so without it the run returns in a second having
written only the panel manifest and leaving the previous switch length in place.
`--force` re-renders mdps + manifests and **never** deletes em/equil/extract/switch
outputs (comment at `prepare_neq.py:120`) — which is exactly why any equil already
completed at the 100 ps setting is preserved and reused: **equil and extract do not
depend on switch length.**

## Step 5b — the gate (do not skip; this has caught a wrong config twice)

The JSON:

```bash
python3 -c "import json; p=json.load(open('results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/neq_prepare.json')); print('switch_ps',p['switch_ps'],'| n_tasks',p['n_tasks'],'|','OK' if p['switch_ps']==500.0 and p['n_tasks']==9 else 'WRONG')"
```

Want `switch_ps 500.0 | n_tasks 9 | OK`.

And the rendered mdp, since that — not the JSON — is what GROMACS reads:

```bash
grep -E 'nsteps|delta.lambda' results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/mdp/nonequil_fwd.mdp
```

Want `nsteps = 250000` (500 ps at 2 fs) and `delta-lambda = 4e-06`. `50000` means
the env var did not take.

## Step 5c — if resuming after an interrupted run, check the surviving equil

A truncated `equil.gro` is worse than a missing one: `run_neq_task.py` treats its
presence as completion and **skips** the unit. After any crash, wall-out, or
filesystem incident:

```bash
for ph in holo apo; do for r in 01 02 03; do for l in 0 1; do
  d=results/analysis/fep_pmx/legs/wt_to_G190E/$ph/rep_$r/neq/eq_lambda$l
  printf "%-4s rep%s l%s: " $ph $r $l
  if [ -f "$d/equil.gro" ]; then echo "$(wc -l < $d/equil.gro) lines"; else echo "MISSING"; fi
done; done; done
```

Every present `equil.gro` must have the same line count (same atom count). Delete
the directory of any outlier so the re-run regenerates it.

## Step 6 — CPU pre-flight (no GPU burned)

`em` runs grompp + a CPU minimization: it validates topology/mdp/structure. If
this passes, the GPU stages reuse the same topology and will run.

```bash
export MANIFEST=results/analysis/fep_pmx/neq_g190e500_manifest.csv
STAGE=em bash ops/slurm/fep/submit_p0_neq.sh
python -m nnrti.fep.audit_neq_panel --manifest $MANIFEST | grep '==='   # want EM (6/6 ok)
```

## Step 7 — chain equil -> extract -> switch

One leg at 500 ps = 12 equil + 24 switch GPU array elements (500 ps bundles 50
snapshots per task, so 4 switch tasks per unit instead of 2).

**Set a per-batch `TASK_ID_FILE` on every stage.** Array elements resolve the
task-id file at *runtime*, so a concurrently-submitted batch that writes the
default unqualified path will corrupt this one mid-flight — this destroyed a full
G190E equil wave on 2026-08-15. See [`OPERATIONS.md`](OPERATIONS.md) §4c.

```bash
export MANIFEST=results/analysis/fep_pmx/neq_g190e500_manifest.csv
export EXCLUDE_NODES=sh03-12n12,sh02-16n06,sh03-13n01

EQUIL=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_equil_g190e500_task_ids.txt \
        STAGE=equil bash ops/slurm/fep/submit_p0_neq.sh | tail -1);              echo "equil=$EQUIL"
EXTRACT=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_extract_g190e500_task_ids.txt \
        STAGE=extract DEPENDENCY=afterok:$EQUIL   bash ops/slurm/fep/submit_p0_neq.sh | tail -1); echo "extract=$EXTRACT"
SWITCH=$(TASK_ID_FILE=$PWD/results/analysis/fep_pmx/neq_switch_g190e500_task_ids.txt \
        STAGE=switch DEPENDENCY=afterok:$EXTRACT bash ops/slurm/fep/submit_p0_neq.sh | tail -1); echo "switch=$SWITCH"

for j in $EQUIL $EXTRACT $SWITCH; do echo "== $j =="; scontrol show job $j | grep -ioE 'Dependency=[^ ]*'; done
```

Each submit must print `file=neq_<stage>_g190e500_task_ids_chunk000.txt`. A bare
`neq_<stage>_task_ids_chunk000.txt` means the env var did not take — **stop**.

Trust `scontrol ... Dependency`, not squeue's transient `REASON` column.

Rough timing: equil ~1h50m/task at 5 ns; switch ~3-9 h/task at 500 ps.

## Step 8 — monitor

```bash
python -m nnrti.fep.audit_neq_panel --manifest $MANIFEST | grep '==='
python -m nnrti.fep.audit_neq_panel --manifest $MANIFEST | grep -E 'FAIL|RUNNING'
squeue -u $USER
```

`RUNNING` = mid-flight, not a failure. Goal state: `SWITCH (12/12 ok)`.

**If an equil element dies**, `afterok` is unsatisfiable and everything
downstream hangs as `DependencyNeverSatisfied`. Fix = cancel downstream, then
resubmit + re-chain exactly as in step 7 (completed units skip in seconds). Full
recipe: `OPERATIONS.md` §4.

## Step 9 — analysis, ON SHERLOCK

The rsync deliberately excludes `dgdl.xvg`, so the Mac has nothing to analyze
until this runs and writes `analysis.json`.

```bash
python -m nnrti.fep.combine_neq --targets G190E --replicates 3
python -m nnrti.fep.qc_neq --legs wt_to_G190E --replicates 3
```

- `combine_neq` takes `--targets` (genotype); `qc_neq` takes `--legs` (leg id).
- **Do NOT pass `--force`.** The other 17 genotypes' raw `dgdl.xvg` was destroyed
  by the 2026-08-13 `git clean -fd`; `--force` re-analyzes everything and would
  fail on missing inputs. Without it, cached `analysis.json` is reused and only
  G190E (which now has none — the stale co-alchemical copies were archived) gets
  analyzed.
- First pass is ~1 min per unit (12 units, so ~10-20 min). Watch with
  `find results/analysis/fep_pmx/legs/wt_to_G190E -name analysis.json | wc -l` (of 6).

Then from the Mac:

```bash
SHERLOCK_USER=rsatija bash ops/sync/rsync_fep_pmx.sh pull
```

## Step 10 — sanity checks on the result

Before it goes on the panel:

1. **Hysteresis is K103N-like, not ion-like.** Measure
   `<W_f> - <W_r>` from `integ_{fwd,rev}.dat` (see the step 11 snippet). Under the
   raw-box protocol K103N sits at **14-19 kcal/mol** holo — elevated versus
   neutral legs (1-8) but far below the co-alchemical ion's ~20-26 with its
   catastrophic per-rep errors. Do NOT expect the "1-3 kcal, neutral-like" figure
   quoted in OPERATIONS §7; that describes neutral mutations, not a delta_q = -1
   leg run raw.
2. **Per-rep BAR errors are sane.** The archived co-alchemical run had an apo
   `bar_err` of 35.49 — that is the non-convergence signature. Single digits or
   below is fine.
3. **BAR vs Jarzynski agree** within ~1 kcal/mol.
4. **Hysteresis should be well below the 100 ps figures.** At 500 ps expect
   roughly 3-4 kcal/mol and separation near 2 sigma, versus 14-19 and 4.2-5.1
   sigma for K103N at 100 ps. That drop is the mechanism this protocol is buying;
   if hysteresis is unchanged, the switch length did not take (re-check step 5b).
5. **SEM is the deliverable.** Under 1 kcal/mol is the bar. The K103N family sat
   at 2.19-2.32 at 100 ps — that is what we are trying to beat.

## Step 11 (OPTIONAL) — the matched 100 ps arm

**Do not run this by default.** 500 ps is the deliverable; this is a contingency.

Run it only if **K103N-500 comes back ambiguous** — hysteresis drops but sigma_DDG
does not, or the two move inconsistently. Then a matched-snapshot pair is the only
way to separate "switch length didn't help" from "re-equilibration hurt", because
K103N's 100 ps and 500 ps arms differ in equilibration too (its original snapshots
were destroyed by the 2026-08-13 `git clean`).

It is otherwise unnecessary: **dissipation is a property of the driving speed, not
of endpoint equilibration**, so a hysteresis drop between archived K103N-100 and
new K103N-500 is attributable to switch length despite that confound.

Cheap when needed — switch stage only, reusing snapshots that equilibration
already paid for. Invert the recipe below: clear the switch outputs, re-prepare
**without** `NEQ_EXTRA_LONG_SWITCH_LEGS` (so it renders 100 ps), and re-submit
`STAGE=switch` alone with its own `TASK_ID_FILE`.

### Background: why the charge legs run long at all

The panel's switch-length-invariance result (V106A: +1.69 +- 0.70 at 100 ps vs
+1.76 +- 0.51 at 500 ps; F227C likewise) was measured on **neutral** legs that
dissipate 1-4 kcal/mol — they had nothing to gain. The charge leg is a different
regime. From the surviving work files:

| leg | switch | hysteresis | separation (hyst / sigma_work) |
|---|---|---|---|
| wt_to_K103N | 100 ps | **14.2-19.4** | **4.2-5.1 sigma** |
| wt_to_Y181C | 100 ps | 7.9-10.8 | 2.3-4.8 |
| wt_to_G190S | 100 ps | 4.5-6.5 | 3.9-5.2 |
| wt_to_V106A | 500 ps | 1.1-4.2 | 1.2-3.2 |
| wt_to_G190A | 100 ps | -0.5-1.5 | small |

At 4-5 sigma the forward/reverse work distributions barely overlap, which is
where BAR becomes both noisy and biased. And K103N's per-replicate ΔG tracks its
per-replicate dissipation almost 1:1:

```
wt_to_K103N holo:  hyst 14.19 -> bar_dg -75.71
                   hyst 18.54 -> bar_dg -77.99
                   hyst 19.36 -> bar_dg -79.79
```

If BAR had removed the dissipation bias those would be independent. So part of
the charge family's sigma_DDG = 3.80 is *varying residual bias*, which longer
switches attack directly: dissipation falls ~linearly with switch time while
work spread falls as its square root, so 100 -> 500 ps should improve separation
by ~sqrt(5) (~4.5 sigma -> ~2 sigma).

**Stakes:** all four K103N genotypes inherit their ~2.2 SEM from the single
`wt_to_K103N` leg (`K103N_to_K103N_M230L`'s own sigma_DDG is 0.59). If 500 ps
works here, rebuilding that one leg fixes four of the panel's five worst points.

### Recipe for the optional 100 ps arm (switch stage only)

Reuses G190E's own equil/extract snapshots, so the two arms differ **only** in
switch length — the controlled comparison K103N cannot provide.

```bash
source ops/slurm/cluster/activate_pmx_env.sh

# preserve the 500 ps result before overwriting the switch outputs
cp -r results/analysis/fep_pmx/legs/wt_to_G190E \
      results/analysis/fep_pmx/_archive/wt_to_G190E_500ps_$(date +%F)

# clear ONLY the switch outputs; em/equil/extract are untouched and reused
rm -rf results/analysis/fep_pmx/legs/wt_to_G190E/*/rep_*/neq/switches/*
find results/analysis/fep_pmx/legs/wt_to_G190E -name analysis.json -delete

# re-render mdps at 100 ps: note NO NEQ_EXTRA_LONG_SWITCH_LEGS here
python -m nnrti.fep.prepare_neq --legs wt_to_G190E \
  --replicates 3 --n-snapshots 100 --force \
  --panel-manifest results/analysis/fep_pmx/neq_g190e100_manifest.csv

# verify before submitting
python3 -c "import json;print(json.load(open('results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/neq_prepare.json'))['switch_ps'])"   # want 100.0

MANIFEST=results/analysis/fep_pmx/neq_g190e500_manifest.csv EXCLUDE_NODES=sh03-12n12 \
  STAGE=switch bash ops/slurm/fep/submit_p0_neq.sh
```

Remember the per-batch `TASK_ID_FILE` on that switch submission too
(OPERATIONS §4c).

100 ps switches bundle 100 snapshots per task, so this arm is 12 array elements at
~1-2 h each — far cheaper than the 500 ps run it is being compared against.

If switches emit `switch.xvg` instead of `dgdl.xvg`, `run_neq_task.py` self-heals;
for older units recover with the `find ... -name switch.gro` loop in `README.md`.

### Reading the result

```bash
python -m nnrti.fep.combine_neq --targets G190E --replicates 3
python -m nnrti.fep.qc_neq --legs wt_to_G190E --replicates 3

# hysteresis + separation, current run vs whatever is archived
python3 - <<'PY'
import statistics as st
from pathlib import Path
KJ = 4.184
def w(root, ph, r, which):
    p = Path(f'{root}/{ph}/rep_{r:02d}/neq/analysis/integ_{which}.dat')
    return [float(l.split()[-1]) / KJ for l in open(p)] if p.exists() else None
import glob
old = sorted(glob.glob('results/analysis/fep_pmx/_archive/wt_to_G190E_*ps_*'))
for tag, root in [('current', 'results/analysis/fep_pmx/legs/wt_to_G190E')] + \
                 ([('archived', old[-1])] if old else []):
    print(f'== {tag}')
    for ph in ('holo', 'apo'):
        for r in (1, 2, 3):
            f, rv = w(root, ph, r, 'fwd'), w(root, ph, r, 'rev')
            if not f: continue
            h = st.mean(f) - st.mean(rv)
            pooled = ((st.stdev(f)**2 + st.stdev(rv)**2) / 2) ** 0.5
            print(f'  {ph:5s} rep{r} hyst={h:7.2f}  sep={h/pooled:5.2f} sigma')
PY
```

**Decision rule.** The quantity that matters is **sigma_DDG** (the spread of the
three per-replicate ΔΔG values), not hysteresis itself — hysteresis is the
mechanism, sigma_DDG is the payoff:

- **sigma_DDG drops materially** (say 3.8 -> under ~2) -> adopt 500 ps for charge
  legs and rebuild `wt_to_K103N` the same way; four genotypes improve at once.
- **sigma_DDG roughly unchanged** while hysteresis clearly falls -> the residual
  bias was not the driver; the spread is conformational, and the lever is
  endpoint equilibration (`NEQ_EQUIL_NS`) plus replicates, per `STATUS.md`.
- **Neither moves** -> stop spending on this leg; go straight to more replicates,
  which always works (SEM = sigma_DDG/sqrt(n)) and needs n ~ 15 for the K103N family.

Caveat on reading it: with n = 3, sigma_DDG carries roughly +-40% relative
uncertainty, so only a large change is interpretable. A drop from 3.8 to 3.0 is
noise; 3.8 to 1.5 is not.

## History — what was changed and why (2026-08-15)

1. **Stale co-alchemical G190E artifacts archived** (ΔΔG 2.50 ± 1.41, apo
   `bar_err` 35.49) to
   `results/analysis/fep_pmx/_archive/wt_to_G190E_coalchemical_500ps_2026-08-15/`.
   Left in place, `combine_neq` without `--force` would have reused them and put a
   non-convergent G190E on the panel silently.
2. **`config.py`: `wt_to_G190E` removed from `_BASE_LONG_SWITCH_LEGS`** (500 → 100 ps)
   to match K103N. That original 500 ps entry was speculative, from the
   co-alchemical era.
3. **First run cancelled** — `scancel 39256102 39256103 39256104` — after the
   task-id file collision (OPERATIONS §4c) killed equil elements 7–11. Elements
   0–6 completed and their equil output is reusable.
4. **Switched back to 500 ps**, this time on evidence rather than assumption:
   K103N is being re-run at 500 ps, so 500 ps is now the matched charge-leg
   configuration. Driven by `NEQ_EXTRA_LONG_SWITCH_LEGS`, leaving the committed
   default at 100 ps until K103N-500 justifies a code change.

The net lesson worth carrying: the 500 ps entry was right by accident originally,
wrong for the reason it was written, and is right again for a different, measured
reason (14–19 kcal/mol dissipation at 4.2–5.1σ separation). Steps 5b and 6 exist
because two of the three configuration errors in this leg's history were caught by
a gate rather than by inspection.
