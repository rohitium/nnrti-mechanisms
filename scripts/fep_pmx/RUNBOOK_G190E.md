# G190E run — K103N-matched protocol

Completes the 18-genotype panel. **G190E, not G190S.**

Run on a Sherlock login node at `/scratch/users/rsatija/nnrti-mechanisms-git`.
Agents cannot reach Sherlock — paste output back.

## Protocol (identical to K103N and the three K103N-compound legs)

| setting | value | source |
|---|---|---|
| charge handling | raw non-neutral box + Rocklin/Hunenberger analytical correction | `USE_COALCHEMICAL_ION = False` |
| switch length | **100 ps** | `config.py` (G190E removed from `_BASE_LONG_SWITCH_LEGS`) |
| snapshots / endpoint / rep | 100 | `NEQ_SNAPSHOTS_DEFAULT` |
| endpoint equilibration | 5 ns | `NEQ_EQUIL_NS` |
| replicates | 3 | |
| warmup | 500 ps C-rescale -> P-R | `NEQ_WARMUP_PS` |
| delta_q | -1 (Gly0 -> Glu-) | `CHARGE_LEG_DELTA_Q` |

Everything is idempotent (SKIP existing) — a re-run costs nothing.

## Step 0 — sync the code change

The G190E switch-length change lives in `config.py` and must reach Sherlock
before prep, or G190E rebuilds at 500 ps.

```bash
cd /scratch/users/rsatija/nnrti-mechanisms-git
git pull
grep -n "_BASE_LONG_SWITCH_LEGS = " scripts/fep_pmx/config.py   # want: {"wt_to_V106A", "wt_to_Y188L"}
```

## Step 1 — activate pmx env (the gotcha that bites every session)

Login-node `python3` is 3.6 and dies on `from __future__ import annotations`.

```bash
source scripts/sherlock/activate_pmx_env.sh
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
  results/md_runs/ rsatija@login.sherlock.stanford.edu:/scratch/users/rsatija/nnrti-mechanisms-git/results/md_runs/
```

## Step 3 — hybrids

```bash
LEGS="wt_to_G190E" REPLICATES="1-3" bash scripts/fep_pmx/prepare_p0_hybrids.sh
```

`REPLICATES` here is a **list/range**, not a count — `REPLICATES=3` means rep 3
only. The script guards this and warns, but use `"1-3"`.

## Step 4 — solvated systems

`build_p0_systems.sh` needs both GROMACS and pmx; load in this order so `GMXLIB`
ends up pointing at pmx's `mutff`.

```bash
source scripts/sherlock/load_gromacs_module.sh && source scripts/sherlock/activate_pmx_env.sh
echo "$GMXLIB"    # must end in .../mutff
LEGS="wt_to_G190E" REPLICATES="1-3" bash scripts/fep_pmx/build_p0_systems.sh
```

G190E is a **growth** mutation (Gly -> Glu): the B-state sidechain is a dummy in
A, so it is never touched by the global A-state `em`. It is relaxed by the per-λ
minimization inside the equil stage (`run_neq_task.py::_run_equil` step 1,
PLAN §4.3) — this is expected, not a missing step.

## Step 5 — NEQ prep

```bash
python3 scripts/fep_pmx/prepare_neq.py --legs wt_to_G190E \
  --replicates 3 --n-snapshots 100 \
  --panel-manifest results/analysis/fep_pmx/neq_g190e_manifest.csv
```

Note the asymmetry: `--replicates 3` here IS a count (reps 1..3), unlike the
shell `REPLICATES`. Same word, opposite meaning.

**Verify the protocol actually matched K103N before spending GPU time:**

```bash
python3 -c "
import json
p = json.load(open('results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/neq_prepare.json'))
k = json.load(open('results/analysis/fep_pmx/legs/wt_to_K103N/holo/rep_01/neq/neq_prepare.json'))
for f in ('switch_ps','equil_ns','n_snapshots','n_tasks'):
    print(f'{f:14s} G190E={p[f]:<8} K103N={k[f]:<8} {\"OK\" if p[f]==k[f] else \"MISMATCH\"}')
"
```

Want `switch_ps=100.0`, `equil_ns=5.0`, `n_snapshots=100`, `n_tasks=7` on both.
`switch_ps=500.0` means step 0 did not land.

## Step 6 — CPU pre-flight (no GPU burned)

`em` runs grompp + a CPU minimization: it validates topology/mdp/structure. If
this passes, the GPU stages reuse the same topology and will run.

```bash
export MANIFEST=results/analysis/fep_pmx/neq_g190e_manifest.csv
STAGE=em bash scripts/fep_pmx/submit_p0_neq.sh
python3 scripts/fep_pmx/audit_neq_panel.py --manifest $MANIFEST | grep '==='   # want EM (6/6 ok)
```

## Step 7 — chain equil -> extract -> switch

One leg = 12 equil + 12 switch GPU array elements, far under the
`MaxSubmitPU=100` cap, so no batching needed.

```bash
export MANIFEST=results/analysis/fep_pmx/neq_g190e_manifest.csv
export EXCLUDE_NODES=sh03-12n12          # repeat offender: hangs equil to the 12 h wall

EQUIL=$(STAGE=equil   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1);              echo "equil=$EQUIL"
EXTRACT=$(STAGE=extract DEPENDENCY=afterok:$EQUIL   bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "extract=$EXTRACT"
SWITCH=$(STAGE=switch   DEPENDENCY=afterok:$EXTRACT bash scripts/fep_pmx/submit_p0_neq.sh | tail -1); echo "switch=$SWITCH"

for j in $EQUIL $EXTRACT $SWITCH; do echo "== $j =="; scontrol show job $j | grep -ioE 'Dependency=[^ ]*'; done
```

Trust `scontrol ... Dependency`, not squeue's transient `REASON` column.

Rough timing: equil ~1h50m/task at 5 ns; switch ~1-2 h/task at 100 ps.

## Step 8 — monitor

```bash
python3 scripts/fep_pmx/audit_neq_panel.py --manifest $MANIFEST | grep '==='
python3 scripts/fep_pmx/audit_neq_panel.py --manifest $MANIFEST | grep -E 'FAIL|RUNNING'
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
python3 scripts/fep_pmx/combine_neq.py --targets G190E --replicates 3
python3 scripts/fep_pmx/qc_neq.py --legs wt_to_G190E --replicates 3
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
SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull
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
4. **Expect SEM > 1 at this configuration.** The K103N family sits at SEM
   2.19-2.32 under exactly this protocol. That is the honest baseline, and the
   input to step 11 — not a reason to discard the run.

## Step 11 — controlled switch-length test (the SEM experiment)

**Only after step 10.** Step 9 gives the panel-consistent, K103N-matched number;
this asks whether 500 ps switches shrink the error bar on a charge leg.

### Why this is worth GPU time (and why the earlier "no" was wrong)

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

### The design

Re-run **only the switch stage**, reusing G190E's own equil/extract snapshots.
Identical endpoint ensembles, only switch length differs — so this isolates
switch length from equilibration, which the V106A test never did.

```bash
source scripts/sherlock/activate_pmx_env.sh

# preserve the 100 ps result before overwriting the switch outputs
cp -r results/analysis/fep_pmx/legs/wt_to_G190E \
      results/analysis/fep_pmx/_archive/wt_to_G190E_100ps_$(date +%F)

# clear ONLY the switch outputs; em/equil/extract are untouched and reused
rm -rf results/analysis/fep_pmx/legs/wt_to_G190E/*/rep_*/neq/switches/*
find results/analysis/fep_pmx/legs/wt_to_G190E -name analysis.json -delete

# re-render mdps at 500 ps + rebuild the manifest (never deletes em/equil/extract)
NEQ_EXTRA_LONG_SWITCH_LEGS=wt_to_G190E \
python3 scripts/fep_pmx/prepare_neq.py --legs wt_to_G190E \
  --replicates 3 --n-snapshots 100 --force \
  --panel-manifest results/analysis/fep_pmx/neq_g190e500_manifest.csv

# verify before submitting
python3 -c "import json;print(json.load(open('results/analysis/fep_pmx/legs/wt_to_G190E/holo/rep_01/neq/neq_prepare.json'))['switch_ps'])"   # want 500.0

MANIFEST=results/analysis/fep_pmx/neq_g190e500_manifest.csv EXCLUDE_NODES=sh03-12n12 \
  STAGE=switch bash scripts/fep_pmx/submit_p0_neq.sh
```

`NEQ_EXTRA_LONG_SWITCH_LEGS` is the env hook that adds a leg to the long-switch
set without a code edit — deliberately used here so the committed default stays
100 ps / K103N-matched. 500 ps switches bundle 50 snapshots per task
(`SWITCH_SNAPSHOTS_PER_TASK_LONG`), so expect ~2x the array elements at ~3-9 h each.

If switches emit `switch.xvg` instead of `dgdl.xvg`, `run_neq_task.py` self-heals;
for older units recover with the `find ... -name switch.gro` loop in `README.md`.

### Reading the result

```bash
python3 scripts/fep_pmx/combine_neq.py --targets G190E --replicates 3
python3 scripts/fep_pmx/qc_neq.py --legs wt_to_G190E --replicates 3

# hysteresis + separation, 500 ps vs the archived 100 ps
python3 - <<'PY'
import statistics as st
from pathlib import Path
KJ = 4.184
def w(root, ph, r, which):
    p = Path(f'{root}/{ph}/rep_{r:02d}/neq/analysis/integ_{which}.dat')
    return [float(l.split()[-1]) / KJ for l in open(p)] if p.exists() else None
import glob
old = sorted(glob.glob('results/analysis/fep_pmx/_archive/wt_to_G190E_100ps_*'))
for tag, root in [('500ps', 'results/analysis/fep_pmx/legs/wt_to_G190E')] + \
                 ([('100ps', old[-1])] if old else []):
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

## What was changed to make this run correct

- `config.py`: `wt_to_G190E` removed from `_BASE_LONG_SWITCH_LEGS` (500 -> 100 ps),
  so it matches K103N. The 500 ps entry predated the co-alchemical ion being
  abandoned.
- Stale co-alchemical G190E artifacts (ΔΔG 2.50 ± 1.41, apo `bar_err` 35.49)
  moved to
  `results/analysis/fep_pmx/_archive/wt_to_G190E_coalchemical_500ps_2026-08-15/`.
  Left in place, `combine_neq` without `--force` would have reused them and put a
  non-convergent G190E on the panel silently.
