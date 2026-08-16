# pmx NEQ FEP — current state

**Living snapshot. Update this when state changes.** A fresh agent should read this first, then
[`OPERATIONS.md`](OPERATIONS.md) for how to act. Last meaningful update: **2026-08-15**.

Sherlock repo: `/scratch/users/rsatija/nnrti-mechanisms-git` (account `rshafer`, QOS `long,normal`).
The human runs everything on Sherlock and pastes output — agents cannot reach it. Activate the pmx env
first every session: `source scripts/sherlock/activate_pmx_env.sh` (login-node python3 is 3.6 and dies
on `from __future__ import annotations`).

---

## What "success" means here (read this before interpreting any number)

The deliverable is a **reproducible pipeline that produces a high-confidence best estimate of
ΔΔG_bind** — an equilibrium binding observable. **Confidence is defined internally to the simulation,
never by agreement with experiment:** replicate reproducibility (SEM across independent reps), estimator
agreement (BAR/CGI/Jarzynski), switch-length invariance, and endpoint-seeding invariance. Overlap
matters only as it bounds the *precision* of the estimate.

Whether ΔΔG_bind tracks experimental fold-change is a **downstream scientific question we test, not a
pass/fail criterion for the pipeline.** In-vitro resistance can arise from mechanisms other than
equilibrium binding (catalysis, processivity, fitness, conformational effects), so a genotype where a
confident, reproducible ΔΔG does *not* match its fold is a *finding about binding-vs-phenotype*, not a
pipeline defect. Do not call such cases "overpredictions" or "false positives."

---

## Panel status (2026-08-15): n = 17 of 18

`panel_ddg_vs_experiment.png` — **no fit line; the weak correlation IS the finding**: Spearman
ρ = 0.375, Pearson R² = 0.088, p = 0.25 (n = 17 with fold). Tiers in `panel_discussion_tiers.csv` +
`FEP_SECTION_NOTES.md` (main_text SEM ≤ 0.6 / show / omit). Per-genotype protocol walkthroughs
`01`–`05` in `results/analysis/fep_pmx/protocol/<genotype>/` (V106A worked example `protocol_v106a/`),
regen `plot_protocol_figures.py`. Replot scatter/tiers only: `combine_neq.py --replot-only`.

Genotypes in panel: F227C, G190A, G190S, V106A, V106I, V106M, Y181C, Y188L, Y318F, A98G+F227C,
V106A+F227L, V106A+L234I, V106A+P225H, V106I+F227C, K103N, K103N+M230L, K103N+P225H, L100I+K103N.
**Only `G190E` remains** to complete the 18-genotype panel.

**Charge legs now RUN** (supersedes the old "excluded, need co-alchemical ion"). Co-alchemical ion
ABANDONED (didn't converge). Charge legs run in a RAW non-neutral box + Rocklin/Hünenberger analytical
net-charge correction (applied in `combine_neq`). K103N + the 3 K103N-compounds are in the panel;
G190E is the last charge leg to run.

---

## DATA LOSS + REBUILD (READ THIS — 2026-08-13)

A `git clean -fd` deleted the **entire raw `legs/` tree** (~15 GB: equil trajectories,
`switches/*/dgdl.xvg`, `gromacs_build`) on Sherlock scratch. **Not recoverable.** What survived: the
light `analysis.json` + `integ_*.dat` (committed in git) for every already-analyzed leg — so all 17
panel results stand and are reproducible from cached BAR. See memory `data-safety-and-sync`.

**To re-run any leg from scratch** (Sherlock; source structures are the MD
`results/md_runs/<geno>/rep_0N/assets/*_start.pdb`, which survive on Mac + Sherlock):
```bash
source scripts/sherlock/activate_pmx_env.sh
LEGS="<leg_id>" REPLICATES="1-3" bash scripts/fep_pmx/prepare_p0_hybrids.sh
source scripts/sherlock/load_gromacs_module.sh && source scripts/sherlock/activate_pmx_env.sh  # GMXLIB must be .../mutff
LEGS="<leg_id>" REPLICATES="1-3" bash scripts/fep_pmx/build_p0_systems.sh
python3 scripts/fep_pmx/prepare_neq.py --legs <leg_id ...> --replicates 3 --n-snapshots 100 --force
bash scripts/fep_pmx/submit_p0_neq_pipeline.sh   # em→equil→extract→switch, afterok, ≤3-leg batches
```
- `REPLICATES` is a LIST/range (`"1-3"`), NOT a count (`REPLICATES=3` = rep 3 only — guarded now).
- Interactive-GPU smoke first (`salloc_neq_gpu.sh` needs `SHERLOCK_CPUS_PER_TASK=4` — dev QOS caps CPU);
  the smoke looks up its task in `neq_panel_manifest.csv`, so run `prepare_neq --legs ...` before it.
- Compound legs (`K103N_to_*`) source from the K103N MD start.pdb.
- The 3 K103N-compounds were rebuilt exactly this way 2026-08-15 (126 tasks em→…→switch, **0 failures**).

**`combine_neq`: run WITHOUT `--force`.** `--force` re-analyzes every leg, but the 15 existing genotypes'
raw dgdl is deleted → it would fail on missing inputs. No `--force` = uses cached `analysis.json` and
auto-analyzes only legs that *lack* one (the freshly-run ones). Verify first with the check in
`data-safety-and-sync` memory / this session's transcript.

---

## SEM: the lever is EQUILIBRATION, not switches or snapshots (2026-08-15 — important)

Goal (human's bar): **all SEM < 1 kcal/mol.** Decomposition of each genotype's SEM (per-rep `bar_dg`
spread `σ_between` vs per-rep `bar_err_boot` `σ_within`) shows the SEM is dominated by **between-replicate
variance everywhere**: σ_between 0.2–4.7 vs σ_within 0.07–0.96, ratio **2.7–20**. Therefore:

- **More snapshots** (`NEQ_SNAPSHOTS`, tightens per-rep BAR) and **longer/slower switches** (reduces
  dissipation/overlap bias) both target the *within-rep* error, which is already tiny — **they will NOT
  reduce the panel SEM.** Making each rep more precise just makes 3 reps confidently disagree.
- The SEM is the reps DISAGREEING — under-equilibrated endpoints (the mutant / λ=1 state is the noisy
  half). Efficient levers, in order: **(1) longer endpoint equilibration** (`NEQ_EQUIL_NS`, currently
  5 ns, env-overridable) to shrink σ_between at the root; **(2) more reps** (SEM = σ_between/√n) — but
  only after (1), because for K103N-family + G190S (σ_between 2.2–4.7) hitting <1 by reps alone needs
  n ≈ 20+.
- By tier: neutrals (σ_between ~1.2–1.5, SEM 0.7–0.87) ≈ there with a modest equil bump + n = 6.

**CORRECTION (2026-08-15, later session) — two errors in the paragraph above.**

1. **G190S is NOT a high-variance case and is a bad test case.** Its σ_between ~4 is *per-phase*;
   holo and apo move together across reps and cancel in the double difference, so **σ_DDG = 0.59**
   — already under the < 1 target. Rank genotypes by **σ_DDG** (the spread of the three per-rep
   ΔΔG values, where SEM = σ_DDG/√n and n needed = σ_DDG²), never by per-phase σ_between: across
   the 18 legs the two are nearly uncorrelated, because holo/apo cancellation ranges from
   near-total (G190S, V106M, K103N+M230L) to actively additive (K103N, V106I, Y181C).
2. **"Switch length will not help" does not hold for the charge leg.** That was measured on
   neutral legs dissipating 1–4 kcal/mol (V106A, F227C), which had nothing to gain. `wt_to_K103N`
   dissipates **14–19 kcal/mol** with forward/reverse work distributions **4.2–5.1σ** apart — the
   regime where BAR is biased, not just noisy — and its per-rep ΔG tracks its per-rep dissipation
   ~1:1 (hyst 14.19/18.54/19.36 → bar_dg −75.71/−77.99/−79.79), which is the signature of varying
   residual bias. Longer switches attack exactly that.

**The whole charge-family SEM problem is ONE leg.** All four K103N genotypes inherit their ~2.2 SEM
from `wt_to_K103N`; `K103N_to_K103N_M230L`'s own σ_DDG is 0.59. Fix that leg → four of the five
worst panel points improve at once.

Test it on G190E (step 11 of [`RUNBOOK_G190E.md`](RUNBOOK_G190E.md)): re-run **only its switch
stage** at 500 ps reusing its own equil/extract snapshots — identical endpoint ensembles, so it
isolates switch length from equilibration in a way the V106A test never did. If σ_DDG drops
materially, rebuild `wt_to_K103N` the same way; if only hysteresis drops, the lever really is
`NEQ_EQUIL_NS` + more reps. With n = 3 every σ here carries ~±40%, so only large changes count.

Required reps for SEM < 1 at current σ_DDG: V106M/G190A/G190S/Y318F/Y188L/V106A already there;
F227C/V106I/A98G+F227C 3–5; V106A+F227L/L234I/P225H, Y181C, V106I+F227C 4–5; **K103N ×4 need 15–17.**

---

## IN FLIGHT (submitted 2026-08-15) — two batches, 60 GPU array elements

| batch | manifest | equil | extract | switch |
|---|---|---|---|---|
| **G190E, 100 ps** (completes the panel) | `neq_g190e_manifest.csv` | ~~39256102~~ | ~~39256103~~ | ~~39256104~~ |
| **K103N, 500 ps** (the SEM experiment) | `neq_k103n500_manifest.csv` | `39257456` done | `39257460` done | `39257464` (24) running |
| **apo MD** y181c / v106i / v106a_l234i | (globbed, not manifest-driven) | — | — | `39258048`–`39258062`, 9 jobs |

Both FEP batches passed `EM (6/6 ok)`. Runbook: [`RUNBOOK_G190E.md`](RUNBOOK_G190E.md).

### G190E was CANCELLED and must be re-submitted — task-id file collision

`scancel 39256102 39256103 39256104` was run. Cause: `submit_p0_neq.sh` writes task
ids to a path with **no batch qualifier**, and array elements read it **at runtime**
— so submitting K103N overwrote the file under G190E's still-queued elements.
G190E equil 0–6 completed (pre-overwrite); 7–11 died in ~10 s running K103N's task
ids against G190E's manifest. **Full explanation and the `TASK_ID_FILE` fix:
[`OPERATIONS.md`](OPERATIONS.md) §4c.** Re-submit with per-batch task-id files —
stages are idempotent so completed units skip.

**Also always pass `--manifest` explicitly** to `audit_neq_panel.py`; a stale
exported `MANIFEST` silently audits the wrong batch. The tell is the switch count:
**G190E = 12, K103N-500 = 24**.

### `$SCRATCH` was DEGRADED (hardware issue) as of the cancel

Login banner reported `- $SCRATCH: Hardware issue`. Filesystem-heavy commands hung
unkillably (D-state). **Do not re-submit until https://status.sherlock.stanford.edu
is clear**, and when it is, verify anything that was mid-write — K103N switch
elements 0–5 and all 9 apo jobs were running through it. A truncated `dgdl.xvg`
parses into a wrong work value rather than erroring. See [`OPERATIONS.md`](OPERATIONS.md) §4d.

When each finishes, analyse **on Sherlock** (rsync excludes `dgdl.xvg`, so the Mac
has nothing to analyse until `analysis.json` exists), and **without `--force`**:

```bash
python3 scripts/fep_pmx/combine_neq.py --targets G190E --replicates 3
python3 scripts/fep_pmx/combine_neq.py --targets K103N --replicates 3
python3 scripts/fep_pmx/qc_neq.py --legs wt_to_G190E wt_to_K103N --replicates 3
```

**What to look for in K103N-500** (baseline table in
`results/analysis/fep_pmx/_archive/wt_to_K103N_100ps_2026-08-15/README.md`):
σ_DDG 3.80 → under ~1.7 puts SEM < 1 and clears four genotypes at once; hysteresis
14–19 → ~3–4; separation 4.2–5.1σ → ~2σ. G190E is expected to come in at SEM > 1 —
that is the K103N-matched baseline, not a failure.

## Replicate work (the other lever) — blocked, partially

Extra reps always work (SEM = σ_DDG/√n) but the leg resolver needs
`results/md_runs/<label>/rep_{n:02d}/` and **only three MD replicates exist**.
`scripts/fep_pmx/seed_extra_replicates.py` creates a rep_04 seed from a run's
`*_md_final.pdb` (validated: atom count, PBC-split, ligand still bound). Dry run by
default.

- **Seedable now:** `wt_to_V106A` +1, `V106A_to_V106A_P225H` +1 (16 elements) →
  clears V106A+P225H 1.13 → 0.98 and V106A+F227L 1.00 → 0.97. New seeds live under
  `results/md_runs/` (rsync territory, not git) so they must be pushed to Sherlock.
- **Blocked:** `wt_to_Y181C`, `wt_to_V106I`, `V106A_to_V106A_L234I` — their apo MD
  never ran (`results/md_runs/apo/{y181c,v106i,v106a_l234i}/rep_0*/` hold only
  `assets/`, on **both** Mac and Sherlock), so there is no 4th apo conformation.
  Options: run ~180 ns of apo MD; derive from a stripped holo frame (introduces a
  within-leg treatment asymmetry — PLAN §6.3); or reuse the apo start with fresh
  velocities (only partially independent, so SEM falls slower than √n).
  **Deferred pending K103N-500** — these genotypes sit at 1.06–1.16, the least
  valuable per GPU-hour, and if 500 ps works it may apply here too.

Full cost model (validated: reproduces every panel SEM to 0.000) —
`SEM_genotype = sqrt(Σ_legs σ_leg²/n_leg)`. All genotypes < 1 by replicates alone
costs **29 leg-reps, of which 19 are `wt_to_K103N`** (≈152 GPU elements vs ≈36 for
the 500 ps route). One leg-rep = 8 GPU elements at 100 ps, 12 at 500 ps.

## Next (human's intent at session end 2026-08-15)

1. **Run G190E** (last panel genotype; charge leg) at the **K103N-matched protocol** — 100 ps
   switches / 100 snapshots / 5 ns equil / 3 reps, raw non-neutral box + analytical charge
   correction. Full runbook: [`RUNBOOK_G190E.md`](RUNBOOK_G190E.md).
   **(It is G190E, NOT G190S — a session confusion the human corrected.)**
   Two fixes were needed first, both done on the Mac:
   - `config.py`: `wt_to_G190E` removed from `_BASE_LONG_SWITCH_LEGS` (was 500 ps — a speculative
     entry from the co-alchemical era; K103N ran at 100 ps). **Must reach Sherlock via `git pull`
     before prep**, or G190E rebuilds at 500 ps.
   - Stale co-alchemical G190E artifacts (ΔΔG 2.50 ± 1.41, apo `bar_err` 35.49) archived to
     `results/analysis/fep_pmx/_archive/wt_to_G190E_coalchemical_500ps_2026-08-15/`. They were
     git-tracked, so `combine_neq` **without `--force`** (now mandatory) would have reused them
     and put a non-convergent G190E on the panel silently.

   Expect **SEM > 1** at this configuration — the K103N family sits at 2.19–2.32 under exactly this
   protocol. That is the honest baseline and the input to step 11, not a failed run.
2. **Panel-wide SEM < 1.** Rank by σ_DDG (see the CORRECTION above). Test the switch-length lever
   on G190E via runbook step 11 (switch stage only, reusing its own snapshots) before committing
   GPU to a `wt_to_K103N` rebuild; if that comes back flat, fall back to `NEQ_EQUIL_NS` + more reps.
3. **Manuscript**: `manuscript/DorDRM-FEP-08-05-26.docx` (untracked working drafts in `manuscript/`).
   The second-wave analysis is done: protocol figures, `modern_md_suite`, `occupancy_stats` (FWER +
   Welch), the DCD-fingerprint MD-timing correction (`src/analysis/md_timing.py`).

---

## Repo sync + pipeline mechanics

- **All three repos (Mac, Sherlock, GitHub) at commit `acb36fc`.** git is the hub (code + light results);
  rsync only for heavy MD artifacts. `.gitignore` whitelists `legs/**/analysis/**` + manifests; heavy
  switch/equil data stays gitignored on Sherlock scratch. Sherlock pushes over SSH.
- **Uncommitted on the Mac** (session ended mid-commit): a widened-panel figure edit
  (`combine_neq.py` figsize 6×5 → 9.5×6.5 + regenerated `panel_ddg_vs_experiment.png`/tiers). Decide
  whether to keep it next session.
- Stages `em→equil→extract→switch` via `submit_p0_neq_pipeline.sh`; em/extract = CPU, equil/switch = GPU;
  whole-array `afterok` (one failed equil element stalls the batch — OPERATIONS §1). Equil ~1h50m/task at
  5 ns; switch ~3–9 h/task (GPU heterogeneity). **Read `OPERATIONS.md` before acting.**
