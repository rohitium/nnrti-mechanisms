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
- By tier: neutrals (σ_between ~1.2–1.5, SEM 0.7–0.87) ≈ there with a modest equil bump + n = 6; the
  charge family + **G190S (σ_between 4.3, a clean neutral outlier — good cheap test case)** need longer
  equilibration (mandatory). This corrects earlier notes that leaned on more-reps / 500 ps switches.

---

## Next (human's intent at session end 2026-08-15)

1. **Run G190E** (last panel genotype; charge leg) **with longer equilibration + more reps** so its SEM
   lands < 1 — at 5 ns / 3-rep it would be > 1 (unacceptable per the human). Build `wt_to_G190E`
   (source = WT `_start.pdb`), then prepare with a longer `NEQ_EQUIL_NS` and n ≈ 6, then submit.
   **(It is G190E, NOT G190S — a session confusion the human corrected.)**
2. **Panel-wide SEM < 1** via longer equilibration (+ modest extra reps), NOT more switches/snapshots.
   Validate the equilibration lever on a high-σ_between genotype first, then scale.
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
