# pmx NEQ FEP — current state

**Living snapshot. Update this when state changes.** A fresh agent should read this first, then
[`OPERATIONS.md`](OPERATIONS.md) for how to act. Last meaningful update: **2026-08-03**.

Sherlock repo: `/scratch/users/rsatija/nnrti-mechanisms-git` (account `rshafer`, QOS `long,normal`).
The human runs everything on Sherlock and pastes output — agents cannot reach it.

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

## Roadmap (P-stages)

- **P0 — pilot (V106A, Y188L): DONE & validated.** ✅
- **P1 — single mutations → Spearman ρ vs fold (a scientific *readout*, not a gate):** in progress.
  - `P1_NEUTRAL_LEGS` (config.py): F227C, G190A, V106I, V106M, Y181C, Y318F.
  - With P0's V106A + Y188L that's the **8 single legs** the current pipeline can run.
  - `P1_CHARGE_LEGS` (K103N, G190E): **excluded from this pipeline** — need the co-alchemical ion /
    double-box charge protocol (PLAN §6.2), not yet implemented. They are *not* part of the ρ readout.
- **P2 — compound genotypes → additivity check:** not started.
- **P3 — full manuscript table + experimental correlation:** not started.

Panels are submitted in **3-leg batches** to stay under the GPU QOS cap (equil 36 + switch 36 = 72 ≤ 100).

---

## Results so far (all 5 ns-seed, the current pipeline)

**P0 (validated pilot):**

| genotype | ΔΔG_bind (kcal/mol) | exp. fold |
|---|---|---|
| V106A | +1.76 ± 0.51 | 9.6 |
| Y188L | +4.52 ± 0.49 | 149 |

Converged: running-BAR plateau, BAR/CGI/Jarzynski agree within ~0.4, replicate SEM < 1, switch-length
invariant (100 ps ≈ 500 ps).

**P1a (complete, analyzed 2026-08-03):**

| genotype | ΔΔG_bind (kcal/mol) | exp. fold | notes |
|---|---|---|---|
| G190A | +0.27 ± 0.17 | 2.7 | tight, healthy overlap (0.14–0.77) |
| F227C | −0.21 ± 0.82 | (none) | **least reliable** — holo reps 2&3 have 0.00 overlap; no fold to compare |
| V106I | +2.27 ± 0.74 | 1.1 | reproducible (rep sd 0.54), overlap in the validated regime |

**The V106I finding (not a bug — a result):** V106I gives a confident, reproducible ΔΔG (~+2.3) whose
overlap is no worse than the trusted P0 legs, so it is *not* a switching-overlap artifact. The method
predicts V106I and V106A within ~0.5 kcal of each other (within error) — i.e. it **cannot resolve Ile
from Ala at position 106** — while experiment separates them strongly (Ile 1.1× vs Ala 9.6×). Two live,
both-interesting explanations: (a) the ΔΔG estimate could still shift under longer endpoint sampling
(the seeding test probes this — overlap measures switch dissipation, not endpoint-ensemble adequacy),
or (b) V106I's near-neutral phenotype is genuinely not binding-mediated. V106I is therefore the **best
target** for the endpoint-seeding test (better than V106A, which is well-behaved).

---

## In flight (as of 2026-08-03)

**Batch P1b — legs V106M, Y181C, Y318F** (manifest `results/analysis/fep_pmx/neq_p1b_manifest.csv`).
Launched after a **GPU-free pre-flight**: ran the CPU `em` stage alone → `EM (18/18 ok)`, confirming
topology/structure/mdp inputs are valid, *then* chained the GPU stages. Live job chain:

- em `37405830` (done, CPU) → equil `37405873` (gpu) → extract `37405881` (normal) → switch `37405887` (gpu).
- Watch for `SWITCH (36/36 ok)` via the audit (OPERATIONS §2).

**P1a is complete** — see results above. Its recovery from a bad-GPU-node stall (node `sh03-12n12`,
`GPU is lost`) is the canonical worked example in OPERATIONS §4.

**P2 legs PREPARED (staged, not yet submitted), 2026-08-04.** 6 neutral compound/extra legs —
`wt_to_G190S`, `F227C_to_A98G_F227C`, `V106A_to_V106A_F227L`, `V106A_to_V106A_L234I`,
`V106A_to_V106A_P225H`, `V106I_to_V106I_F227C` — now have hybrids + solvated systems + neq inputs
(`results/analysis/fep_pmx/neq_p2_manifest.csv`, 6 legs). Submit in 3-leg batches (em pre-flight +
em→…→switch) once P1b frees GPUs. This brings the runnable panel to **13 genotypes** (8 singles with
fold + 5 neutral compounds). See OPERATIONS §6 for how they were built (the `_start.pdb` push / openmm
gotchas). Env note: the Sherlock pmx venv now has **openmm 8.1.1** and was pinned back to **numpy<2**.

**Proline mutations FIXED (2026-08-04, commit c21415b).** `V106A_to_V106A_P225H` (Pro→His, fold 153, a
positive control) now builds. Root cause was proline HG/HD naming (`HG2/HG3, HD2/HD3` vs pmx's
`HG1/HG2, HD1/HD2`); `normalize_openmm_for_pmx` now renames it PRO-scoped + idempotent, verified through
mutate **and** pdb2gmx. No legs are deferred for proline anymore.

**Charge protocol: co-alchemical ion ABANDONED — analytical correction next (2026-08-05).** The explicit co-alchemical ion (`coalchemical_ion.py`) ran end-to-end but does **not converge**: charge legs dissipate ~20–26 kcal/mol vs ~1–3 for neutral (~10×) → near-zero forward/reverse overlap, SEM ~1.4, BAR–Jarz disagreement up to ~3.7. Diagnostic pinned it to the shared Cl⁻ decoupling (dissipation ~constant across K103N ΔG≈9 and G190E ΔG≈36). **Decision:** skip the co-alchemical ion; run charge legs raw (non-neutral box) + apply the Rocklin/Hünenberger analytical net-charge correction (zero added perturbation). K103N, G190E, and the 3 K103N-compound genotypes await this. `coalchemical_ion.py` kept but not used for production. See OPERATIONS §7 for the lesson. (First co-alchemical numbers, not trusted: K103N −0.72±1.38, G190E +2.50±1.41.)

**F227C 500 ps switch rerun IN FLIGHT (2026-08-08, job `38228510`, array 0-23%16).** F227C's
100 ps holo reverse switches dissipate erratically (fwd/rev gap 2.3→7.7→11.6 across reps, 0.00
overlap on rep 3) — regrowing the Phe ring into the collapsed Cys pocket catches variable clashes.
Its ~2 kcal per-rep ΔG spread is the source of the imprecision in **A98G+F227C** (fold 93), which
is built on this exact leg (= WT→F227C + F227C→A98G+F227C); the A98G second leg is tight (holo SEM
0.22, apo 0.16). **Correction (2026-08-10):** V106I+F227C (fold 105) does NOT use this leg — it
reaches F227C via a separate `V106I_to_V106I_F227C` leg (same Phe→Cys, V106I background) that was
NOT bumped and likely needs its own 500 ps rerun; check its QC in the panel combine. So this rerun
tightens **two** genotypes (F227C, A98G+F227C), not three. `wt_to_F227C` was added to
`LONG_SWITCH_LEGS` (commit cf7bd45) and its **switch stage
only** was rerun at 500 ps from an **isolated** manifest (`neq_f227c500_manifest.csv`) with its own
`TASK_ID_FILE` (so it can't clobber the concurrently-pending P2b switch array). equil + snapshots
reused; only `legs/wt_to_F227C/*/rep_*/neq/{switches,analysis}` were deleted before the rerun.
**When 38228510 finishes** (24 tasks × ~15 h): re-run combine/qc for all three dependents and pull:
```bash
source scripts/sherlock/activate_pmx_env.sh
python3 scripts/fep_pmx/combine_neq.py --targets F227C --replicates 3
python3 scripts/fep_pmx/combine_neq.py --targets A98G+F227C V106I+F227C --replicates 3
python3 scripts/fep_pmx/qc_neq.py --targets F227C A98G+F227C V106I+F227C --replicates 3
# then, on the Mac: SHERLOCK_USER=rsatija bash scripts/rsync_fep_pmx.sh pull
```

**F227C 500 ps RESULT — noise is endpoint-limited, not switch-limited (2026-08-10).** The rerun
finished (24/24) and was analysed. Verdict: **switch length is not the lever.** Per-rep ΔG at 500 ps
is essentially identical to 100 ps (holo −0.18/−1.72/−2.15 vs −0.1/−1.8/−2.2); dissipation gaps shrank
~20 % (holo r3 11.6→9.2) but BAR already corrects for that, so the estimate didn't move — clean
switch-length invariance. The imprecision is **across-rep endpoint scatter**: F227C leg ΔΔG per rep =
+0.51, +0.51, −1.93 (reps 1&2 agree; rep 3 is a ~2.4 kcal outlier whose whole work distribution is
shifted/more dissipative → it equilibrated into a different pocket conformation). QC panel-wide shows
overlap is **decoupled from precision** (28/30 units low_overlap, yet SEM spans 0.16–1.61; V106M all
6 reps low_overlap but SEM 0.16, fully converged) — trust replicate SEM + BAR–Jarz agreement, NOT the
overlap flag. Panel ρ = **0.374 (n=13)**, up from 0.188. V106M (+6.10, fold 3.4, tight) is a genuine
converged overprediction and the main ρ drag. **Consequences:** (a) do NOT bump `V106I_to_V106I_F227C`
to 500 ps — same Phe→Cys, same endpoint problem, switch length won't help; (b) the levers for the
F227C family are MORE REPLICATES (SEM∝1/√n, cheapest at 100 ps since ≈500 ps) and ultimately the
endpoint-seeding test (the root fix). Mechanistic echo: M66I/lenacapavir capsid paper (side-chain
reorganization free energy) — see manuscript/fep_resistance_references.md §5d.

**ALL 18 GENOTYPES NOW PREPARED (2026-08-04).** The 3 K103N-compound 2nd legs
(`K103N_to_K103N_M230L`, `K103N_to_K103N_P225H`, `K103N_to_L100I_K103N` — manifest
`neq_kcompound_manifest.csv`) are staged (hybrids + systems + neq), so nothing in the manuscript panel
remains to *prepare* — only to run + analyze. Run-status: P0/P1a done; P1b in switch tail; charge legs
(K103N, G190E) running + validating; P2 (6 legs) and the K-compound legs (3 legs) staged, awaiting GPU.
Submit remaining in 3-leg batches as GPU frees; the K-compound genotypes' ΔΔG also needs the K103N
charge leg (running).

---

## Sequenced next steps (order matters — this is the human's explicit decision)

1. **Finish P1b, then compute the full 8-leg ρ** — the deliverable:
   ```bash
   # on Sherlock after P1b SWITCH 36/36 (dgdl live there; the Mac rsync excludes them):
   python3 scripts/fep_pmx/combine_neq.py --targets V106M Y181C Y318F --replicates 3
   python3 scripts/fep_pmx/combine_neq.py --targets V106A Y188L F227C G190A V106I V106M Y181C Y318F --replicates 3
   ```
   This is the complete 5 ns-seed panel + Spearman ρ vs fold. **Do this before touching endpoint seeding.**

2. **THEN the endpoint-sampling experiment** (deferred until step 1 is in hand). First extend WT apo MD
   10 ns → 100 ns (`./scripts/fep_pmx/submit_wt_apo_md.sh` — resumes each rep's checkpoint to 100 ns
   total, verified in `src/md/openmm/md_protocol.py`; the wrapper bakes in `MD_FORCE_RERUN=1` +
   `SKIP_IF_AT_TARGET=1`; each 12 h SLURM job advances the checkpoint, rerun until at target). Then
   re-seed switch snapshots from decorrelated 100 ns frames and test whether ΔΔG moves — **on V106I**
   (and V106A as control). Purpose is *confidence/robustness of the estimate*, not matching fold. If it
   moves, escalate to enhanced endpoint sampling (REST2/HREX/metadynamics). Manuscript §8.2.

3. **Manuscript:** `manuscript/DorDRM-FEP-07-30-26.docx` — update with the P1 panel + ρ once step 1
   lands, framed per "What success means" above (lead with reproducibility/convergence; treat the fold
   comparison as a hypothesis test, noting mechanism may differ from binding). A draft collaborator
   email still needs a fact-check fix (separate the switch-length argument from the µs-equilibration
   argument; keep MM/GBSA comparison claims supported).
