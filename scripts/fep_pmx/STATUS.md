# pmx NEQ FEP — current state

**Living snapshot. Update this when state changes.** A fresh agent should read this first, then
[`OPERATIONS.md`](OPERATIONS.md) for how to act. Last meaningful update: **2026-08-02**.

Sherlock repo: `/scratch/users/rsatija/nnrti-mechanisms-git` (account `rshafer`, QOS `long,normal`).
The human runs everything on Sherlock and pastes output — agents cannot reach it.

---

## Roadmap (P-stages)

- **P0 — pilot (V106A, Y188L): DONE & validated.** ✅
- **P1 — single mutations → Spearman ranking gate:** in progress.
  - `P1_NEUTRAL_LEGS` (config.py): F227C, G190A, V106I, V106M, Y181C, Y318F.
  - With P0's V106A + Y188L that's the **8 single legs** for the ρ gate.
  - `P1_CHARGE_LEGS` (K103N, G190E): **deferred** — need the co-alchemical ion / double-box charge
    protocol (PLAN §6.2), not yet implemented.
- **P2 — compound genotypes → additivity check:** not started.
- **P3 — full manuscript table + experimental correlation:** not started.

Panels are submitted in **3-leg batches** to stay under the GPU QOS cap (equil 36 + switch 36 = 72 ≤ 100).

---

## P0 results (validated pilot)

| genotype | ΔΔG_bind (kcal/mol) | exp. fold | signs/ranking |
|---|---|---|---|
| V106A | **+1.76 ± 0.51** | 9.6 | positive (resistance), ranked below Y188L ✓ |
| Y188L | **+4.52 ± 0.49** | 149 | positive, correctly largest ✓ |

Converged: running-BAR plateau, BAR/CGI/Jarzynski agree within ~0.4, replicate SEM < 1, switch-length
invariant (100 ps ≈ 500 ps). Marginal forward/reverse overlap inflates per-leg error but does not bias
ΔΔG. (Earlier README quotes V106A +1.69 ± 0.70 from the very first run; +1.76 ± 0.51 is the current number.)

---

## In flight (as of 2026-08-02)

**Batch P1a — legs F227C, G190A, V106I** (manifest `results/analysis/fep_pmx/neq_p1a_manifest.csv`):
just recovered from a bad-GPU-node incident and re-submitted. Live job chain:

- `37288420` equil (gpu) → `37288421` extract (normal) → `37288424` switch (gpu), chained `afterok`.
- Watch for `SWITCH (36/36 ok)` via the audit (OPERATIONS §2).

What happened (so the next agent doesn't re-investigate): node **sh03-12n12** had a lost/dead GPU;
equil array elements 23/24/25 (V106I holo rep1 λ0/λ1, G190A apo rep3 λ1) failed there with
`0 detected device(s)` / `GPU is lost`. That stalled all 36 extract + 36 switch via whole-array
`afterok`. Fixed with the standard §4 recovery: cancelled the old extract/switch (`36839194`,
`36839195`), re-ran equil (33 skipped, 3 recomputed), re-chained. The bad node self-invalidated
(`inval`), so no manual exclusion was needed. **This is the canonical example of OPERATIONS §4.**

**Batch P1b (and beyond) — HELD.** Remaining neutral legs (V106M, Y181C, Y318F) not yet submitted.
Hold until (a) P1a's GPUs free up, and (b) the seeding decision below is resolved.

---

## Open decisions / pending work

1. **FEP endpoint seeding (agreed, not yet implemented).** Current pipeline seeds FEP endpoints from
   `*_start.pdb` + a 5 ns hybrid re-equilibration — **not** from the 100 ns plain-MD end frames the
   collaborators asked about. Plan: seed switch snapshots from decorrelated frames of the long MD
   trajectories instead, and run a **sensitivity test** (5 ns-seed vs 100 ns-seed on V106A apo) — same
   logic as the switch-length test, applied to endpoint sampling. If ΔΔG moves, escalate to enhanced
   sampling (REST2/HREX/metadynamics). See manuscript §8.2.

2. **Apo WT 100 ns extension (scripted, ready to launch).** WT apo MD only ran to **10 ns** (5M steps ×
   3 reps; MM/GBSA didn't need more), while holo ran to 100 ns. WT apo is the shared endpoint for every
   single-mutation leg, so to enable 100 ns-seeded FEP we extend it first. Wrapper committed:
   [`submit_wt_apo_md.sh`](submit_wt_apo_md.sh) — resumes each rep's checkpoint to 100 ns total (not
   +100 ns; verified in `src/md/openmm/md_protocol.py`). Needs `MD_FORCE_RERUN=1` (the 10 ns runs are
   `status=ok` and would otherwise skip) + `SKIP_IF_AT_TARGET=1`, both baked into the wrapper. **Launch
   after Batch P1a finishes.** Each SLURM job is 12 h and 90 ns won't fit in one — rerun the same
   command after each batch; it resumes and skips reps already at 100 ns.

3. **Manuscript:** `manuscript/DorDRM-FEP-07-30-26.docx` (collaborator update; MM/GBSA→FEP pivot
   rationale, protocol, P0 results, limitations, next steps). To update with real P1 numbers once the
   ranking gate completes. A draft collaborator email still needs a fact-check fix (separate the
   switch-length argument from the µs-equilibration argument; keep MM/GBSA comparison claims supported).

---

## Definition of done for P1

`combine_neq --targets V106A Y188L F227C G190A V106I V106M Y181C Y318F --replicates 3` yields a
Spearman ρ vs experimental fold across the 8 single mutations, with QC (`qc_neq`) passing. That ρ is
the P1 ranking gate; clearing it unlocks P2.
