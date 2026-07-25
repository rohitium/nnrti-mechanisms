# Plan: truncated-site Doravirine mutation FEP

Status: **design only** (no implementation in this folder yet).  
Branch target: `truncated-fep` off `main`.

---

## 1. Scientific goal

Estimate **relative Doravirine (DOR) resistance** for manuscript RT mutations:

\[
\Delta\Delta G_\text{bind} = \Delta G_\text{mut}^\text{holo} - \Delta G_\text{mut}^\text{apo}
\]

where \(\Delta G_\text{mut}\) is the alchemical free energy to transform the **wild-type side
chain into the mutant side chain** at one RT site, in either the DOR-bound complex (**holo**)
or the ligand-free pocket (**apo**).

Compare \(\Delta\Delta G_\text{bind}\) (kcal/mol) to experimental **DOR fold-change in
susceptibility** (manuscript supplementary table / `dor_fold_reduction` in analysis CSVs).
We do **not** require exact MCPRO reproduction; we require a **Jorgensen-shaped** cycle that
is cheap enough to run for the full mutation panel.

### Observable for the manuscript

Primary metric per target genotype:

| Quantity | Definition |
| --- | --- |
| \(\Delta G_\text{mut}^\text{holo}\) | WT → mutant alchemical \(\Delta G\) in DOR-bound truncated pocket |
| \(\Delta G_\text{mut}^\text{apo}\) | Same mutation in apo truncated pocket |
| \(\Delta\Delta G_\text{bind}\) | holo − apo (positive ⇒ weaker DOR binding upon mutation) |
| Experimental reference | \(\mathrm{RT}\ln(\mathrm{fold\_change})\) or ranked order vs WT |

Compound mutants (e.g. `L100I+K103N`) are **sums of sequential single-mutation legs**, reusing
`TargetPlan` / `MutationLeg` logic from `scripts/fep_jorgensen/mutations.py`.

---

## 2. Why this protocol (lessons from the Perses pilot)

The merged `jorgensen-fep` work tested **full-protein Perses hybrids + fixed-λ MD + MBAR**
(`scripts/fep_jorgensen/`). V106A holo/apo windows finished on Sherlock; analysis showed:

| Issue | Implication for truncated plan |
| --- | --- |
| ~200k-atom PME systems | Too expensive; wrong geometry for Jorgensen comparison |
| λ ≥ 0.8 endpoint clashes | Ghost methyls at high λ → MBAR re-eval at λ=0 blows up; **not fixed by longer sampling** |
| Own-state energy “drift” | Barostat artifact; ignore; use **cross-state spread** for QC |
| \(\Delta G_\text{mut} \sim -18\) kcal/mol | Non-physical magnitudes; do not trust full-path MBAR |
| Auth **V106** = PDB **103** | All prep must map manuscript residue IDs → PDB (`prepare_backend.json` pattern) |

**Conclusion:** match Jorgensen’s ** reduced model** (truncated pocket, fixed backbone, local
flexibility, small solvent), not Perses full-protein AREX.

### Literature anchors (what we copy vs what we omit)

| Source | What to copy | What we omit (no license / too heavy) |
| --- | --- | --- |
| [Rizzo/Jorgensen 2000 (JACS)](https://doi.org/10.1021/ja003113r) | Thermodynamic cycle; side-chain mutation in bound vs unbound model; inhibitor-relative reporting | IMPACT, MCPRO, CM1P-OPLS, ε=4r, MC sampling |
| [Smith et al. 2007 (BMCL)](https://doi.org/10.1016/j.bmcl.2007.12.033) | ~123-residue binding site; 22 Å water cap; fixed backbone; 10 λ windows; holo + apo | MCPRO, OPLS charges on ligand |
| Perses / Zhang JCTC 2023 | Endpoint-aware diagnostics (cross-state spread, forward/reverse) | Full-protein AREX as default engine |

Exact MCPRO invariants remain documented in `scripts/fep_jorgensen/exact_protocol.py` for
boundary checking only.

---

## 3. System model

### 3.1 Starting structures (reuse manuscript MD assets)

| Phase | Source (existing) |
| --- | --- |
| **Holo** | `results/md_runs/{genotype}/rep_01/assets/*_md_rep01_start.pdb` (DOR-bound) |
| **Apo** | `results/md_runs/apo/{genotype}/rep_01/assets/*_apo_md_rep01_start.pdb` |

Genotypes follow `MANUSCRIPT_PLANS` in `scripts/fep_jorgensen/mutations.py` (19 targets,
19 unique single-mutation legs for WT→single mutants, plus compound legs).

**Do not** rebuild from raw PDB if equilibrated rep_01 assets exist.

### 3.2 Truncated pocket

Extract residues with any atom within **15 Å** of DOR heavy atoms (Jorgensen used ~15 Å around
the NNRTI site; Smith et al. ~123 residues — we derive the set per structure, not hard-code
2000 numbering).

Output per (genotype, phase):

- `truncated.pdb` — pocket + DOR (holo) or pocket only (apo)
- `residue_map.json` — manuscript auth ID, PDB resSeq, chain, mutation site flag
- `flexible_residues.txt` — side chains allowed to move (see §3.3)

### 3.3 Flexibility and restraints (Jorgensen-shaped)

Map Jorgensen 2000 partitioning onto our pocket:

| Region | Treatment in OpenMM |
| --- | --- |
| **Protein backbone** (N, C, Cα, carbonyl O) | **Fixed** (`CustomExternalForce` or `ZeroMagnitudeForce` on backbone atoms) |
| **Side chains within 10 Å of DOR geometric center** | **Fully flexible** (rotamer + χ sampling via MD) |
| **Side chains 10–12 Å shell** | **Harmonic positional restraints** (k ≈ 10–100 kcal/mol/Å²; tune in pilot) |
| **Pocket residues beyond 12 Å** | Fixed (already excluded by truncation) |
| **DOR** | Fully flexible (translation/rotation/χ); holo only |

Apo uses the same protein restraints; no ligand.

### 3.4 Solvent

**Primary (pilot):** explicit **spherical cap** ~22 Å radius centered on the mutation site /
ligand centroid, TIP3P, **no PME** (cutoff + reaction field or GB neck — decide in pilot).

**Rationale:** matches MCPRO water-cap count (~850 waters) and avoids full-box barostat noise
that dominated the Perses pilot.

**Fallback:** implicit solvent (OBC2/GBn2) for apo leg only if cap proves unstable — document
if used; holo/apo must stay consistent within a given run series.

### 3.5 Force field

**Pragmatic choice (match existing MD):**

- Protein: AMBER ff14SB  
- Water: TIP3P  
- DOR: OpenFF 2.0.0 (same as `approx_protocol.py`)

This is **not** OPLS/CM1P. Document as an approximation vs Jorgensen; prioritize consistency
with the manuscript’s existing OpenMM MD structures.

---

## 4. Alchemical transformation

### 4.1 Strategy

**Single-topology nonbonded scaling** on mutation-site side-chain atoms (reuse design from
`scripts/fep_jorgensen/alchemical.py` + `worker.py` `nonbonded-scaling` path).

Avoid Perses hybrid topologies (delete/insert atoms) in v1 — they caused endpoint pathology.

For each leg (e.g. WT→V106A):

1. Identify side-chain atoms of **old** residue at the mutation site in the truncated PDB.
2. Scale their LJ/charges by λ: fully on at λ=0 (WT), fully off at λ=1 (mutant endpoint).
3. **Mutant endpoint geometry:** mutate PDB in place (delete side-chain atoms beyond Cβ for
   Ala, etc.) for **analysis endpoints only**; for alchemical path, use **single hybrid
   coordinate set** with annihilation (Jorgensen mutates large→small).

Charge-changing mutations (none in initial pilot except possibly Y181C, K103N context): flag
for longer windows / separate validation.

### 4.2 λ schedule

| Parameter | Value |
| --- | --- |
| Windows | **11** (λ = 0.0, 0.1, …, 1.0) — same spacing as Smith overlap sampling |
| Transform direction | **Annihilate WT side chain** (large → small), matching MCPRO “mutate larger to smaller” |

Optional v2: 14-window double-wide spacing near λ=0 and λ=1 if variance high.

### 4.3 Phases (holo + apo)

Each mutation leg requires **two independent FEP strings**:

| Phase | System |
| --- | --- |
| holo | truncated pocket + DOR |
| apo | truncated pocket, no ligand |

\[
\Delta\Delta G_\text{bind} = \Delta G_\text{mut}^\text{holo} - \Delta G_\text{mut}^\text{apo}
\]

---

## 5. Sampling protocol

Target: **≤ ~1 GPU-day per mutation leg** (both phases), not hundreds of ns × 200k atoms.

| Step | Setting (initial) |
| --- | --- |
| Integrator | LangevinMiddle, 2 fs, H-bonds constrained |
| Temperature | 300 K |
| Minimization | 500 steps (STeepestDescent or L-BFGS) |
| Equilibration | 50 ps restrained (match Jorgensen MD equil order-of-magnitude) |
| Production per window | **2–5 ns** (pilot: 2 ns; extend if MBAR uncertainty > 1 kcal/mol) |
| Energy samples | every 5 ps → 400–1000 samples/window |
| Platform | Sherlock CUDA (`py-openmm` module); local M3 Max for prep/QC only |

**No independent fixed-λ MBAR re-eval across Perses hybrids.** Use either:

- **MBAR** from fixed-λ windows with **clash filter**: discard samples where
  `(max_k u_k - min_k u_k) > 50 kJ/mol`, or  
- **BAR** between adjacent windows if MBAR unstable.

**Do not** trust λ=0.8–1.0 windows if cross-state spread exceeds threshold — for scaling FEP,
monitor spread from λ=0 onward.

### Forward / reverse

For pilot mutations, run **reverse leg** (mutant→WT annihilation reversed) or duplicate
half-window set to check internal consistency (Smith et al. criterion: forward ≈ −reverse
within ~1 kcal/mol).

---

## 6. Implementation plan

### 6.1 Directory layout (this package)

```
scripts/fep_truncated/
  PLAN.md                 ← this file
  README.md
  config.py               ← TruncatedFEPConfig dataclass
  pocket.py               ← 15 Å extraction, residue maps
  restraints.py           ← backbone fixed, shell restraints
  solvent.py              ← water cap builder
  alchemical.py           ← atom selection + NB scaling
  prepare.py              ← build holo/apo OpenMM systems → XML
  worker.py               ← one λ window MD
  analyze.py              ← MBAR/BAR + ΔΔG_bind
  panel.py                ← manuscript legs manifest
  convergence.py          ← cross-state spread, clash fraction
  tests/
```

Results:

```
results/analysis/fep_truncated/
  legs/{leg_id}/
    config.json
    holo/  apo/
      system.xml, truncated.pdb, schedule.json
      windows/state_XX_energies.csv
  targets/{genotype}/summary.json
```

### 6.2 Reuse from `fep_jorgensen`

| Module | Reuse |
| --- | --- |
| `mutations.py` | Import `MANUSCRIPT_PLANS`, `MutationLeg`, `Mutation.parse` |
| `equilibrate.py` | Optional pre-MD on full structure before truncation |
| `analyze.py` cycle logic | Adapt ΔΔG_bind aggregation |
| `rsync_fep_jorgensen.sh` | Template for `rsync_fep_truncated.sh` |
| Sherlock submit scripts | Clone pattern with smaller resource requests |

Do **not** import Perses / `perses_hybrid.py` in v1.

### 6.3 Implementation phases

| Phase | Deliverable | Exit criterion |
| --- | --- | --- |
| **0** | This plan + branch | Reviewed |
| **1** | `pocket.py` + `prepare.py` for V106A holo | Truncated PDB ~120–180 residues; `residue_map` shows auth 106 → PDB 103 |
| **2** | Restraints + water cap + holo/apo systems | Stable 50 ps MD, no clashes |
| **3** | `worker.py` 11 windows × holo V106A | Cross-state spread < 50 kJ/mol all windows |
| **4** | Apo leg + `analyze.py` | \(\Delta\Delta G_\text{bind}\) with uncertainty |
| **5** | Pilot panel: V106A, Y188L, K103N | Rank-order vs experimental fold-change |
| **6** | Batch panel + manuscript table export | 19 targets or documented failures |

---

## 7. Validation criteria

### Per-leg QC (automated)

- [ ] All λ windows complete (≥400 samples)
- [ ] Max cross-state energy spread < **50 kJ/mol** (stricter than Perses pilot’s failed threshold)
- [ ] Clash fraction (spread > 500 kJ/mol) = **0%**
- [ ] Forward/reverse \(\Delta G_\text{mut}\) within **1 kcal/mol** (pilot set)

### Scientific (pilot trio)

| Mutation | Experimental DOR fold-change (manuscript) | Pass if |
| --- | --- | --- |
| V106A | (from table) | Sign of \(\Delta\Delta G_\text{bind}\) matches resistance direction |
| Y188L | (from table) | Same |
| K103N | (from table) | Same |

Absolute kcal/mol agreement within **±1–2 kcal/mol** is aspirational for v1; **ranking and sign**
are the hard gates before full panel.

---

## 8. Compute estimate (order of magnitude)

Per single-mutation leg (holo + apo):

- 11 windows × 2 phases = **22 windows**
- ~2–5 ns × ~500 atoms in cap ≈ **minutes–tens of minutes per window** on GPU
- **~22 × 15 min ≈ 5–6 GPU-hours** per leg (rough)

Full panel: 19 unique legs → **~100–120 GPU-hours** total (vs thousands for full-protein Perses).

---

## 9. Non-goals (v1)

- Exact MCPRO / IMPACT reproduction
- Perses hybrid topologies or DELETE/INSERT atoms
- Full-protein PME alchemical MD
- Alchemical replica exchange (AREX) unless pilot fails QC
- Charge-change special protocols (defer K103N deep dive until charge-neutral pilots pass)
- Free-energy differences **between inhibitors** (Jorgensen Table 2 multi-drug) — only DOR for now

---

## 10. Open decisions (resolve in Phase 1 pilot)

1. **Water cap vs implicit apo** — start explicit cap for both phases.  
2. **Mutant endpoint coordinates** — minimize mutant side chain after annihilation vs use existing mutant MD PDB as λ=1 reference.  
3. **Restraint spring constants** — scan 10 vs 50 kcal/mol/Å² on shell.  
4. **Production length** — 2 ns default; extend to 5 ns if MBAR error > 0.5 kcal/mol.  
5. **Whether to quench/reheat** before FEP (Jorgensen 2000 MD protocol) — optional 30 min wall-time test on V106A holo.

---

## 11. References

1. Rizzo RC, Wang DP, Tirado-Rives J, Jorgensen WL. JACS 2000. [doi:10.1021/ja003113r](https://doi.org/10.1021/ja003113r)  
2. Smith RH Jr et al. Bioorg Med Chem Lett 2007. [doi:10.1016/j.bmcl.2007.12.033](https://doi.org/10.1016/j.bmcl.2007.12.033)  
3. Repo: `docs/Jorgensen-FEP-protocol.md`, `scripts/fep_jorgensen/README.md` (Perses pilot findings)
