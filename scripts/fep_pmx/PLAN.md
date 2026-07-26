# Plan: DOR resistance ΔΔG_bind via pmx + GROMACS NEQ

Status: **design** — implementation on branch `fep-pmx` off `main`.  
Supersedes the truncated-sphere / Perses fixed-λ plan.

**Method:** pmx hybrid topology + GROMACS **non-equilibrium switching** (Crooks/BAR) on
**full solvated** holo and apo systems. Open-source analogue of Schrödinger FEP+ for protein
mutation binding affinity.

Background: [`APPROACHES.md`](APPROACHES.md) · Perses post-mortem:
[`scripts/fep_jorgensen/README.md`](../fep_jorgensen/README.md)

---

## 1. Scientific goal

For every manuscript genotype \(G\):

\[
\Delta\Delta G_\text{bind}(G) = \Delta G_\text{mut}^\text{holo} - \Delta G_\text{mut}^\text{apo}
\]

- **Holo:** alchemical WT→mut in **DOR-bound** full RT (rhombic dodecahedron, PME).  
- **Apo:** same mutation in **ligand-free** full RT (identical box/restraint policy).  
- **Positive ΔΔG_bind:** resistance direction (weaker binding).

**Experimental reference:**
`results/analysis/dor_susceptibility_bar_chart/tables/dor_susceptibility_values.csv`
(`dor_fold_reduction`). Primary gate: **Spearman rank** vs experiment (P1 ρ ≥ 0.7); absolute
ΔΔG is secondary given apo pocket uncertainty (§6).

**Compound genotypes:** additive leg sum from `MANUSCRIPT_PLANS` in
`scripts/fep_jorgensen/mutations.py` (19 unique legs → 19 targets).

---

## 2. Why this method (expert consensus)

### 2.1 Direct precedent

| Paper | What it proves for us |
| --- | --- |
| **Aldeghi, Gapsys & de Groot**, ACS Cent. Sci. 2018 ([doi:10.1021/acscentsci.8b00717](https://doi.org/10.1021/acscentsci.8b00717)) | **134 mutations, 17 proteins, 27 ligands** — ΔΔG_bind upon *protein* mutation; pmx hybrids + NEQ; RMSE **1.2 kcal/mol** (0.8 on reproducible subset) |
| **Aldeghi et al.**, ACS Cent. Sci. 2019 ([doi:10.1021/acscentsci.9b00590](https://doi.org/10.1021/acscentsci.9b00590)) | Same NEQ/pmx applied to **kinase inhibitor resistance** (Abl, 144 ΔΔG) — closest published analogue to NNRTI resistance ranking |
| **Gapsys et al.**, JCC 2015 ([pmx](https://doi.org/10.1002/jcc.23804)) | Automated hybrid topology for all amino-acid pairs |
| **Aldeghi**, Methods Mol Biol 2019 | NEQ workflow tutorial (Crooks/BAR) |
| **King et al.**, JMB 2019 | Co-alchemical ion/water for charge-changing mutations (RMSE ~1.2 kcal/mol) |
| **Gapsys et al.**, JCIM 2020 | pmx double-system/single-box for charge mutations |

Jorgensen MCPRO (HIV RT V106A/Y181C) remains useful as **physics intuition** for the holo/apo
cycle, but is not the implementation path.

### 2.2 Why NEQ fixes the Perses endpoint failure

Equilibrium fixed-λ FEP must **sample and converge** unphysical intermediates. Perses DELETE
path at λ≥0.8: ghost atoms clash when re-evaluated at low λ → MBAR blowups (repo pilot).

NEQ restructures the problem:

1. **Equilibrium only at λ=0 and λ=1** — physical WT and mutant ensembles.  
2. **Short driven switches** (~100 ps) from endpoint snapshots; work recorded; pathological
   transitions = **outliers**, not poisoned equilibrated states.  
3. **Forward + reverse** work distributions → BAR via **Crooks fluctuation theorem**.

This is the same failure-class fix FEP+ uses in practice; NEQ is the open pmx implementation.

### 2.3 Why full system, not truncated sphere

Truncation was a **cost workaround** for expensive equilibrium λ-window FEP. NEQ is
embarrassingly parallel (hundreds of independent ~100 ps array jobs) and makes full PBC
affordable (~1,700 GPU-h for 19 legs). Truncation adds boundary artifacts (Genheden 2012,
Huang 2016) with no compensating benefit under NEQ budgeting.

**Dropped:** sphere radius, shell restraints, nonperiodic cap, P0-geom bake-off.

---

## 3. Manuscript panel

Unchanged from prior plan — full 19 targets, 19 unique legs. See
`scripts/fep_jorgensen/mutations.py`.

### 3.1 Pilot (P0)

Run **first**, before batching:

| Leg | Why |
| --- | --- |
| `wt_to_V106A` | Small neutral; should just work |
| `wt_to_Y188L` | Largest perturbation; strongest signal (149×); 500 ps switches |

If both pass Crooks overlap + correct sign → method proven; remaining 17 legs are bookkeeping.

### 3.2 Rollout

| Stage | Scope |
| --- | --- |
| **P0** | V106A + Y188L (holo + apo each) |
| **P1** | All 11 WT single-residue legs → Spearman vs CSV |
| **P2** | 8 compound legs → additivity check |
| **P3** | Full 19-target table for manuscript |

### 3.3 Residue numbering

Manuscript **auth** IDs (V106) ≠ PDB **resSeq** (VAL103). Every leg emits `residue_map.json`.
Reuse auth/PDB logic from `scripts/fep_jorgensen/prepare_backend.json` patterns.

---

## 4. Concrete protocol

Per **leg × phase** (holo or apo), per replicate:

### 4.1 Force field & box

| Setting | Value |
| --- | --- |
| Protein | AMBER **ff14SB** (match existing MD) |
| Ligand (holo) | **OpenFF 2.0.0** via OpenMM → export for GROMACS (ACPYPE/interchange or pre-parametrized `.itp`) |
| Water | TIP3P |
| Box | **Rhombic dodecahedron**, ~1.0 nm buffer (match existing solvation) |
| Ions | 0.15 M NaCl; **no ions in protein interior** (Aldeghi 2018 reproducibility fix) |

### 4.2 Hybrid topology (pmx)

1. Start from equilibrated **physical** PDB per endpoint (WT or mutant; holo or apo).  
2. `pmx mutate` → hybrid structure; `gentop` → dual-state `.top`.  
3. **pmx `develop` branch** for P225H (proline ring opening) — test early, not at P2.  
4. His protonation: **neutral HID/HIE** at pH 7 unless evidence otherwise (P225H is formally charge-neutral).

### 4.3 Endpoint equilibration

Existing plain-MD trajectories **cannot** be used directly (no dummy atoms in hybrid top).

**What they provide:** decorrelated starting conformations per genotype.

| Step | Setting |
| --- | --- |
| Replicates | **3** |
| Seed frames | Decorrelated snapshots from `results/md_runs/{genotype}/rep_*/` (and apo analogues) |
| Hybrid relaxation | **5 ns** production per endpoint per replicate (λ pinned at 0 or 1) |
| Total equil per leg-phase | 3 reps × 2 endpoints × 5 ns = **30 ns** |

Aldeghi 2018 optimal split: ~equal time in equilibrium and switching. Our budget targets ~60 ns
equil + ~60 ns switching per leg-phase.

Each endpoint runs a **three-step** relaxation at fixed λ before the 5 ns production
trajectory that snapshots are drawn from (`run_neq_task.py::_run_equil`):

1. **Per-λ minimization** (`em_fep.mdp`, free-energy + gapsys soft-core at the endpoint λ).
   The global `em` stage minimizes only the A-state, so B-state atoms that are dummies
   in A — i.e. the *grown* sidechain of a growth mutation such as **G190E at λ=1** — are
   otherwise never relaxed. This pass relieves those forces before any dynamics.
2. **C-rescale warmup** (`npt_warmup.mdp`, 0.5 ns). Starting Parrinello–Rahman directly
   from a minimized structure with generated velocities can blow up a ~200k-atom box;
   the stochastic-cell-rescaling barostat relaxes the box first and is a valid NPT
   ensemble ([Bernetti & Bussi, JCP 2020](https://doi.org/10.1063/5.0020514)).
3. **Parrinello–Rahman production** (`npt_eq.mdp`, 5 ns), `continuation = yes` from the
   warmup checkpoint. Sampling trajectory.

Temperature coupling uses **separate `Protein` / `non-Protein` baths** (not a single
`System` group) throughout warmup, production, and switching to avoid the
hot-solvent/cold-solute artifact.

### 4.4 Non-equilibrium switching

| Setting | Default | Exceptions |
| --- | --- | --- |
| Snapshots per endpoint per replicate | **100** | — |
| Switch length | **100 ps** | **500 ps** for Y188L, G190E |
| Direction | Forward (λ 0→1) + reverse (λ 1→0) from each snapshot | — |
| Replicates | 3 | — |
| Soft-core | `sc-function = gapsys` (GROMACS ≥2022) | — |
| Gapsys defaults | `sc-gapsys-scale-linpoint-q=0.3`, `sc-gapsys-scale-linpoint-lj=0.85` | — |

**Total switching per leg-phase:** 3 × 100 × 100 ps × 2 directions ≈ **60 ns** (600 switches).

### 4.5 Analysis

- Pool work values; **BAR** via Crooks (pmx `analyse -m bar` or equivalent).  
- Per-replicate ΔG with error bars; combine replicates.  
- **ΔΔG_bind** = (ΔG_mut^holo − ΔG_mut^apo) for each leg.

### 4.6 QC (replaces Perses cross-state spread)

| Metric | Pass |
| --- | --- |
| Crooks forward/reverse work distribution overlap | Visual + automated overlap test |
| Work outlier fraction | Flag if >5% extreme outliers; investigate, optionally discard |
| BAR vs Jarzynski | Agreement within 1 kcal/mol |
| Sign check (P0) | V106A, Y188L positive ΔΔG_bind vs CSV |
| Replicate spread | Per-leg SEM reported |

---

## 5. Compute budget (Sherlock)

Assumption: ~40 ns/GPU-day for ~200k-atom system.

| | Per leg-phase | Per leg (holo+apo) | Panel (19 legs) |
| --- | --- | --- | --- |
| Equilibrium | ~30 ns | 60 ns | ~1.1 µs |
| Switching | ~30 ns | 60 ns | ~1.1 µs |
| **GPU-hours** | ~18 | ~36 | **~1,700** |

| Milestone | GPU-h | Wall clock @ 20 GPUs |
| --- | --- | --- |
| **P0** (2 legs) | ~110 | Few days + queue |
| **Full panel** | ~1,700 | ~2 weeks + queue |

NEQ array shape: **600 switches × 3 reps × 2 phases × 19 legs** ≈ 68k independent ~100 ps
jobs — ideal for `sbatch --array`.

---

## 6. Panel-specific traps

### 6.1 P225H (PRO→HIS)

Classic pmx exclusion (ring opening = bond breaking). **BioExcel `develop` branch** now
supports this — install develop and **test in P0 extension** before P2 compound legs
(`*_P225H`).

### 6.2 K103N and G190E (net charge change)

Only true charge-changing legs in the panel. Pick **one** route, use consistently:

| Option | Ref | Notes |
| --- | --- | --- |
| **Co-alchemical ion/water** | King JMB 2019; FEP+ implementation | RMSE ~1.2 kcal/mol on charge muts |
| **Double-system / single-box** | Gapsys JCIM 2020; pmx-native | 30 Å separation |

Do **not** mix schemes within the panel. Run charge legs after neutral P0 passes.

### 6.3 Apo leg — primary scientific risk

The NNRTI pocket is **created by inhibitor binding**. Unliganded RT collapses on timescales
longer than 5–10 ns equilibration.

**Mitigations (both required in manuscript):**

1. **Rank over absolute:** P1 gate is Spearman ρ ≥ 0.7 — apo error may largely cancel as a
   systematic offset across the panel if treatment is identical.  
2. **Identical weak pocket restraints** in every apo leg (same force constant, same anchor
   atoms, same duration) — document explicitly if reporting absolute ΔΔG.  
3. **Never** vary apo equilibration length between legs.

Holo leg is the well-defined physics; apo is the honest uncertainty.

---

## 7. Software stack

### 7.1 Sherlock (verified 2026-03 on gpu partition)

| Role | Load / tool | Notes |
| --- | --- | --- |
| **GPU MD** | `source scripts/sherlock/load_gromacs_module.sh` | Sets `GMX_MDRUN=gmx_cuda` |
| **Module** | `chemistry gromacs/2023.1` | Auto-pulls `cuda/11.2.0` — do not force cuda/12.6.1 |
| **Prep** | `gmx` (grompp, solvate, editconf) | CPU binary in same module |
| **Production** | **`gmx_cuda mdrun`** | `GPU support: CUDA`, cuFFT; tested on Tesla V100 |
| **Skip** | `gromacs/2025.1` | CPU-only `gmx`; no `gmx_cuda` in that install |

**`.mdp` soft-core** (not in `-h`; validated at grompp):

```
sc-function               = gapsys
sc-gapsys-scale-linpoint-q  = 0.3
sc-gapsys-scale-linpoint-lj = 0.85
```

**Do not** `module load py-openmm` in the same shell as GROMACS — CUDA stacks conflict
(openmm → cuda/12.4; gromacs/2023.1 → cuda/11.2).

SLURM GPU jobs: `#SBATCH -p gpu -G 1`; use `gmx_cuda -nb gpu` (or `-gpu_id 0`).

### 7.2 pmx + parametrization (not Sherlock modules)

| Component | Install | Where |
| --- | --- | --- |
| **pmx** | `bash scripts/fep_pmx/setup_pmx_env.sh` — [deGrootLab/pmx](https://github.com/deGrootLab/pmx) **`develop`**; **not** `pip install pmx` | **Mac** (prep/analysis); optional on Sherlock login |
| **pmx analyze** | same venv | Mac or login |
| **OpenMM / OpenFF DOR `.itp`** | existing Mac env | **Mac** — not mixed into GROMACS GPU jobs |
| **Schrodinger FEP+** (optional) | `module load chemistry schrodinger/2024-1` | External P0 benchmark only |

**Mac M3 Max:** pmx topology generation, `.mdp` validation, all analysis. Not production MD.

**Rejected for v1:**

| Approach | Reason |
| --- | --- |
| Truncated sphere | Unnecessary under NEQ budget; adds artifacts |
| Perses fixed-λ full protein | Endpoint DELETE failure (pilot) |
| Perses AREX | Complexity without benefit once NEQ chosen |
| OpenMM alchemical scaling | No mature protein-mutation NEQ in OpenMM |
| Amber TI (`pmemd.cuda`) | Licensing friction |
| OpenFE/feflow | Conda-centric; immature protein mutation support |

---

## 8. Implementation layout

```
scripts/fep_pmx/
  PLAN.md, README.md, APPROACHES.md
  setup_pmx_env.sh       # Mac: conda env pmx + deGrootLab/pmx develop
  config.py              # ff, box, switch lengths, sc-function params
  panel.py               # wrap mutations.py legs
  prepare_hybrid.py      # pmx mutate + gentop from md_runs PDBs
  mdp/                   # equil, prod_lambda0/1, nonequil forward/reverse templates
  submit_equil.sh        # 3 reps × 2 endpoints
  submit_switch.sh       # SLURM array: 100 switches × reps × directions
  analyze.py             # BAR/Crooks, per-leg ΔG, ΔΔG_bind, QC plots
  charge_protocol.py     # Tier C: co-ion or double-box for K103N, G190E
  tests/

scripts/sherlock/
  load_gromacs_module.sh # chemistry gromacs/2023.1 + gmx_cuda (see PLAN §7.1)
  setup_pmx_env.sh       # optional Sherlock login install (python/3.9)
  connect.sh, remote.sh

results/analysis/fep_pmx/
  legs/{leg_id}/{holo,apo}/rep_{01,02,03}/
  targets/{genotype}/summary.json
  panel_qc.csv
```

Reuse: `scripts/fep_jorgensen/mutations.py`, `results/md_runs/` paths, Sherlock rsync patterns.

---

## 9. Validation gates

### P0 (V106A + Y188L)

- [ ] Crooks overlap clean (forward/reverse)  
- [ ] BAR/Jarzynski agree  
- [ ] Sign(ΔΔG_bind) matches CSV (both positive resistance)  
- [ ] Replicate error bars < 1 kcal/mol (target, not hard fail)

### P1 (11 singles)

- [ ] Spearman ρ ≥ 0.7 vs `dor_fold_reduction`  
- [ ] V106I, K103N, Y181C near bottom of rank (negative controls)

### P2 (8 compound legs)

- [ ] Additive sum within 1 kcal/mol for ≥6/8 targets

### P3

- [ ] Manuscript-ready table: ΔΔG_bind ± SEM per genotype, experimental fold, rank

---

## 10. References (primary)

1. Aldeghi M, Gapsys V, de Groot BL. ACS Cent. Sci. **2018**, 4, 1708–1718.  
2. Aldeghi M, Gapsys V, de Groot BL. ACS Cent. Sci. **2019**, 5, 515–522 (kinase resistance).  
3. Gapsys V et al. J. Comput. Chem. **2015**, 36, 348–354 (pmx).  
4. Aldeghi M et al. Methods Mol. Biol. **2019** (NEQ tutorial).  
5. King NM et al. J. Mol. Biol. **2019** (co-alchemical water).  
6. Gapsys V, de Groot BL. J. Chem. Inf. Model. **2020** (double-system).  
7. Smith MB et al. Bioorg. Med. Chem. Lett. **2007** (Jorgensen HIV RT — cycle precedent).  
8. Repo Perses pilot: `scripts/fep_jorgensen/README.md`.
