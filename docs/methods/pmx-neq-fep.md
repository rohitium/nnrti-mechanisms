# pmx + GROMACS NEQ FEP — current plan

**Status (Jul 2026):** P0 pilot on Sherlock (`wt_to_V106A`, `wt_to_Y188L`; holo + apo; 3 replicates).  
Interactive GPU smoke test before batch `sbatch`.  
Implementation: `scripts/fep_pmx/` · Design detail: `scripts/fep_pmx/PLAN.md`.

---

## 1. Quantity of interest

For each resistance mutation leg (WT → mutant), estimate the change in **binding free energy**:

\[
\Delta\Delta G_\text{bind} = \Delta G_\text{mut}^\text{holo} - \Delta G_\text{mut}^\text{apo}
\]

- **Holo:** alchemical mutation in DOR-bound HIV-1 RT (full solvated periodic system).
- **Apo:** same mutation in ligand-free RT (same box/ion/restraint policy).
- Positive ΔΔG_bind → weaker DOR binding (resistance direction).

Thermodynamic cycle for relative resistance: mutate the protein in bound and unbound states ([Rizzo et al., JACS 2000](https://doi.org/10.1021/ja003113r); HIV-RT mutation energetics: [Smith et al., BMCL 2008](https://doi.org/10.1016/j.bmcl.2007.12.033)).

**P0 experimental targets:** V106A ≈ RT·ln(9.6) = **1.35 kcal/mol**; Y188L ≈ RT·ln(149) = **2.98 kcal/mol** (`dor_fold_reduction` CSV).

**Validation:** Spearman rank vs experiment (P1 ρ ≥ 0.7); absolute kcal/mol is secondary given apo-pocket uncertainty (§8).

---

## 2. Method: hybrid topology + non-equilibrium switching

### 2.1 Hybrid Hamiltonian (pmx)

Side-chain mutation uses **pmx hybrid residues** — shared backbone/common atoms, dummy atoms only where WT and mutant side chains differ. This is **not** dual topology (both full side chains present simultaneously). Automated parameterization: [Gapsys et al., J. Comput. Chem. 2015](https://doi.org/10.1002/jcc.23804).

Force field: AMBER ff14SB protein (pmx `amber14sbmut`), OpenFF 2.0.0 ligand (holo), TIP3P water, 0.15 M NaCl, rhombic dodecahedron ~1 nm padding — matched to existing OpenMM MD prep.

**OpenMM → GROMACS validation (required once per build):** single-point energy comparison on one frame, per component (protein, ligand, solvent), OpenMM vs `gmx mdrun -rerun`, agreeing to **< 0.1 kJ/mol**. Silent parameter loss in Interchange export is common and invisible to downstream QC.

### 2.2 Why NEQ instead of equilibrium λ-windows

**Equilibrium FEP/TI** requires adequate sampling at every λ where the Hamiltonian is unphysical. Our Perses fixed-λ pilot on full protein failed MBAR convergence when ghost-atom intermediates clashed at high λ.

**Non-equilibrium switching (NEQ)** avoids long equilibration in intermediates ([Jarzynski, Phys. Rev. Lett. 1997](https://doi.org/10.1103/PhysRevLett.78.2690); [Crooks, Phys. Rev. E 1999](https://doi.org/10.1103/PhysRevE.60.2721)):

1. **Equilibrate only at λ = 0 and λ = 1** (physical WT and mutant endpoints).
2. **Drive** λ linearly over a short time (100–500 ps); record the **work** \(W\) on the hybrid coordinates.
3. Draw **N decorrelated snapshots from each endpoint ensemble** → **N forward** switches (λ 0→1, starting from λ=0 frames) **+ N reverse** switches (λ 1→0, starting from λ=1 frames). Do not run reverse switches from λ=0 snapshots.
4. Estimate ΔG from work distributions via **Bennett acceptance ratio (BAR)** ([Bennett, J. Comput. Phys. 1976](https://doi.org/10.1016/0021-9991(76)90078-4); [Shirts et al., JCP 2003](https://doi.org/10.1063/1.1618247)).

Benchmark for protein–ligand ΔΔG_bind upon protein mutation: pmx + NEQ, 134 mutations, RMSE ~1.2 kcal/mol ([Aldeghi, Gapsys & de Groot, ACS Cent. Sci. 2018](https://doi.org/10.1021/acscentsci.8b00717)); kinase resistance analogue ([Aldeghi et al., ACS Cent. Sci. 2019](https://doi.org/10.1021/acscentsci.9b00590)). Tutorial: [Aldeghi, Methods Mol. Biol. 2018](https://doi.org/10.1007/978-1-4939-8736-8_2).

Soft-core decoupling: Gapsys linearized point charges ([Gapsys et al., JCTC 2012](https://doi.org/10.1021/ct3000814); GROMACS `sc-function = gapsys`).

### 2.3 Why full solvated RT (not truncated sphere)

Jorgensen MCPRO used ~15 Å binding-site spheres and water caps ([Rizzo et al., 2000](https://doi.org/10.1021/ja003113r)). Truncation reduces cost but introduces boundary artifacts ([Genheden & Ryde, JCTC 2012](https://doi.org/10.1021/ct200853g)). NEQ switches are short and **embarrassingly parallel**, so full PBC systems are affordable.

---

## 3. Computational protocol (per leg × phase × replicate)

| Stage | Engine | What | Setting |
|-------|--------|------|---------|
| **Build** | pmx + GROMACS | Hybrid top, solvate, ions; OpenMM energy check | From OpenMM MD start structures |
| **EM** | GROMACS CPU | Minimize hybrid system (A-state) | — |
| **Equil** | GROMACS GPU | Per-λ min → C-rescale warmup → P-R production, λ fixed at 0 or 1 | **5 ns** production / endpoint / rep (2 endpoints × 3 reps) |
| **Extract** | GROMACS CPU | Decorrelated frames + velocities for switches | **100 snapshots** / endpoint / rep (skip first 100 ps); `.trr` required |
| **Switch** | GROMACS GPU | NEQ MD, linear λ ramp | 100 ps default; **500 ps** for Y188L, G190E |
| **Analysis** | pmx BAR | ΔG per leg-phase-rep; then ΔΔG_bind | Crooks overlap QC |

### 3.1 Snapshot count — non-negotiable for P0

P0 gate requires **Crooks overlap** and **correct resistance sign**. With **5 work values per direction** there is no work distribution to histogram, and BAR error on n=5 is ~1.5–3 kcal/mol — larger than the V106A signal (1.35 kcal/mol).

**Minimum ~50 snapshots/endpoint; P0 uses 100.** Switching is the cheap half of NEQ (~40 ns/GPU-day, ~200k atoms):

| | 5 snapshots | 100 snapshots |
|---|---|---|
| Equilibration (2 × 5 ns) | 6.0 GPU-h | 6.0 GPU-h |
| Switching (100 ps each) | 0.6 GPU-h | 12 GPU-h |

P0 total: ~40 → ~145 GPU-h (~one extra queue day) for an interpretable pilot vs an uninterpretable one. Prior cut to 5 saved **job count**, not compute — solve job count by **bundling snapshots per SLURM array task** (100 for V106A, 50 for Y188L; see §7.2), not by discarding statistics.

If P0 must be cheaper: cut **3 → 2 replicates** or **5 → 3 ns equilibration**. **Never cut snapshots below ~50.**

### 3.2 Error bars

- **Point estimate:** pool work values within each replicate for BAR ΔG.
- **Uncertainty:** report spread **across the 3 replicates** (e.g. SEM of replicate ΔG values). Pooling all work into one BAR error **underestimates** uncertainty because replicates share correlated endpoint ensembles.

Plain MD trajectories seed **conformations** only — hybrid topologies require fresh endpoint equilibration ([Aldeghi 2018](https://doi.org/10.1021/acscentsci.8b00717)).

### 3.3 Endpoint relaxation (`_run_equil`)

Each endpoint runs three fixed-λ steps before the 5 ns production trajectory:

1. **Per-λ minimization** (`em_fep.mdp`) with free-energy + gapsys soft-core at the endpoint λ. The global `em` stage minimizes only the A-state; B-state atoms that are dummies in A (the grown sidechain of **G190E at λ=1**, etc.) are unrelaxed until this pass.
2. **C-rescale warmup** (`npt_warmup.mdp`, 0.5 ns) to relax the box and generate velocities — starting Parrinello–Rahman from a bare minimized structure risks box blow-up on a ~200k-atom system ([Bernetti & Bussi 2020](https://doi.org/10.1063/5.0020514)).
3. **Parrinello–Rahman production** (`npt_eq.mdp`, 5 ns, `continuation = yes`) — the sampling trajectory.

Temperature coupling uses separate **`Protein` / `non-Protein`** baths (not `System`) in all dynamics stages. Before batch extract, confirm `equil.trr` time range with `gmx check` (production clock may or may not include the 500 ps warmup offset).

### 3.4 P0 pilot results & caveats (first run: 100 ps V106A / 500 ps Y188L, 3 reps)

| genotype | ΔΔG_bind (kcal/mol) | exp. fold | sign | overlap (r1,r2,r3) |
|---|---|---|---|---|
| V106A | **+1.69 ± 0.70** | 9.6 | ✅ + | 0.42, 0.06, 0.08 |
| Y188L | **+4.52 ± 0.49** | 149 | ✅ + | holo 0.49/0.24/0.01 |

**What passed:** both ΔΔG_bind positive (resistance direction); ranking correct (Y188L ≫ V106A, matching 149 ≫ 9.6); replicate SEM < 1 kcal/mol; BAR `Conv ≈ 0` and BAR/CGI/Jarzynski agree within ~1 kcal/mol.

**The caveat — marginal Crooks overlap.** Per-rep forward/reverse overlap is 0.01–0.53 (most < 0.3). Forward works are tight and reproducible; **reverse (λ=1 / mutant) works are broad and their means scatter across replicates** (e.g. V106A reverse mean 15–19 kcal/mol vs a stable forward ~21). That asymmetry drags overlap down. It is not disqualifying — BAR converges and the ΔΔG signs/ranking/SEM hold — but per-phase absolute ΔG is less certain than the SEM alone implies, worst for **apo** legs (pocket instability, §6.3).

**Levers (which actually help overlap):**
- **Switch length** — longer switches → less dissipation → tighter fwd/rev gap. Env-overridable via `NEQ_LONG_SWITCH_PS` / `NEQ_EXTRA_LONG_SWITCH_LEGS`.
- **Endpoint equilibration** — targets the noisy λ=1 ensemble. `NEQ_EQUIL_NS` (env-overridable; applies to both endpoints).
- **Snapshots do NOT help overlap** — they only shrink ΔG statistical error on the *same* distributions. The P0 SEM is already < 1 kcal/mol, so snapshots are not the bottleneck.

### 3.5 Switch-length test — ΔΔG is switch-length-invariant

V106A was re-run at **500 ps** (5× the 100 ps default), 3 reps, to test whether longer switches fix the marginal overlap and/or move ΔΔG:

| | 100 ps | 500 ps |
|---|---|---|
| ΔΔG_bind (kcal/mol) | +1.69 ± 0.70 | +1.76 ± 0.51 |
| overlap (holo r1/r2/r3) | 0.42 / 0.06 / 0.08 | 0.57 / 0.12 / 0.08 |

**ΔΔG is statistically identical**; overlap barely improved. Interpretation: the marginal overlap reflects **intrinsic dissipation of the NNRTI pocket mutation, not under-sampling** — so it does not bias ΔΔG (confirmed by 100 ps ≡ 500 ps, plus BAR/CGI/Jarzynski agreement and replicate agreement). **The panel therefore runs at 100 ps** — longer switches cost 5× for no gain in the ranking quantity (P1 ≈ 450 vs ≈ 2,000 GPU-h). The Crooks-overlap gate is treated as **conservative** here: report it, but marginal overlap with switch-length-invariant, estimator-consistent ΔΔG is acceptable for the ranking POC.

**Overlap-metric note:** compare `W_f` vs `W_r` **directly** — pmx already stores the reverse work in the forward frame (`integ_rev` crosses `integ_fwd` at ΔG). `qc_neq.py` initially negated `W_r`, which produced false near-zero overlap for large-ΔG legs (V106A) and false-high overlap for ΔG≈0 legs (Y188L holo); fixed.

---

## 4. Panel rollout

| Stage | Legs | Gate |
|-------|------|------|
| **P0** | `wt_to_V106A`, `wt_to_Y188L` | Crooks overlap; correct resistance sign vs CSV |
| **P0.5** | `wt_to_K103N` | Charge-changing protocol pilot (see §5) before P1 ranking |
| **P1** | 11 WT single-residue legs | Spearman ρ ≥ 0.7 vs experiment |
| **P2** | 8 compound targets (background-dependent legs) | Epistatic additivity (see below) |
| **P3** | 19 targets | Manuscript table |

19 unique mutation legs from `scripts/fep_jorgensen/mutations.py`. Compound genotypes are **paths**, not sums of singles:

| Target | Leg 1 | Leg 2 (mutant background) |
|--------|-------|----------------------------|
| `L100I+K103N` | `wt_to_K103N` | `K103N_to_L100I_K103N` |
| `K103N+P225H` | `wt_to_K103N` | `K103N_to_K103N_P225H` |
| … | … | … |

**P2 additivity gate:** for each compound, compare ΔΔG_bind of **leg 2** (e.g. `K103N_to_L100I_K103N`) against the **single-residue leg on WT background** (e.g. `wt_to_L100I`). Deviation = epistasis. Do **not** compare a summed path to itself.

**P0 extension (parallel with P0 batch):** test pmx **`develop`** branch on `P225H` (proline ring opening) — mainline excludes PRO; three P2 legs contain P225H.

---

## 5. Residue-specific protocol requirements

### 5.1 Net charge change (K103N, G190E)

Plain PME gives nonsense ΔG when net charge changes. Pick **one** scheme, apply identically to holo and apo, **before P1**:

| Option | Reference | Notes |
|--------|-----------|-------|
| Co-alchemical ion/water | Clark et al., J. Mol. Biol. 2019 ([doi:10.1016/j.jmb.2019.02.003](https://doi.org/10.1016/j.jmb.2019.02.003)) | FEP+ implementation; RMSE ~1.2 kcal/mol on charge muts |
| Double-system / single-box | Gapsys et al., JCC 2015; Bhatt et al., JCTC 2021 | pmx-native; ~30 Å separation ([doi:10.1021/acs.jctc.0c01045](https://doi.org/10.1021/acs.jctc.0c01045)) |

**G190E** is the only leg with a true net charge change (+1). **K103N** is formally neutral but sits in a charged entrance — run in **P0.5** to validate the chosen protocol before it enters the P1 ranking gate.

### 5.2 Proline (P225H)

pmx mainline excludes proline; **`develop` branch** handles ring opening. Validate on P225H during P0, not at P2.

### 5.3 Histidine protonation

**P225H:** explicit decision — **HID or HIE at pH 7** (charge-neutral); do not leave to `pdb2gmx` defaults. Document in leg `config.json`.

---

## 6. Apo leg QC

The NNRTI pocket collapses without inhibitor. **Before batching apo legs:** plot pocket volume vs time during apo equilibration (same metric as `compute_nnbp_pocket_volume.py`). If volume has not plateaued by 5 ns, work distributions depend on equilibration length and the leg is not reproducible — extend equilibration or apply identical weak pocket restraints across **all** apo legs (document force constant and anchor atoms).

Holo leg is well-defined physics; apo is the primary scientific uncertainty. P1 gate is rank (Spearman), not absolute ΔΔG.

---

## 7. Operational workflow (Sherlock)

**Rule:** validate on **interactive GPU (`salloc` on `dev`)** before **batch (`sbatch` on `gpu`)**.

```
salloc -p dev → smoke_neq_em → smoke_neq_task (equil → extract → switch)
             → audit_neq_panel.py
             → submit_p0_neq_pipeline.sh
```

### 7.1 GPU partition QOS — the hard wall

```
Partition gpu:  QoS=gpu
QOS gpu:  MaxSubmitPU=100   MaxTRESPU=cpu=128,gres/gpu=16   MaxWall=48:00:00
```

- **100 submitted GPU jobs** (queued + running), **16 concurrent GPUs**, **128 CPUs**, **48 h wall**.
- Array elements count individually against `MaxSubmitPU` — `--array=0-227` is rejected regardless of `%` throttle.
- `TRESBillingWeights: GRES/gpu=15.0` — sustained 16-GPU use depresses fairshare; budget queue wait accordingly.

**Per-job shape** (must hit 16 GPU × 128 CPU simultaneously):

| `--cpus-per-task` | Concurrent tasks | GPUs used |
|---|---|---|
| 8 | 16 | **16** ✓ |
| 16 | 8 | 8 |
| 4 | 16 (CPU-limited at 64) | 16 ✓ |

Default GPU batch request: `--gres=gpu:1 --cpus-per-task=8 --mem=32G`.

**Target ≤ 90 array tasks** per submission (headroom for `salloc`). Bundle switches per task; never 1 switch = 1 array element.

### 7.2 P0 switch layout (36 tasks)

2 legs × 2 phases × 3 reps × 2 directions = 24 units × 100 snapshots:

| Leg | Snapshots/task | Tasks | ~Wall/task |
|-----|----------------|-------|------------|
| V106A (100 ps) | 100 | 12 | ~6 h |
| Y188L (500 ps) | 50 | 24 | ~15 h |
| **Total** | | **36** | |

Submit: `--array=0-35%16`, one job. ~500 GPU-h, ~1.5 days wall clock.

Full panel (~2,050 GPU-h) requires several 48 h arrays — **chain with `afterany`** so only one switch array is queued at a time (partial failures must not stall the chain; `.ok` markers + audit handle correctness).

### 7.3 Partition assignment

| Stage | Partition | Why |
|-------|-----------|-----|
| EM, extract, grompp, BAR, audit | **`normal`** | QOS: 2000 submits, 512 CPUs — free |
| Equil, switch | **`gpu`** | Only stages that need CUDA |
| Smoke / debug | **`dev`** | Separate QOS (2 GPUs, 2 h); does not eat the 100-submit quota |

### 7.4 Resumability (mandatory)

Switch tasks run 6–15 h against a 48 h cap with preemption possible. Per-switch marker files (`dgdl.xvg` + `status.json`) make requeue free.

- `sbatch`: `--requeue`, `--open-mode=append`
- `mdrun`: `-cpt 5` (equil + switch)
- Bundled switch tasks skip completed snapshots individually

### 7.5 Other

- Logs: `logs/pmx_neq/{em,equil,extract,switch}/`.
- **`prepare --force` must not delete completed stage outputs** (regenerates mdps/manifests in place).
- Implementation: `scripts/fep_pmx/submit_p0_neq.sh`, `config.py` (`switch_snapshots_per_task`).

---

## 8. What we are *not* doing

| Approach | Why dropped |
|----------|-------------|
| Perses / OpenMM fixed-λ FEP (full protein) | MBAR failure at high λ; pilot in `scripts/fep_jorgensen/README.md` |
| Truncated-sphere NEQ | Unnecessary under NEQ parallelism; boundary risk ([Genheden & Ryde 2012](https://doi.org/10.1021/ct200853g)) |
| Jorgensen MCPRO reproduction | Proprietary; approximated in `scripts/fep_jorgensen/` for comparison only |
| Dual-topology mutation | pmx uses hybrid residues, not full dual side chains |

---

## 9. Key references

DOIs verified via `https://doi.org/api/handles/` (Jul 2026).

1. Jarzynski, C. *Phys. Rev. Lett.* **78**, 2690 (1997). [doi:10.1103/PhysRevLett.78.2690](https://doi.org/10.1103/PhysRevLett.78.2690)  
2. Crooks, G. E. *Phys. Rev. E* **60**, 2721 (1999). [doi:10.1103/PhysRevE.60.2721](https://doi.org/10.1103/PhysRevE.60.2721)  
3. Bennett, C. H. *J. Comput. Phys.* **22**, 245 (1976). [doi:10.1016/0021-9991(76)90078-4](https://doi.org/10.1016/0021-9991(76)90078-4)  
4. Gapsys, V. et al. *J. Comput. Chem.* **36**, 348 (2015) — pmx hybrid topologies. [doi:10.1002/jcc.23804](https://doi.org/10.1002/jcc.23804)  
5. Aldeghi, M.; Gapsys, V.; de Groot, B. L. *ACS Cent. Sci.* **4**, 1708 (2018). [doi:10.1021/acscentsci.8b00717](https://doi.org/10.1021/acscentsci.8b00717)  
6. Aldeghi, M. et al. *ACS Cent. Sci.* **5**, 1468 (2019) — kinase inhibitor resistance NEQ. [doi:10.1021/acscentsci.9b00590](https://doi.org/10.1021/acscentsci.9b00590)  
7. Gapsys, V. et al. *JCTC* **8**, 2373 (2012) — soft-core electrostatics. [doi:10.1021/ct3000814](https://doi.org/10.1021/ct3000814)  
8. Rizzo, R. C. et al. *JACS* **122**, 12898 (2000) — HIV-RT resistance FEP cycle (MCPRO). [doi:10.1021/ja003113r](https://doi.org/10.1021/ja003113r)  
9. Smith, M. B. et al. *Bioorg. Med. Chem. Lett.* **18**, 969 (2008) — HIV-RT mutation energetics. [doi:10.1016/j.bmcl.2007.12.033](https://doi.org/10.1016/j.bmcl.2007.12.033)  
10. Aldeghi, M. *Methods Mol. Biol.* (2018) — NEQ tutorial. [doi:10.1007/978-1-4939-8736-8_2](https://doi.org/10.1007/978-1-4939-8736-8_2)  
11. Genheden, S.; Ryde, U. *JCTC* **8**, 1449 (2012) — system truncation artifacts. [doi:10.1021/ct200853g](https://doi.org/10.1021/ct200853g)  
12. Clark, A. J. et al. *J. Mol. Biol.* **431**, 1481 (2019) — co-alchemical ion for charge-changing mutations. [doi:10.1016/j.jmb.2019.02.003](https://doi.org/10.1016/j.jmb.2019.02.003)  
13. Gapsys, V. et al. *J. Comput. Chem.* **36**, 348 (2015) — double-system/single-box setup (Fig. 3). [doi:10.1002/jcc.23804](https://doi.org/10.1002/jcc.23804); NEQ + pmx implementation: Bhatt, V. et al. *JCTC* (2021). [doi:10.1021/acs.jctc.0c01045](https://doi.org/10.1021/acs.jctc.0c01045)
