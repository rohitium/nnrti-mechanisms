# Approach decision record

**Decision (2026-03):** pmx + GROMACS NEQ on **full solvated systems**.  
Supersedes truncated-sphere / OpenMM / Perses fixed-λ bake-off.

---

## Why NEQ won

| Problem | Equilibrium fixed-λ (Perses pilot) | NEQ (pmx) |
| --- | --- | --- |
| Intermediate states | Must converge λ=0.1…0.9 | **Only λ=0 and λ=1** equilibrated |
| Endpoint clashes | DELETE ghosts → MBAR blowups at λ≥0.8 | Bad switch = **outlier work**, discardable |
| Throughput | 11 windows × long prod each | **Array of 100 ps jobs** |
| Precedent for ΔΔG_bind | Zhang 2023 (PPI, AREX) | **Aldeghi 2018: 134 protein mutations, 27 ligands** |

---

## Key papers (most direct to our problem)

### Aldeghi, Gapsys & de Groot, ACS Cent. Sci. 2018

- **134 mutations, 17 proteins, 27 ligands** — exactly ΔΔG_bind upon protein mutation.  
- pmx hybrids + nonequilibrium work + BAR/Crooks.  
- RMSE **1.2 kcal/mol** full set; **0.8** on high-reproducibility subset.  
- Optimal protocol: ~**equal** equilibrium and switching time; chosen calibration used
  **5×3 ns equil** + **150×80 ps** switches (108 ns total per ΔΔG); production used 216 ns
  with 10 equil repeats for reproducibility.  
- Ion placement near protein hurt reproducibility — exclude protein interior from ion seeding.

### Aldeghi et al., ACS Cent. Sci. 2019 (kinase resistance)

- Same NEQ/pmx on **Abl:TKI** resistance (144 ΔΔG, 31 mutations, 8 inhibitors).  
- Direct analogue: predict which mutations weaken inhibitor binding.  
- Compared to Rosetta flex_ddg and ML; physics-based NEQ competitive for resistance ranking.

---

## Why truncated sphere was dropped

Truncation (Genheden 2012, Huang 2016, QresFEP-2) was a **cost workaround** for expensive
equilibrium multi-window FEP. Expert budget for NEQ full system:

- P0: ~110 GPU-h  
- Full panel: ~1,700 GPU-h (~2 weeks @ 20 GPUs)

At that cost, truncation saves little and introduces boundary/restraint artifacts. Apo pocket
collapse is the real scientific risk — truncation doesn't fix it.

---

## Alternatives considered (brief)

| Approach | Verdict |
| --- | --- |
| **Schrödinger FEP+** | Industry standard; Stanford site license on Sherlock — optional P0 benchmark |
| **Jorgensen MCPRO** | HIV RT precedent; licensed OPLS/MC — not implementable |
| **QresFEP-2** | Best cost/accuracy in Q; OPLS; different engine |
| **Perses AREX** | Valid but heavy; NEQ simpler for same endpoint problem |
| **OpenMM hybrid scaling** | No mature NEQ protein-mutation pipeline |
| **Amber TI** | pmemd.cuda licensing |
| **OpenFE** | Conda-centric; immature for protein muts |

---

## Panel traps (unchanged by method choice)

1. **P225H** — pmx `develop` branch for PRO ring opening.  
2. **K103N, G190E** — co-alchemical ion/water OR double-box; pick one.  
3. **Apo leg** — collapsing NNRTI pocket; identical weak restraints; rank > absolute ΔΔG.

---

## Perses pilot lesson (retained)

Full-protein fixed-λ Perses failed at λ≥0.8 due to DELETE-path ghost clashes on MBAR re-eval.
This is a **known failure class** of equilibrium alchemical FEP, not a one-off bug. NEQ avoids
sampling those states entirely.

See `scripts/fep_jorgensen/README.md` for diagnostic details.
