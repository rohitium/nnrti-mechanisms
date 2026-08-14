# FEP Results talking points (incomplete panel)

Source numbers: `results/analysis/fep_pmx/panel_ddg.csv`,
`panel_discussion_tiers.csv`, `panel_ddg_vs_experiment.png`.
Regen scatter/tiers: `python scripts/fep_pmx/combine_neq.py --replot-only`.

## Frame

- FEP panel is **incomplete**. Confidence = replicate SEM + BAR/CGI/Jarzynski
  agreement, **not** experimental match (`scripts/fep_pmx/STATUS.md`).
- n = 14 scatter vs log10(fold): **Pearson R² = 0.09, p = 0.28** (no fitted line
  on the figure). Weak correlation is a finding, stated once — not a success claim.
- Protocol trust comes from the V106A walkthrough
  (`results/analysis/fep_pmx/protocol/V106A/` / `protocol_v106a/`).

## Main text (SEM ≤ 0.6 kcal/mol)

| genotype | ΔΔG_bind | SEM | fold |
|---|---:|---:|---:|
| G190A | +0.27 | 0.17 | 2.7 |
| G190S | +2.01 | 0.34 | 5.2 |
| V106A | +1.76 | 0.51 | 9.6 |
| Y318F | +1.41 | 0.45 | 11 |
| Y188L | +4.52 | 0.49 | 149 |
| V106M | +6.10 | 0.16 | 3.4 |

V106M: tight and reproducible, but binding ΔΔG is much larger than the modest fold
implies — a **binding-vs-phenotype** finding (do not call it a pipeline
“overprediction” failure).

## Show on scatter; do not over-interpret

V106I, V106A+F227L/L234I/P225H, A98G+F227C, V106I+F227C, Y181C, F227C
(SEM ~0.7–1.2 or endpoint-limited).

## Omit from main-text point estimates

- K103N (SEM 2.19; charge leg)
- G190E (charge protocol)
- Missing entirely: K103N+M230L, K103N+L100I, K103N+P225H — say so.

## MM/GBSA

One paragraph + SI table only. See
`results/analysis/binding_energy/MMGBSA_METHOD_AND_RECOMPUTE.md`
(“Manuscript role”). Strip hybrid-topology / 5 ns language from that section.
