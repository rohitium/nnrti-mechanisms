# Baseline: wt_to_K103N at 100 ps switches (the panel result as of 2026-08-15)

**This is a COPY, not a retirement.** The originals remain live under
`results/analysis/fep_pmx/legs/wt_to_K103N/` and are what `panel_ddg.csv` reports.
Copied here because the 500 ps re-run overwrites `analysis.json` in place, and
this is the arm we compare against.

## The result being preserved

ΔΔG_bind = **+0.92 ± 2.19 kcal/mol** (n = 3), experimental DOR fold = 1.2.

Per-replicate ΔΔG_bind: **+5.09, −0.01, −2.33** → σ_DDG = **3.80**, SEM = 2.19.

| phase | rep | bar_dg | boot err | cgi_dg | jarz_dg | hysteresis | separation |
|---|---|---|---|---|---|---|---|
| holo | 1 | −75.71 | 0.36 | −74.76 | −75.84 | 14.19 | 4.24 σ |
| holo | 2 | −77.99 | 1.27 | −76.52 | −78.32 | 18.54 | 4.15 σ |
| holo | 3 | −79.79 | 0.37 | −78.12 | −79.80 | 19.36 | 5.07 σ |
| apo | 1 | −80.80 | 0.52 | −79.89 | −80.91 | 15.38 | 5.05 σ |
| apo | 2 | −77.98 | 0.59 | −78.16 | −78.18 | 14.71 | 4.15 σ |
| apo | 3 | −77.46 | 0.27 | −77.67 | −77.40 | 15.21 | 3.65 σ |

`hysteresis` = ⟨W_f⟩ − ⟨W_r⟩ in kcal/mol; `separation` = hysteresis / pooled work
sd — the dimensionless Crooks overlap measure. BAR is comfortable at ≲ 2 σ; this
leg sits at 3.7–5.1 σ, where BAR is both noisy and biased.

Protocol: 100 ps switches, 100 snapshots/endpoint/rep, 5 ns equilibration, 3 reps,
raw non-neutral box (delta_q = −1) + Rocklin/Hünenberger analytical correction.

## Why this leg is worth re-running

`wt_to_K103N` is the single most expensive item in the panel's error budget. All
four K103N genotypes (K103N, K103N+M230L, K103N+P225H, L100I+K103N) sum this leg,
so they all inherit its σ_DDG = 3.80 — their own second legs are quiet by
comparison (`K103N_to_K103N_M230L` σ_DDG = 0.59).

Getting every genotype under SEM 1 by replicates alone costs 29 leg-reps panel-wide;
**19 of those 29 are this leg** (3 → 22 reps ≈ 152 GPU array elements). A 500 ps
re-run at 3 reps costs ~36 elements. If it works, it is ~4× cheaper and fixes four
genotypes at once.

## The hypothesis being tested

Per-replicate ΔG tracks per-replicate dissipation almost 1:1 in the holo phase:

```
hyst 14.19 -> bar_dg -75.71
hyst 18.54 -> bar_dg -77.99
hyst 19.36 -> bar_dg -79.79
```

If BAR had fully removed the dissipation bias these would be independent. So part
of σ_DDG = 3.80 looks like *varying residual bias* rather than pure conformational
sampling. Dissipation falls ~linearly with switch time while work spread falls as
its square root, so 100 → 500 ps should improve separation by ~√5 (≈4.5 σ → ≈2 σ),
which is where BAR becomes trustworthy.

Note this contradicts the panel's earlier switch-length-invariance finding — that
was measured on **neutral** legs (V106A, F227C) dissipating 1–4 kcal/mol, which had
nothing to gain. This leg dissipates 14–19.

## OUTCOME (2026-08-17): the 500 ps re-run worked

| | 100 ps (this archive) | 500 ps |
|---|---|---|
| ddG_bind | +0.92 +- 2.19 | **+0.54 +- 0.23** |
| per-rep ddG | +5.09, -0.01, -2.33 | +0.09, +0.82, +0.71 |
| sigma_DDG | 3.80 | **0.40** |
| hysteresis holo | 14.19 / 18.54 / 19.36 | 12.04 / 12.01 / 8.68 |
| hysteresis apo | 15.38 / 14.71 / 15.21 | 11.47 / 12.63 / 9.01 |
| separation | 3.65-5.07 sigma | 3.83-5.35 sigma |
| Crooks overlap | ~0 | ~0 (all 6 units still flagged) |

All four K103N genotypes dropped below SEM 1: K103N 2.19 -> 0.23,
K103N+M230L 2.22 -> 0.41, K103N+P225H 2.32 -> 0.79, L100I+K103N 2.28 -> 0.68.

### The predicted mechanism was WRONG; the real one is holo/apo cancellation

The prediction was that dissipation would fall ~5x and separation improve by
~sqrt(5) to ~2 sigma. Hysteresis fell only ~35% and **separation did not improve
at all** (marginally higher). The work distributions narrowed roughly *in
proportion* to the hysteresis, so their width is not dissipation-driven -- that
would give sigma ~ sqrt(W_diss) and a visible separation gain. It is set by
conformational heterogeneity of the endpoint ensembles, which longer switches do
not touch. Hence overlap stayed ~0.

What actually changed is that the residual dissipation bias became **correlated
between holo and apo**, so it cancels in the double difference:

```
holo - apo hysteresis, per replicate
100 ps:  -1.20, +3.83, +4.16   -> sd 3.00
500 ps:  +0.57, -0.61, -0.33   -> sd 0.62
```

That difference-spread fell 4.8x; sigma_DDG fell 9.6x. Same effect seen in G190S,
where per-phase variance is large (~4) but ddG variance is tiny (0.59).

**So the lever is holo/apo correlation, not overlap or dissipation magnitude.**
Across 14 legs, sd(holo-apo hysteresis) correlates with sigma_DDG at r = 0.61.

Also validated: the Rocklin/Hunenberger charge correction came out at ~1e-5
kcal/mol per replicate, exactly as `charge_correction.py` argues it should -- the
finite-size terms cancel between holo and apo in matched boxes.

### Where 500 ps is worth it (and where it is not)

The `git clean` left these legs with no surviving equil, so a 500 ps re-run means
a full rebuild: ~36 GPU array elements per leg. Compare against replicates
(8 elements per leg-rep, SEM = sigma_DDG/sqrt(n)):

- **K103N family** needed 19 leg-reps (~152 elements) -- 500 ps at 36 was a bargain.
- **The remaining 5 genotypes** need only 7 leg-reps (~56 elements); 500 ps on
  their 5 legs would cost ~180. Use replicates there, which are also *certain*
  rather than resting on an r = 0.61 predictor.

## What the re-run cannot tell us

The `git clean -fd` of 2026-08-13 destroyed this leg's equil snapshots, so the
500 ps arm must re-equilibrate from scratch. Switch length and endpoint ensemble
therefore change **together** — a drop in σ_DDG cannot be attributed to switch
length alone. That attribution is what the G190E test in
`scripts/fep_pmx/RUNBOOK_G190E.md` step 11 is for: G190E's equil is being built
now, so its 500 ps switches can reuse the identical snapshots.

Also: comparing two σ estimates each on 2 degrees of freedom needs a ratio near
19× to be significant, so σ would have to fall below ~0.9 to be *statistically*
distinguishable from 3.80. The practical bar is different and is the one that
matters — does the new SEM land under 1.
