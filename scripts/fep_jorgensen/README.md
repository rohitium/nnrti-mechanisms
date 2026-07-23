# Jorgensen FEP workflows

This directory now separates two scientifically different calculations.

## Exact Rizzo/Jorgensen (2000) reproduction

`exact_protocol.py` is the machine-readable contract for the protocol described
in `docs/Jorgensen-FEP-protocol.md`. It fixes all reported details:

- the published truncated HIV-1 RT binding-site residue set;
- CM1P-augmented OPLS-AA and the distance-dependent dielectric `epsilon=4r`;
- IMPACT c1.00, ten minimization cycles, 1 fs Verlet integration, SHAKE,
  Berendsen coupling at 0.2 ps, 3 ps at 100 K, 50 ps at 300 K, and six 4 ps
  quench blocks from 300 K to 50 K;
- a 22 A TIP4P water cap (about 850 waters);
- the published rigid/flexible MCPRO residue partitions, fixed protein
  backbone, and fully flexible inhibitor;
- 1,000,000 solvent-only equilibration, 10,000,000 full-equilibration, and
  10,000,000 averaging configurations per FEP window;
- protein-side-chain mutation in each inhibitor complex; and
- inhibitor-relative normalization to Sustiva.

Write the immutable protocol manifest with:

```bash
python - <<'PY'
from pathlib import Path
from scripts.fep_jorgensen.exact_protocol import ExactJorgensenProtocol
ExactJorgensenProtocol().write(Path("results/analysis/fep_jorgensen/exact_protocol.json"))
PY
```

The original calculation requires licensed/historical IMPACT c1.00 and MCPRO
1.65 executables plus the CM1P-augmented OPLS-AA parameter set and the original
inhibitor models. Those assets are not distributed in this repository. OpenMM,
Perses, Amber ff14SB, OpenFF, PME, TIP3P, Langevin dynamics, MBAR, or an apo leg
are not exact substitutes and must not be reported as an exact reproduction.

After MCPRO produces one protein-mutation free energy for every inhibitor,
provide a CSV with:

```text
inhibitor,delta_g_mutation_kcal_mol,uncertainty_kcal_mol
sustiva,...,...
nevirapine,...,...
mkc-442,...,...
9-cl-tibo,...,...
```

Then apply the paper's thermodynamic cycle:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze_exact mcpro_legs.csv \
  --mutation V106A --output v106a_relative_to_sustiva.csv
```

For inhibitor `i`, the reported value is
`DeltaG_mutation(i) - DeltaG_mutation(Sustiva)`. Sustiva is therefore 0.00 by
definition. The code propagates independent leg uncertainties in quadrature.

## OpenMM/Perses approximation retained from the source branch

The pre-existing `prepare.py`, `worker.py`, `panel.py`, and `analyze.py` workflow
is retained so earlier work is not destroyed. It computes explicit-solvent,
periodic holo-minus-apo mutation free energies with Perses/OpenMM and MBAR,
normalized to WT. It is mutation-agnostic and useful as a modern approximation,
but it does **not** exactly follow Jorgensen's protocol.

To run that approximation, generate its panel with:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.panel
```

Prepare hybrid systems locally with Perses and OpenMMTools, run the generated
OpenMM manifest on Sherlock, and analyze locally with:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --all-targets
```
