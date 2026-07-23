# Jorgensen-inspired approximate mutation FEP

This workflow consolidates the earlier Perses/OpenMM panel work and the exact-protocol
specification into a single branch:

**OpenMM/openmmtools MD-equilibrated MCMC alchemical mutation FEP, Jorgensen-cycle inspired.**

That is *not* an exact MCPRO reproduction.  It keeps the paper's scientific shape:

1. start from the 19 manuscript DOR-bound complexes;
2. equilibrate each complex with OpenMM MD using a Jorgensen-inspired schedule;
3. build **holo-only** protein-side-chain alchemical legs;
4. sample each λ window with fixed-lambda OpenMM workers or optional openmmtools MCMC;
5. estimate `ΔG_mutation` with MBAR;
6. normalize mutation free energies relative to a reference system (default: WT).

## Branch

All FEP work lives on **`jorgensen-fep`**.  Older branches
(`agent/exact-jorgensen-fep`, `codex/fep-jorgensen-v106a`) are retired.

## Three execution tiers

| Tier | Purpose | Entry point |
| --- | --- | --- |
| Exact MCPRO | Only valid exact reproduction | `exact_protocol.py`, `analyze_exact.py` |
| Approximate OpenMM | Default merged workflow | `prepare.py`, `worker.py`, `analyze.py` |
| Optional MCMC | Local/Sherlock path with openmmtools replica exchange | `mcmc_sample.py` |

## Approximate workflow

Write the protocol manifest:

```bash
python - <<'PY'
from pathlib import Path
from scripts.fep_jorgensen.approx_protocol import ApproxJorgensenProtocol
ApproxJorgensenProtocol().write(Path("results/analysis/fep_jorgensen/approx_protocol.json"))
PY
```

Plan the full 19-leg panel:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.panel
```

This writes:

- `equilibrate_all.sh` — Jorgensen-inspired MD on all 19 complexes
- `prepare_all.sh` — equilibration + Perses hybrid setup (or use `--skip-equilibration` if reusing MD start structures)
- `worker_manifest.csv` — 209 holo λ-window tasks (19 legs × 11 states)

Prepare one leg locally (`nnrti-fep` env: Perses + OpenEye required):

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A
```

This runs, in order:

1. `equilibrate.py` — Jorgensen-inspired minimization, 100 K → 300 K equilibration, quench blocks;
2. Perses hybrid construction for the bound-complex mutation leg;
3. serialization of `holo/hybrid_system.xml` and `holo/schedule.json` for Sherlock workers.

Run fixed-λ windows on Sherlock:

```bash
bash scripts/sherlock/submit_fep_jorgensen_windows.sh
```

Analyze locally:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --all-targets
```

Outputs:

- `manuscript_panel_summary.csv` — absolute `ΔG_mutation` per target
- `manuscript_panel_relative.csv` — values relative to WT (WT = 0 by definition)

Optional integrated openmmtools MCMC on one prepared leg:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.mcmc_sample --mutation V106A
```

## Exact MCPRO boundary

The exact Rizzo/Jorgensen (2000) contract remains in `exact_protocol.py` and
`docs/Jorgensen-FEP-protocol.md`.  After licensed MCPRO produces inhibitor legs, use:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze_exact mcpro_legs.csv \
  --mutation V106A --output v106a_relative_to_sustiva.csv
```

Do not label OpenMM/Perses/OpenFF results as exact MCPRO reproduction.

## What changed from `codex/fep-jorgensen-v106a`

- Removed default holo-minus-apo thermodynamic cycle.
- Added Jorgensen-inspired MD equilibration before alchemical setup.
- Analysis now reports `ΔG_mutation` in the inhibitor-bound complex.
- Added reference normalization mirroring the paper's inhibitor-relative cycle, with WT
  as the manuscript reference instead of Sustiva.
