# Jorgensen-inspired approximate mutation FEP

**OpenMM/openmmtools MD-equilibrated MCMC alchemical mutation FEP, Jorgensen-cycle inspired.**

## Prepare backends

| Backend | Physics | When to use |
| --- | --- | --- |
| **`perses` (default)** | Hybrid topology; bonded + nonbonded λ staging | Production FEP |
| **`scaling`** | Single topology; nonbonded scaling only | Fallback without Perses |

## One-time env setup

Perses requires `numpy<2.4` for openmmtools/numba compatibility and AmberTools GAFF
(`ambertools` from conda-forge). Ligand/residue handling uses an in-repo OpenEye shim
(`openeye_shim.py`) backed by RDKit + OpenFF — no licensed OpenEye install required.
Endstate validation is skipped because the shim objects are not deep-copyable.

If you have a user-site NumPy in `~/.local`, prefix commands with `PYTHONNOUSERSITE=1`:

```bash
bash scripts/fep_jorgensen/setup_perses_env.sh
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A
```

On macOS without OpenEye, hybrid prep still works via the RDKit/OpenFF shim:

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A
```

If you only need a quick local approximation, use scaling:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A --backend scaling
```

Or create an isolated env:

```bash
conda env create -f envs/nnrti-fep.yml
conda activate nnrti-fep
```

## Environments

| Step | Env |
| --- | --- |
| `prepare.py` (Perses hybrid) | `nnrti-prep` after setup script |
| `analyze.py` | `nnrti-prep` or `nnrti-openmm` |
| `worker.py` on Sherlock | `nnrti-openmm` |

Perses hybrid prep extracts unsolvated protein+ligand coordinates from your existing MD
start PDBs, then builds a solvated hybrid system with proper topology changes. It does **not**
reuse the serialized MD `system.xml` directly — that file is endpoint-specific, not hybrid.

## Workflow

See **[docs/fep_jorgensen_sherlock.md](../docs/fep_jorgensen_sherlock.md)** for the full
local prep → Sherlock GPU pilot → batch → local analysis runbook (WT→V106A pilot included).

```bash
bash scripts/fep_jorgensen/setup_perses_env.sh   # once
PYTHONPATH=. python -m scripts.fep_jorgensen.panel
bash results/analysis/fep_jorgensen/prepare_all.sh
bash scripts/sherlock/submit_fep_jorgensen_windows.sh
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --all-targets
```

Single leg (V106A):

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.panel --mutation V106A
# rsync legs/wt_to_V106A/ to Sherlock, then:
#   bash scripts/sherlock/salloc_fep_jorgensen_gpu.sh
#   bash scripts/sherlock/run_fep_jorgensen_pilot.sh
#   ./scripts/sherlock/submit_fep_jorgensen_v106a.sh
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --target V106A
```

## Scaling fallback

The `scaling` backend copies the start genotype MD `system.xml` and scales side-chain
nonbonded interactions. Faster to set up, but ghost-atom bonded artifacts remain. See
`alchemical.py`.

## Three execution tiers

| Tier | Purpose | Entry point |
| --- | --- | --- |
| Exact MCPRO | Only valid exact reproduction | `exact_protocol.py`, `analyze_exact.py` |
| Perses hybrid FEP | Default production path | `prepare.py --backend perses` |
| Scaling fallback | MD-asset nonbonded approximation | `prepare.py --backend scaling` |

## Exact MCPRO boundary

See `exact_protocol.py` and `docs/Jorgensen-FEP-protocol.md`. Do not label OpenMM/Perses
results as exact MCPRO reproduction.
