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

PRO→X mutations (e.g. P225H) require a runtime patch (`perses_patches.py`) because
Perses omits PRO from its amino-acid registry but still asserts on it during charge
handling.

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
| `analyze.py` | `nnrti-prep` (local) |
| `worker.py` on Sherlock | `module load chemistry py-openmm/8.1.1_py312` |

Perses hybrid prep uses your existing MD assets:

| Phase | Start structure | Endpoint validation |
| --- | --- | --- |
| **holo** | `results/md_runs/{start}/.../*_md_rep01_start.pdb` | `{end}` holo start PDB |
| **apo** | `results/md_runs/apo/{start}/.../*_apo_md_rep01_start.pdb` | `{end}` apo start PDB |

`prepare.py --phase all` (default) builds both Perses hybrids. Analysis reports
**ΔΔG_bind = ΔG_mut(holo) − ΔG_mut(apo)** when apo windows exist.

## Sampling tiers

| Tier | Sampler | Where | Sherlock deps |
| --- | --- | --- | --- |
| **A. Fixed-λ MD windows** | `worker.py` | Sherlock GPU | `module load py-openmm` only ✓ |
| **B. Replica-exchange MCMC** | `mcmc_sample.py` | **Local Mac** (`nnrti-prep`) | perses + openmmtools — **not** on Sherlock without a custom venv |
| **C. Exact MCPRO** | licensed MCPRO | N/A | Documented only |

Sherlock policy discourages conda. The validated production path is **Tier A** on Sherlock
(prep + analyze local, λ windows on GPU). Tier B MCMC requires the full Perses stack; we do
not currently support running it on Sherlock's `py-openmm` module alone. If you need MCMC
later, options are: local GPU with `nnrti-prep`, or a `$GROUP_HOME` venv (untested on Sherlock).

Check convergence before trusting MBAR:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.convergence_cli --leg-dir results/analysis/fep_jorgensen/legs/wt_to_V106A
```

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
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A --phase all
PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.panel --mutation V106A
# manifest now has holo + apo tasks (22 windows/leg). rsync full leg dir to Sherlock.
#   SHERLOCK_USER=rsatija bash scripts/rsync_fep_jorgensen.sh push V106A
#   ./scripts/sherlock/submit_fep_jorgensen_v106a.sh          # holo tasks 0-10
#   ./scripts/sherlock/submit_fep_jorgensen_v106a_apo.sh      # apo tasks 11-21
PYTHONNOUSERSITE=1 PYTHONPATH=. ~/miniconda3/envs/nnrti-prep/bin/python -m scripts.fep_jorgensen.analyze --target V106A
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
