# Jorgensen-inspired approximate mutation FEP

**OpenMM/openmmtools MD-equilibrated MCMC alchemical mutation FEP, Jorgensen-cycle inspired.**

## Why this is not hard anymore

The manuscript MD pipeline already produced, for every genotype:

- an equilibrated holo `*_md_rep01_start.pdb`
- a matching serialized `*_md_rep01_system.xml`
- often 100 ns production trajectories

FEP does **not** need to rebuild those systems. It needs an alchemical path between the
**start** and **end** genotype in the bound complex. Prepare therefore:

1. copies the existing start-genotype OpenMM system;
2. diffs the start/end holo PDBs to locate the mutated residue;
3. uses Amber ff14SB residue templates to pick side-chain atoms for nonbonded scaling.

No Perses, OpenEye, or RDKit is required. RDKit is for small-molecule chemistry; protein
mutation mapping uses OpenMM plus the endpoint structures we already prepared.

## Environments

| Step | Env |
| --- | --- |
| `prepare.py`, local analysis | `nnrti-prep` |
| `worker.py` on Sherlock | `nnrti-openmm` |

## Workflow

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.panel
bash results/analysis/fep_jorgensen/prepare_all.sh   # or one leg below
bash scripts/sherlock/submit_fep_jorgensen_windows.sh
PYTHONPATH=. python -m scripts.fep_jorgensen.analyze --all-targets
```

Single leg:

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.prepare --mutation V106A
```

Optional extra Jorgensen-inspired equilibration (usually unnecessary because MD starts are
already equilibrated):

```bash
PYTHONPATH=. python -m scripts.fep_jorgensen.equilibrate --mutation V106A
```

## Alchemical strategies

| Strategy | When | System used |
| --- | --- | --- |
| `annihilate_wt_sidechain` | Most mutations | start genotype MD system |
| `annihilate_mutant_sidechain` | Side-chain growth (e.g. GLY→ALA) | end genotype MD system, sign flipped in analysis |
| `annihilate_shared_sidechain` | Rare parameter-only fallback | start genotype MD system |

## Three execution tiers

| Tier | Purpose | Entry point |
| --- | --- | --- |
| Exact MCPRO | Only valid exact reproduction | `exact_protocol.py`, `analyze_exact.py` |
| Approximate OpenMM | Default merged workflow | `prepare.py`, `worker.py`, `analyze.py` |
| Optional MCMC | openmmtools replica exchange (requires Perses hybrid rebuild) | `mcmc_sample.py` (legacy) |

## Exact MCPRO boundary

See `exact_protocol.py` and `docs/Jorgensen-FEP-protocol.md`. Do not label OpenMM results as
exact MCPRO reproduction.
