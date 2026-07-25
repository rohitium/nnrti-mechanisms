# Jorgensen-inspired approximate mutation FEP

This document describes the default workflow on the **`jorgensen-fep`** branch.
It supersedes the earlier holo-minus-apo Perses panel.

## Scientific label

**OpenMM/openmmtools MD-equilibrated MCMC alchemical mutation FEP, Jorgensen-cycle inspired.**

## Protocol

1. **Systems** — 19 manuscript DOR-bound RT complexes from `results/md_runs/`.
2. **MD equilibration** — OpenMM relaxation inspired by Rizzo/Jorgensen (2000):
   minimization, 3 ps at 100 K, 50 ps at 300 K, six 4 ps quench blocks from 300 K to 50 K,
   then reheat to 300 K.  Force field: AMBER ff14SB + OpenFF-2.0.0 + TIP3P/PME.
3. **Alchemical leg** — Perses point-mutation hybrid in the **bound complex only**.
4. **Sampling** — default: fixed-λ OpenMM windows with multistate energy reevaluation and MBAR;
   optional: openmmtools GHMC replica exchange via `mcmc_sample.py`.
5. **Analysis** — `ΔG_mutation` per leg/target with MBAR; compound mutants sum sequential legs.
6. **Normalization** — relative to WT by default, analogous to the paper's Sustiva reference.

## Exact reproduction boundary

Exact MCPRO/IMPACT execution remains documented in `docs/Jorgensen-FEP-protocol.md` and
enforced by `scripts/fep_jorgensen/exact_protocol.py`.  OpenMM results must not be reported
as exact reproduction unless those licensed tools and parameter sets are actually used.

## References

- [OpenMM alchemical free-energy tutorial](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Alchemical_free_energy_calculations.html)
- [openmmtools alchemical transformations](https://openmmtools.readthedocs.io/en/stable/gettingstarted.html)
- [openmmtools MCMC framework](https://openmmtools.readthedocs.io/en/stable/mcmc.html)
