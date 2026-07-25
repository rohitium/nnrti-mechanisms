# Jorgensen FEP Protocol

Reference: Rizzo, R. C., Wang, D. P., Tirado-Rives, J., & Jorgensen, W. L. (2000). Validation of a model for the complex of HIV-1 reverse transcriptase with sustiva through computation of resistance profiles. Journal of the American Chemical Society, 122(51), 12898-12900.

For the merged OpenMM/openmmtools approximation used in this repository, see
[`docs/Jorgensen-approx-FEP-protocol.md`](Jorgensen-approx-FEP-protocol.md) and
`scripts/fep_jorgensen/README.md`.

# FEP Protocol for Table 2 Relative Fold Resistance Energies

## Overview

The protocol computes relative fold resistance (ΔΔG) using a thermodynamic cycle in which the **protein side chain is mutated** (V106→A or Y181→C) in the presence of each of two drugs, rather than mutating the drugs. Values in Table 2 are normalized to Sustiva, so the reported ΔΔG for each inhibitor is the difference between that inhibitor's mutation free energy and Sustiva's.

## Step 1 — Construction of Binding-Site Models

- Sustiva model: docked into the NNRTI site using the MATADOR program, starting from the 2.55‑Å MKC‑442/HIVRT crystal structure (pdb 1rt1) with MKC‑442 removed; only residues within ~15 Å of MKC‑442 were retained.
- The comparison complexes (nevirapine, MKC‑442, 9‑Cl TIBO) were built analogously from their X‑ray structures: nevirapine (1vrt), HEPT/MKC‑442 (1rti/1rt1), 9‑Cl TIBO (1rev).
- Protein residues in the binding‑site model: 91–110A, 161–205A, 222–242A, 316–321A, 343–349A, 381–383A, and 134–140B.
- Initial Sustiva geometry and nonbonded energies came from the **CM1P‑augmented OPLS‑AA** force field; nonbonded energies were stored on a spherical grid, and a distance‑dependent dielectric (ε = 4r) was used.

## Step 2 — Molecular Dynamics Equilibration (IMPACT)

Each docked complex was relaxed by MD to allow backbone and side chains to adjust, using the CM1P‑augmented OPLS‑AA force field in the **IMPACT program (Version c1.00, Schrödinger)**:

- Ten cycles of gradient‑based energy minimization preceded MD.
- **Restraint scheme:** residues within ~10 Å of the binding site moved freely (95–107A, 172A, 177–182A, 188–192A, 198A, 227A, 229A, 234–236A, 318–319A, 321A, 135–139B); residues in the 10–12 Å shell were harmonically restrained (94A, 108A, 175–176A, 183A, 187A, 225A, 237–239A, 317A, 320A, 349A, 382–383A, 134B, 140B); all others were fixed at their conjugate‑gradient‑minimized positions.
- **Integrator:** Verlet algorithm, time step 0.001 ps.
- **Thermostat:** Berendsen bath, coupling/relaxation parameter 0.2 ps.
- Bond lengths constrained by SHAKE; distance‑dependent dielectric ε = 4r.
- **Schedule:** 3 ps initial equilibration at 100 K → 50 ps equilibration at 300 K → quenching from 300 K to 50 K over 6 blocks of 4 ps each.
- The identical MD protocol was applied to the nevirapine, MKC‑442, and 9‑Cl TIBO complexes; the resulting structures were passed to the MC simulations.

## Step 3 — Monte Carlo / Free Energy Perturbation (MCPRO)

MC/FEP simulations were run with **MCPRO Version 1.65 (Jorgensen, Yale)**, CM1P‑augmented OPLS‑AA force field:

- Each complex was briefly energy‑minimized before MC (distance‑dependent dielectric ε = 4r).
- **Solvation:** 22‑Å water cap containing ~850 TIP4P water molecules.
- **System partitioning** into rigid residues (91–94A, 109–110A, 116–178A, 184–185A, 192–197A, 199–205A, 222–224A, 230–232A, 240–242A, 316–317A, 320–321A, 343–349A, 381–383A, 134–135B, 137B, 140B) and flexible residues (95–108A, 179–183A, 186–191A, 198A, 225–229A, 233–239A, 318–319A, 136B, 138B).
- All HIVRT side chains within ~10 Å of the water‑cap center were sampled; the **protein backbone was fixed**, and each inhibitor was fully flexible.
- **Sampling per FEP window:** 1 million configurations of solvent‑only equilibration, 10 million of full equilibration, and 10 million of averaging.

## Step 4 — Free Energy Calculation via Thermodynamic Cycle

The mutation was carried out as an alchemical perturbation of the wild‑type side chain into the mutant side chain (Val106→Ala or Tyr181→Cys) in the presence of each drug (Figure 3 cycle):

- For two inhibitors A and B, ΔG_WT and ΔG_MUT are the B‑vs‑A binding free‑energy differences with wild‑type and mutant protein; ΔG_A and ΔG_B are the mutant‑vs‑wild‑type changes for each drug.
- The identity used is: **ΔG_MUT − ΔG_WT = ΔG_B − ΔG_A = ΔΔG**, the experimentally observable difference in fold resistance (RT ln FR_B − RT ln FR_A).
- Table 2 reports ΔΔG for each inhibitor **normalized to Sustiva** (Sustiva = 0.00), for both the Y181C and V106A mutations, with associated statistical uncertainties (±0.3–0.5 kcal/mol).