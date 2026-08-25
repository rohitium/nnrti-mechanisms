## Non-equilibrium alchemical free energy calculations

For each genotype, the change in DOR binding free energy was computed from a
thermodynamic cycle in which the p66 substitution was introduced alchemically in both
the DOR-bound and the ligand-free enzyme, ∆∆G_bind = ∆G_mut^holo − ∆G_mut^apo, so that
positive values correspond to weaker DOR binding. Multiply substituted genotypes were
decomposed into sequential single-residue legs (e.g. WT to K103N, then K103N to
K103N+M230L), whose values were summed and standard errors combined in quadrature.
Wild-type and mutant end states were connected by pmx hybrid residues,47 which share a
common backbone and carry dummy atoms only where the two side chains differ, built with
the amber14sbmut mutation force field to match the ff14SB parameters used for MD; DOR
parameters were carried over from the same OpenFF 2.0.0 assignment and exported to
GROMACS format through OpenFF Interchange. Hybrid systems were built in GROMACS 2023.1,48
solvated in a rhombic dodecahedral TIP3P box with 1.0 nm padding and neutralized at
0.15 M NaCl to match the MD prep, using particle-mesh Ewald electrostatics with 1.0 nm
cutoffs, LINCS constraints on bonds to hydrogen, a 2 fs timestep, and Gapsys linearized
soft-core alchemical interactions.49 Each replicate was seeded from the corresponding
equilibrated MD starting structure.

Non-equilibrium switching50,51 requires sampling only at the two physical end states.
Each end state (lambda = 0, wild type; lambda = 1, mutant) was minimized, minimized again
with the free-energy Hamiltonian active at that lambda so that dummy atoms relax, warmed
for 500 ps under a C-rescale barostat,52 and then sampled for 5 ns at 300 K and 1 bar
with velocity-rescaling temperature coupling and a Parrinello-Rahman barostat. From each
end-state trajectory, 100 evenly spaced frames were extracted after discarding the first
100 ps; each seeded an independent simulation in which lambda was driven linearly to the
opposite end state, giving 100 forward and 100 reverse work values per phase per
replicate. Switching times were 100 ps, or 500 ps for legs with widely separated forward
and reverse work distributions (V106A, Y188L, and the charge-changing WT to K103N leg
inherited by the K103N-containing genotypes); control runs on neutral legs gave ∆G values
invariant to switching time. Free energies were estimated with the Bennett acceptance
ratio53 as implemented in pmx (100 bootstrap resamples, 300 K), with the Crooks Gaussian
intersection and Jarzynski estimators computed as consistency checks alongside the
forward/reverse work overlap and hysteresis. Three independent replicates were run per
leg and phase, and the reported uncertainty is the standard error across replicates. For
the two legs that change protein net charge (WT to K103N, WT to G190E), the leading
Rocklin/Hunenberger finite-size periodicity correction was applied analytically post hoc;54
because holo and apo box volumes differ by less than 1%, this correction cancels to
below 10^-3 kcal/mol in ∆∆G_bind.
