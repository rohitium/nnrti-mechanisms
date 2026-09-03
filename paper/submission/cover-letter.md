[Date]

Prof. Kenneth M. Merz, Jr.
Editor-in-Chief, *Journal of Chemical Information and Modeling*
Department of Chemistry, Michigan State University

Dear Prof. Merz,

Please consider the enclosed manuscript, "Molecular Simulation of Doravirine
Resistance in HIV-1 Reverse Transcriptase," for publication as an Article in the
Pharmaceutical Modeling section of *JCIM*.

Doravirine is a second-generation non-nucleoside reverse transcriptase inhibitor
designed to tolerate the substitutions that defeat earlier drugs in its class. It
is now part of the standard of care, but phenotypic susceptibility data exist for
only a fraction of the genotypes clinicians encounter, which complicates its use
in people with pre-existing NNRTI resistance. We asked a question that this gap
makes practical rather than academic: do physics-based binding free energy
calculations track measured drug susceptibility well enough to help fill it?

To answer it we simulated wild-type reverse transcriptase and 18
doravirine-associated genotypes in explicit solvent, three independent 100 ns
replicates each, and computed the effect of every substitution on drug binding by
two independent routes — non-equilibrium alchemical free energy perturbation and
wild-type-referenced MM/GBSA — then compared both against median phenotypic
fold-changes curated from the Stanford HIV Drug Resistance Database.

Three things make this a fit for *JCIM* rather than a virology journal:

**A matched benchmark, not a case study.** Previous simulation studies of NNRTI
resistance examine up to about four substitutions, rarely in combination. Our
panel of 18 includes the multi-mutation patterns that dominate clinical doravirine
resistance, so the correspondence between a computed quantity and a measured one
can be quantified across a dataset rather than illustrated on a favorable
example.

**A result reported as it came out.** ΔΔG_bind correlates only weakly with
log-transformed fold-change (R² = 0.26, p = 0.07 across the 13 genotypes with
established phenotypes, and more weakly over the full panel), and MM/GBSA
identifies resistant genotypes but not susceptible ones. We conclude that neither
method is yet a proxy for in vitro susceptibility testing. A bounded negative
result on a clinically important target seems to us more useful than another
selected success, and it is the kind of finding a methods journal is the right
place to publish.

**Structural mechanism from the same trajectories.** The simulations rationalize
the drug's known resistance profile: the hydrogen bond doravirine accepts from the
residue 103 main-chain carbonyl is preserved across all four K103N-containing
genotypes even as the side chain at that position moves about 3 Å closer,
while substitutions at Val106 displace the drug in opposite directions depending
on whether the side chain is shortened or enlarged.

We note the Journal's exclusion of straightforward docking applications to a
single target without experimental validation. This work uses alchemical free
energy and end-point methods rather than docking, and every computed value is
tested against measured phenotypic data.

All code and derived data needed to reproduce every figure and table are openly
available at https://github.com/rohitium/nnrti-mechanisms under the MIT license,
and per-replicate values for every reported quantity are given in the Supporting
Information.

This manuscript is original, has not been published previously, and is not under
consideration elsewhere. All authors have read and approved the submission and
declare no competing financial interest.

Thank you for your consideration.

Sincerely,

Rohit Satija
Division of Infectious Diseases, School of Medicine
Stanford University, Stanford, California 94305, United States
rsatija@stanford.edu

on behalf of Kaiming Tao and Robert W. Shafer
