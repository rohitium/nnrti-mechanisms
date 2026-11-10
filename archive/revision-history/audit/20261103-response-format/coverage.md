# Response coverage audit — 3 November 2026

Baseline: current shared Dropbox Word files matched the local submission package before editing. Compared revised manuscript and Supplementary Text with the actual submitted DOCX files. Original reviewer wording and all existing native citation fields and drawings were preserved. Only Response-to-Reviewers.docx was edited in submission-ready/.

Locations in the response use section/subsection names and exact paragraph openings rather than page numbers, which can change in Word. Paragraph indices below are internal, zero-based indices in python-docx, including empty paragraphs.

## Manuscript changes

| Paragraph | Revised passage | Response coverage |
|---|---|---|
| 9 | Results: Using a cutoff of ΔΔGbind = +0.5 kcal/mol, approximately the thermal energy (kBT) at 300 K, FEP correctly predicted littl | R3.1 |
| 19 | Figure 1: Crystal structure of Doravirine (DOR) bound to wild-type (WT) Reverse Transcriptase (RT). (A) 3D structure of WT RT in c | R1.4 |
| 23 | The starting holo structure was the crystal structure of HIV-1 RT bound to doravirine (DOR; PDB ligand residue 2KW) from PDB 4NCG  | R1.2 |
| 25 | Equilibrium MD simulations were performed with OpenMM 8.5.2.36 Protein parameters were assigned from the Amber14 OpenMM force-fiel | R2.5 |
| 26 | Each RT-DOR complex was simulated in three independent 100 ns production replicates. Systems underwent a short pre-production mini | R1.2 |
| 31 | Ring planes were obtained by singular value decomposition of the ring atom coordinates, the interplanar angle between two rings wa | R1.2 |
| 32 | Per-residue Cα root-mean-square fluctuation (RMSF) was computed after superposing every frame on the Cα atoms of the sixteen NNIBP | R1.4 |
| 33 | For experimental structure comparisons, sixteen NNIBP-lining Cα atoms were superposed between the DOR-bound structures (Protein Da | R1 general |
| 34 | Convergence tests | R1.2 / R2.1 |
| 35 | To test whether the sampled DOR-bound state is recovered starting from gradually extreme perturbations to the starting configurati | R1.2 |
| 36 | For a series of N successive frames of a structural observable, with values xt and mean x̄, sampled at interval Δt, the normalized | R2.1 |
| 38 | For each replicate trajectory, MM/GBSA energy components were computed from snapshots of the explicit-solvent simulations without  | Other revisions (reference numbering, SI inventory and cross-references) |
| 42 | Here, ∆EvdW and ∆Eelec are the van der Waals and electrostatic molecular-mechanics contributions evaluated by selectively zeroing  | R3.3 |
| 43 | Replicate means are reported from 100 snapshots spaced evenly across the post-equilibration portion (final 75 ns) of each producti | R3.3 (retained protocol; paragraph split) |
| 56 | where kB is the Boltzmann constant, T is the Temperature, and H(r, p, λ) the Hamiltonian of the system at coupling parameter λ. Wi | R2.5 |
| 58 | with λ = 0 the wild type and λ = 1 the mutant residue. Hybrid systems were built in GROMACS 2023.1,53 solvated in a rhombic dodeca | R2.4 |
| 59 | Non-equilibrium switching55,56 requires sampling only at the two physical end states. Each end state was minimized with the free-e | Other revisions (reference numbering, SI inventory and cross-references) |
| 61 | Switching times were 100 ps, or 500 ps for legs with widely separated forward and reverse work distributions (e.g. charge-changing | R2.1 / Other revisions (estimator explanation) |
| 63 | while the Crooks fluctuation theorem56 relates the forward and reverse work distributions PF(W) and PR(W): | Other revisions (reference numbering, SI inventory and cross-references) — citation numbering only |
| 65 | The maximum-likelihood estimate of ∆G from the above theorem is the Bennett acceptance ratio.58 For equal number of forward and re | Other revisions (reference numbering, SI inventory and cross-references) |
| 67 | Free energies were estimated from the above equation as implemented in pmx (100 bootstrap resamples, T = 300 K), with the Crooks G | R2.1 / R2.4 |
| 69 | Associations between the computed energies and phenotype were assessed by ordinary least-squares regression against log₁₀-transfor | R3.1 |
| 70 | Differences between genotypes in a structural observable were quantified with the replicate trajectory as the unit of comparison.  | R2.3 |
| 83 | For each system, we ran three independent 100ns all-atom, explicit-solvent production simulations. Since RT undergoes domain rearr | R2.1 / Other revisions (figure cross-reference) |
| 94 | The geometric observables varied across genotypes (see Table 3, Figures 3-5, and detailed results in Supplementary Table 4). We fi | R2.3 |
| 97 | Burial of the chlorocyanophenyl ring, counted as RT heavy-atom contacts within 4.0 Å, is 19.5 ± 0.5 in WT (Table 3). The largest r | R2.3 |
| 98 | Figure 3: Changes in RT–DOR contacts across the genotype panel, referenced to WT. Each cell gives the shift of a structural metric | R1.4 |
| 100 | G190A leaves the DOR pose essentially unperturbed, with chlorocyanophenyl ring burial (20.0 ± 1.2), stacking angle with Tyr188 sid | R2.3 / R1.3 |
| 101 | Figure 4: Trajectory-level changes in DOR contacts and position. (A) Burial of the DOR chlorocyanophenyl ring, defined as the numb | R1.4 / Other revisions (original Figure 3 -> 4) |
| 102 | One of the most dramatic observations is in simulations of V106 substitutions, that displace DOR in opposite directions depending  | R2.3 |
| 103 | We further tracked RT conformational remodeling in both holo and apo configurations, including pocket volume, side-chain rotamers  | R1.2 |
| 105 | Figure 5: Mutation-induced changes in RT's own conformation, as the difference between the bound and unliganded states, referenced | R1.4 |
| 106 | In WT RT simulated here, binding DOR rotates Tyr181 by −18° and Tyr188 by +53°, from unliganded values of −68° and −97° (Supplemen | R1.2 |
| 107 | Protein-level fluctuations and secondary structure are summarized in Supplementary Section 8 and Supplementary Figure 5. In WT, bi | R1.2 |
| 108 | In K103N+M230L, the unliganded primer grip is markedly more mobile than WT, 0.772 ± 0.041 Å against 0.570 ± 0.043 Å, and binding D | R1.2 |
| 111 | We asked whether explicit-solvent MD simulations of RT-DOR complex could provide an in silico framework to complement in vitro sus | Other revisions (reference numbering, SI inventory and cross-references) |
| 114 | Loss of aromatic stacking at Tyr181 and Tyr188 is the classical explanation for NVP and DLV resistance,10,11,78 but in our traject | R3.2 / Other revisions |
| 115 | G190A resists first-generation NNRTIs by introducing a steric bulge into a compact region of the pocket,6 but that region has room | Other revisions (reference numbering, SI inventory and cross-references) |
| 116 | In all four V106A-containing genotypes, DOR slides 1.3–1.5 Å toward Ser105, out of its crystallographic pose, and loses roughly a  | Other revisions (reference numbering, SI inventory and cross-references) |
| 117 | Similarly, for K103N+M230L, our calculations suggest negligible binding penalty (∆∆Gbind = 0.59 ± 0.41 kcal/mol) despite a 36-fold | R1.3 |
| 118 | The apo comparisons suggest that a mutation can alter the conformational response of RT to DOR binding even when contacts in the b | R1.3 |
| 119 | There are three important limitations of this work. First, each system was simulated in triplicate from a single seeded pose, so t | R1.1 |
| 120 | Second, we simulate a binary RT–DOR complex with no template/primer and no incoming dNTP, so a resistance mechanism acting through | R1.1 / R2.2 |
| 121 | The third limitation is in the simulation forcefields. The two energy calculation methods in this study use different ligand charg | R1.1 / R2.5 |
| 122 | More generally, resistance can act through large-scale allosteric change, through the conformational equilibrium of the unliganded | R1.1 |
| 126 | Supplementary Text: comparison of independent RT–DOR structures (Supplementary Figure 1); ligand charge-model sensitivity; perturb | Other revisions (complete Associated Content inventory) |
| 130 | Supplementary Table 4: Structural observables in the bound and unliganded states, with per-replicate values, WT-referenced shifts, | Other revisions (complete Associated Content inventory) |
| 131 | Supplementary Table 5: WT-referenced MM/GBSA energies under the Gasteiger and AM1-BCC ligand charge models, and the per-atom parti | Other revisions (complete Associated Content inventory) |
| 132 | Supplementary Table 6: WT-referenced MM/GBSA energies under two snapshot relaxation protocols, with their agreement and their corr | Other revisions (complete Associated Content inventory) |
| 133 | Supplementary Table 7: Comparison of the doravirine-bound NNIBP between PDB 4NCG and the independently determined structure 7Z2G ( | Other revisions (complete Associated Content inventory) |
| 134 | Supplementary Table 8: Placement of Tyr181 and Tyr188 in RT bound to four NNRTIs (XLSX) | Other revisions (complete Associated Content inventory) |

## Supplementary Text changes

| Paragraph | Revised passage | Response coverage |
|---|---|---|
| 0 | Supporting Information | Other revisions (SI organization and contents) |
| 1 | Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase | Other revisions (SI organization and contents) |
| 2 | Rohit Satija, Kaiming Tao, Robert W. Shafer | Other revisions (SI organization and contents) |
| 3 | Division of Infectious Diseases, School of Medicine, Stanford University, Stanford, CA, USA | Other revisions (SI organization and contents) |
| 4 | Corresponding author: rsatija@stanford.edu | Other revisions (SI organization and contents) |
| 6 | Contents | Other revisions (SI organization and contents) |
| 7 | 1 Is the RT-DOR crystal structure representative of the true RT-DOR structural ensemble? | Other revisions (SI organization and contents) |
| 8 | 2 Sensitivity of the MM/GBSA energies to the ligand charge model | Other revisions (SI organization and contents) |
| 9 | 3 The crystallographic pose is recovered from perturbed starting poses | Other revisions (SI organization and contents) |
| 10 | 4 Convergence of the molecular dynamics ensemble | Other revisions (SI organization and contents) |
| 11 | 5 Sensitivity of the MM/GBSA energies to the snapshot relaxation protocol | Other revisions (SI organization and contents) |
| 12 | 6 Non-equilibrium switching work distributions | Other revisions (SI organization and contents) |
| 13 | 7 Finite-size corrections for the charge-changing alchemical legs | Other revisions (SI organization and contents) |
| 14 | 8 Protein-level conformational analysis | Other revisions (SI organization and contents) |
| 15 | 9 Conformation of the unliganded binding pocket | Other revisions (SI organization and contents) |
| 17 | 1. Is the RT-DOR crystal structure representative of the true RT-DOR structural ensemble? | R1 general |
| 18 | We address this question in two steps: first by comparing the DOR pose in the crystal structure (PDB ID: 4NCG) with a second, inde | R1 general |
| 19 | To calibrate these alignment results, we applied the identical superposition to three further complexes of wild-type RT with an NN | R1 general |
| 21 | Supplementary Figure 1. The DOR-NNIBP pose is reproducible across independently determined RT-DOR structures. (A) Doravirine in th | R1 general (Supplementary Figure 1) |
| 22 | 2. Sensitivity of the MM/GBSA energies to the ligand charge model | R2.5 |
| 23 | DOR parameters for the equilibrium MD simulations were generated with Gasteiger partial charges, as opposed to AM1-BCC for the FEP | R2.5 |
| 24 | 3. The crystallographic pose is recovered from perturbed starting poses | R1.2 |
| 25 | To test whether representative systems from independent starting poses can recover the same DOR-bound state used in the simulation | R1.2 |
| 27 | Supplementary Figure 2: Perturbed doravirine poses against the wild-type sampled basin. Axes are the pocket-frame coordinates afte | R1.2 (Supplementary Figure 2) |
| 28 | 4. Convergence of the molecular dynamics ensemble | R2.1 |
| 29 | Three 100 ns replicates are long relative to local relaxation of the DOR pose but short relative to RT domain rearrangements and t | R2.1 |
| 33 | Supplementary Figure 3. DOR Root Mean Squared Deviation (RMSD) and distance between centers of mass (COM) of DOR and RT indicate 1 | R2.1 |
| 35 | 5. Sensitivity of the MM/GBSA energies to the snapshot relaxation protocol | R3.3 |
| 36 | Each snapshot is energy-minimized before its energy is evaluated. We compared hydrogen-only relaxation, with heavy atoms restraine | R3.3 |
| 37 | 6. Non-equilibrium switching work distributions | R2.1 / Other revisions (existing switching diagnostic) |
| 39 | Supplementary Figure 4. Non-equilibrium switching work distributions for all alchemical transformations. Forward (0→1, blue) and r | R2.1 / Other revisions (renumbered existing legend) |
| 40 | 7. Finite-size corrections for the charge-changing alchemical legs | R2.4 |
| 41 | The leading periodicity (net-charge self-energy) term of the Rocklin scheme is +0.0407 kcal/mol per phase and cancels to below 0.0 | R2.4 |
| 42 | 8. Protein-level conformational analysis | R1.4 |
| 44 | Supplementary Figure 5. Protein-level conformational analysis. Per-residue Cα RMSF after superposition on the NNIBP-lining residue | R1.4 |
| 49 | 9. Conformation of the unliganded binding pocket | R1.2 |
| 50 | Three 100 ns apo replicates were run for WT and all nineteen simulated mutant genotypes, including the F227C alchemical intermedia | R1.2 |
| 52 | Supplementary Figure 6. Apo and holo NNIBP volume distributions. Per-frame pocket volume for the bound state over full production, | R1.2 |

## Figures tables and references

- Figure 1C: mapped mutation positions — R1.4; experimental structure comparison — R1 general and Methods quote.
- Main Figures 3 and 5: new WT-referenced contact and apo/holo summaries — R1.4; original Figure 3 retained as Figure 4 — R1.4 and Other revisions.
- Main Table 1: phenotype values unchanged; citation caches renumbered — Other revisions.
- Main Table 2: all text/data cells unchanged from submission — no new result claimed.
- Main Table 3: Y318F Val106-distance SEM displayed as 0.004 instead of 0.00 — Other revisions; remaining cells unchanged.
- Supplementary Table 4: both-state values and replicate-level uncertainty/pooled tests — R1.2, R1.4, R2.3; Tables 5–6 sensitivity checks — R2.5, R3.3; Tables 7–8 independent structures — R1 general.
- Supplementary Figures 1–2 and 5–6: new analyses — R1 general, R1.2 and R1.4. Supplementary Figure 3 expands submitted convergence plots — R2.1. Supplementary Figure 4 retains the submitted switching-work diagnostic with new numbering — R2.1 and Other revisions.
- New bibliography entries 44–46 (experimental structures) and 49–51 (minimization methods), plus reference/cache renumbering — Methods quotes in R1 general and R3.3, and Other revisions.
- TOC page and separate SI cover page — Other revisions; editorial typography/cross-references — R3.2 and Other revisions.

## Verification

- All 44 original reviewer/comment heading paragraphs preserved verbatim and formatted bold.
- All response prose and location labels black/plain; all manuscript/SI excerpts black/italic.
- 52 quotation-paragraph checks passed against the current manuscript/SI, including existing multi-paragraph excerpts and native equations.
- All 25 existing citation instruction blocks retained; additional cited excerpts copied natively from the manuscript.
- Manuscript, SI, tables, cover page and graphic files left unchanged.

This is a reporting/coverage audit, not a new numerical reanalysis or a claim that the simulations establish convergence of all degrees of freedom.
