# Notes removed from the supplementary workbooks

The Notes tabs were deleted on 01-Nov-2026: the columns now carry their own
units, so definitions are no longer needed. These lines are what was deleted.
Most were column definitions and are now redundant. The ones marked KEEP are
method or limitation statements that no column name can express -- they have
nowhere to live except the Supplementary Text.

## Supplementary Table 5

    Supplementary Table 5. Sensitivity of MM/GBSA energies to the snapshot relaxation protocol.
    Sheet 'Table 5': WT-referenced ddE in kcal/mol for all 19 mutant genotypes, mean +/- SEM over three replicates, under two protocols.
      h_relax    hydrogens relaxed, heavy atoms restrained (k = 10000 kJ/mol/nm^2)
      full_min   all atoms free; the protocol used for Table 2 of the main text
    Sheet 'Table 5b agreement': Spearman correlation between the two protocols across the 19 genotypes, per energy component.
    Sheet 'Table 5b phenotype': each protocol against log10 fold reduction in susceptibility, over the 13 genotypes with an established phenotype.
    Sampling is identical across protocols and matches the main-text run: 100 snapshots spread evenly over the post-equilibration region of each 100 ns replicate, 0.25 discard fraction, 2000 iteration cap, CPU platform, double precision.
    Single-point evaluation, with no relaxation, was tested on six genotypes and is described in the Supplementary Text. It is not tabulated here: its mean is set by a small number of snapshots carrying steric contacts rather than by the protocol, so the values are not comparable with the two columns above.

## Supplementary Table 6

    Supplementary Table 6. Sensitivity of the MM/GBSA energies to the ligand partial charge model.
    Sheet 'Table 6': WT-referenced ddE in kcal/mol, mean +/- SEM over three replicates, for each energy component under both charge models.
      gast     Gasteiger charges, as used throughout the manuscript
      am1bcc   AM1-BCC charges assigned by openff-toolkit
    The same trajectories and the same snapshots are scored under both models, so the difference isolates the charge assignment.
    Sheet 'DOR partial charges': per-atom charges under both models. AM1-BCC is the more polarised assignment (mean |q| 0.228 e against 0.138 e; RMS difference 0.165 e, maximum 0.390 e), with the largest differences on the triazolinone and pyridinone carbonyl and amide atoms.
    Sheet 'Per replicate AM1-BCC': the individual replicate values behind the means, alongside the wild-type reference.
    Genotypes re-scored under AM1-BCC so far: V106A, Y188L, each against the wild-type reference. Re-scoring is a per-replicate MM/GBSA calculation on the existing trajectories and does not require new molecular dynamics.

## Supplementary Table 7

    Supplementary Table 7. Unliganded (apo) and bound (holo) NNIBP volume.
    Volumes are the convex-hull proxy over the 16 NNIBP-lining Calpha atoms, in cubic angstroms.
    Sheet 'Table 7': per genotype, the apo volume over the matched window beside the bound-state volume of Table 3 of the main text.
      apo_volume           mean over three replicates; sd is across frames
      holo_volume          bound state, over full production
      induced_fit_opening  holo minus apo: the opening DOR has to pay for
      chi1_Y181, chi1_Y188 mean side-chain chi1 in the apo ensemble, degrees
      dVolume_vs_WT        apo volume relative to wild-type apo
    The wild-type apo trajectories are 10 ns and set the comparison window, so every apo trajectory is cut to a matched 10 ns window. The bound-state volumes come from full 100 ns production, so the two columns are not matched in sampling; the apo column is the shorter of the two.
    Sheet 'Crystal reference': the same proxy evaluated on crystal structures. The unliganded crystals 1DLO and 1HMV bracket the simulated apo ensemble, and the DOR-bound 4NCG sits well above both.
    The apo pockets reported here were equilibrated without relocating bulk water into the cavity, which leaves the cavity under-solvated and biases the volume low. The induced-fit openings are therefore upper bounds.

## Supplementary Table 9

    Supplementary Table 9. Placement of Tyr181 and Tyr188 in RT bound to four NNRTIs.
    Sheet 'Table 9': one row per deposited structure. All structures were superposed on the
    sixteen NNIBP-lining Calpha atoms of 4NCG (15 in p66, Glu138 in p51), so side-chain
    differences are reported inside a common pocket frame and are not whole-enzyme differences.
      pocket_ca_rmsd_vs_ref_A   Calpha RMSD of the 16 pocket atoms after that fit
      chi1_deg, chi2_deg        side-chain dihedrals, N-CA-CB-CG and CA-CB-CG-CD1
      dchi1_vs_ref_deg          signed circular difference in chi1 from 4NCG
      ring_shift_vs_ref_A       displacement of the six-membered ring centroid from 4NCG
      min_dist_to_ligand_A      minimum heavy-atom distance to that structure's OWN bound
                                inhibitor; no superposition enters this column
      nnibp_volume_A3           pocket volume by the routine used on the trajectories
    Tyr181 occupies a different rotamer in the two doravirine structures than in the
    rilpivirine, efavirenz and nevirapine complexes, by 90-95 degrees in chi1 and 4.2-4.9 A at
    the ring centroid. Tyr188 is in the same rotamer throughout, differing by 13-20 degrees and
    1.3-2.6 A. The two independent doravirine structures agree with each other to 10 degrees and
    0.5 A, so the placement is a property of the drug rather than of one crystal form.
    Sheet 'Atoms': coordinates of the backbone, ring and hydroxyl atoms of both tyrosines after
    the pocket superposition, for anyone wishing to reproduce the overlay.
