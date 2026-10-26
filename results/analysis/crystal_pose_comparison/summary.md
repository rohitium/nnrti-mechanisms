# 4NCG vs 7Z2G: is the DOR pose structure-specific?

4NCG: X-ray 2.58 A, RT + DOR, no nucleic acid.
7Z2G: cryo-EM 3.65 A, RT + DOR + 38-mer DNA aptamer, independent group.
Ligand is chemical component 2KW in both; 29 heavy atoms matched by
name, up to the 6 automorphisms of its bond graph (the CF3 fluorines);
3 atom labels were permuted to resolve those.
923 residues resolved in both (527 p66, 396 p51).

## Enzyme

| frame   |   n_fit_residues |   fit_ca_rmsd_A |   p66_ca_rmsd_A |   p51_ca_rmsd_A |   nnibp_ca_rmsd_A |
|:--------|-----------------:|----------------:|----------------:|----------------:|------------------:|
| whole   |              923 |           1.483 |           1.458 |           1.516 |             1.178 |
| palm    |              293 |           1.57  |           1.389 |           2.234 |             0.954 |
| nnibp   |               16 |           0.49  |           1.6   |           2.549 |             0.49  |

Where the two differ, with the palm as the reference frame:

| subdomain                      |   n_residues |   ca_rmsd_palm_frame_A |
|:-------------------------------|-------------:|-----------------------:|
| p66 fingers                    |          103 |                  1.833 |
| p66 palm                       |          112 |                  1.208 |
| p66 thumb                      |           81 |                  1.622 |
| p66 connection                 |          108 |                  0.766 |
| p66 RNase H                    |          123 |                  1.369 |
| p51                            |          396 |                  2.234 |
| p66 NNIBP lining (15 residues) |           15 |                  0.975 |

## DOR

| frame                                         |   n_ligand_atoms |   ligand_rmsd_A |   ligand_centroid_shift_A |
|:----------------------------------------------|-----------------:|----------------:|--------------------------:|
| whole                                         |               29 |           1.256 |                     1.08  |
| palm                                          |               29 |           1.02  |                     0.856 |
| nnibp                                         |               29 |           0.619 |                     0.299 |
| ligand-on-ligand (internal conformation only) |               29 |           0.44  |                     0     |

Ligand heavy-atom RMSD in the pocket's own frame: 0.62 A over 29 atoms.
Internal conformation alone (ligand fitted on ligand): 0.44 A.

### By moiety, in the NNIBP frame

| moiety            |   n_atoms |   rmsd_nnibp_frame_A |   centroid_shift_A |
|:------------------|----------:|---------------------:|-------------------:|
| chlorocyanophenyl |         8 |                0.776 |              0.609 |
| pyridinone        |        11 |                0.472 |              0.361 |
| triazolinone      |         8 |                0.542 |              0.42  |
| linker            |         2 |                0.866 |              0.58  |

### Rotatable dihedrals

| torsion   | atoms           |   4NCG_deg |   7Z2G_deg |   delta_deg |   WT_MD_circmean_deg |   WT_MD_circsd_deg |   4NCG_z_in_WT_MD |   7Z2G_z_in_WT_MD |
|:----------|:----------------|-----------:|-----------:|------------:|---------------------:|-------------------:|------------------:|------------------:|
| tau1      | C12-C7-O-C2     |      -56.2 |      -60.2 |        -4   |                -62.3 |               10   |              0.61 |              0.21 |
| tau2      | C7-O-C2-C1      |      -58.2 |      -69.6 |       -11.4 |                -43.1 |               12.3 |             -1.23 |             -2.16 |
| tau3      | C10-N11-C17-C18 |       87.4 |       77   |       -10.4 |                 78.1 |                7.9 |              1.18 |             -0.13 |
| tau4      | N11-C17-C18-N22 |       95.5 |      110.3 |        14.8 |                 94.2 |               10.6 |              0.12 |              1.51 |

## NNIBP residues, in the NNIBP frame

| subunit   |   resid | resname   |   ca_deviation_A |   all_heavy_rmsd_A |   sidechain_rmsd_A | terminal_atoms_swapped   |   n_atoms_compared |   chi1_4NCG_deg |   chi1_7Z2G_deg |   delta_chi1_deg |   delta_chi2_deg | chi2_two_fold   |   min_dist_to_DOR_4NCG_A |   min_dist_to_DOR_7Z2G_A |   delta_min_dist_A |
|:----------|--------:|:----------|-----------------:|-------------------:|-------------------:|:-------------------------|-------------------:|----------------:|----------------:|-----------------:|-----------------:|:----------------|-------------------------:|-------------------------:|-------------------:|
| p66       |     100 | LEU       |            0.242 |              0.394 |              0.415 | False                    |                  8 |           -54.1 |           -59.4 |             -5.3 |             16   | False           |                     3.83 |                     3.97 |               0.14 |
| p66       |     101 | LYS       |            0.426 |              1.793 |              2.336 | False                    |                  9 |           -77.9 |           -59.7 |             18.2 |            -50.2 | False           |                     3.27 |                     3.47 |               0.2  |
| p66       |     103 | LYS       |            0.215 |              0.427 |              0.505 | False                    |                  9 |           -76.5 |           -68.8 |              7.7 |              0.2 | False           |                     2.79 |                     2.38 |              -0.41 |
| p66       |     106 | VAL       |            0.376 |              0.794 |              0.977 | False                    |                  7 |           158   |          -166.9 |             35.1 |            nan   | False           |                     3.45 |                     3.31 |              -0.14 |
| p66       |     107 | THR       |            0.596 |              0.635 |              0.802 | False                    |                  7 |           -59.1 |           -58.1 |              1   |            nan   | False           |                     5.05 |                     5.41 |               0.36 |
| p66       |     108 | VAL       |            0.482 |              0.713 |              0.807 | False                    |                  7 |           166.7 |           163.4 |             -3.3 |            nan   | False           |                     3.83 |                     3.87 |               0.04 |
| p66       |     179 | VAL       |            0.466 |              0.665 |              0.881 | False                    |                  7 |          -167.6 |           165.2 |            -27.2 |            nan   | False           |                     3.3  |                     3.14 |              -0.16 |
| p66       |     181 | TYR       |            0.273 |              0.573 |              0.626 | False                    |                 12 |           -96   |           -86.4 |              9.6 |            -26.1 | True            |                     3.12 |                     3.17 |               0.05 |
| p66       |     188 | TYR       |            0.242 |              0.271 |              0.284 | False                    |                 12 |           -49.4 |           -46.1 |              3.3 |              2.8 | True            |                     3.15 |                     2.81 |              -0.34 |
| p66       |     189 | VAL       |            0.22  |              0.515 |              0.695 | False                    |                  7 |          -146.3 |           175.8 |            -37.9 |            nan   | False           |                     3.1  |                     3    |              -0.1  |
| p66       |     190 | GLY       |            0.221 |              0.368 |            nan     | False                    |                  4 |           nan   |           nan   |            nan   |            nan   | False           |                     3.05 |                     2.96 |              -0.09 |
| p66       |     227 | PHE       |            1.08  |              1.419 |              1.601 | True                     |                 11 |           -60.5 |           -49.7 |             10.8 |            -41.8 | True            |                     3.5  |                     3.75 |               0.25 |
| p66       |     229 | TRP       |            0.77  |              1.2   |              1.145 | False                    |                 14 |          -178.7 |          -161.1 |             17.6 |            -33.4 | False           |                     3.27 |                     3.15 |              -0.12 |
| p66       |     234 | LEU       |            0.304 |              0.805 |              1.086 | False                    |                  8 |           -54.9 |           -77   |            -22.1 |            -61.4 | False           |                     3.48 |                     3.28 |              -0.2  |
| p66       |     318 | TYR       |            0.327 |              0.355 |              0.349 | False                    |                 12 |           -66.8 |           -70.2 |             -3.3 |             -7.1 | True            |                     3.49 |                     3.34 |              -0.15 |
| p51       |     138 | GLU       |            0.649 |              0.731 |              0.742 | False                    |                  9 |           175.8 |           164.8 |            -11   |             20.8 | True            |                     6.64 |                     6.29 |              -0.35 |

Side-chain RMSD: median 0.80 A, max 2.34 A (LYS101).
|delta chi1| > 30 deg: 2 of 15 residues.

## Superposition-free checks

Ligand-to-residue minimum heavy-atom distances agree to a median of 0.15 A (max 0.41 A); Pearson r = 0.974.

| gate      |   4NCG_min_heavy_dist_A |   7Z2G_min_heavy_dist_A |   delta_A |
|:----------|------------------------:|------------------------:|----------:|
| K101-Y188 |                   10.24 |                    8.89 |     -1.35 |
| V106-Y181 |                    8.28 |                    8.87 |      0.59 |
| L100-F227 |                    7.9  |                    8.84 |      0.94 |
| K103-P225 |                    7.83 |                    8.57 |      0.74 |

| structure   | method         | nucleic_acid   |   pocket_volume_A3 |
|:------------|:---------------|:---------------|-------------------:|
| 4NCG        | X-ray 2.58 A   | no             |              198.3 |
| 7Z2G        | cryo-EM 3.65 A | 38-mer DNA     |              194.1 |

## Caveats

7Z2G is 3.65 A, so its side-chain positions carry more uncertainty than
4NCG's. Within 15 A of the NNIBP Ca centroid, 7Z2G leaves
0 residues unmodelled
and truncates the side chain of 0 more.
Those residues are excluded atom-by-atom from the per-residue comparison; the
pocket-volume figures are not corrected for them, so 7Z2G's volume is
an upper bound.
