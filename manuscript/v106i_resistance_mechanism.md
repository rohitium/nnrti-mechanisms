# V106I Doravirine Susceptibility Mechanism

## Summary
`V106I` is supported as a polymorphism rather than a canonical doravirine (`DOR`) resistance mutation. By itself, `V106I` usually causes little or no loss of `DOR` susceptibility. The structural explanation is that `V106I` preserves a branched hydrophobic side chain at a residue that contacts the `DOR` pyridone core, and the current MD dataset indicates that it also preserves the local `105-106` pocket-wall geometry better than `V106A`.

## Mechanistic Model
1. In the wild-type RT-`DOR` structure, `V106` contacts the `DOR` pyridone core.
2. Substituting valine with isoleucine preserves a hydrophobic branched side chain at this position.
3. The current MD dataset indicates that the local `105-106` wall is perturbed less than in `V106A`, because `DOR` remains close to residue `106` and does not shift toward `SER105`.
4. `V106I` alone usually remains phenotypically susceptible, but it can participate in resistant combination pathways.

## Evidence Specific To Doravirine
- In the 2024 JID study, site-directed `V106I` remained below the `DOR` biological cutoff, unlike `V106A`, `V106M`, and `Y188L` ([Giammarino et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38206187/)).
- In clinically derived `V106I` viruses, median fold change was `1.2` in subtype B and `1.8` in non-B viruses ([Giammarino et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38206187/)).
- In the large clinical-isolate analysis, single unique `V106I` had median fold change `0.8` (`0.6–1.3`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).
- Martin et al. conclude that the totality of available data supports `V106I` as a polymorphism rather than a `DOR` resistance-associated substitution ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

## Structural Basis
Smith et al. describe `V106` as a residue that interacts with the `DOR` pyridone core ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). That same contact is consistent with the different phenotypic behavior of `V106I` and `V106A`: isoleucine preserves a bulky hydrophobic side chain, whereas alanine removes it.

## Evidence From The Current MD Dataset
The local `V106I-DOR` trajectories are consistent with preservation of the `105-106` wall geometry:

- Direct residue-`106` sidechain-to-`DOR` distances, measured from the trajectories, remain close to or slightly tighter than wild type:
  - `WT`: mean `3.36 Å`
  - `V106I`: mean `3.22 Å`
- `SER105` does not move closer to compensate:
  - `SER105` sidechain-to-`DOR`
  - `WT`: mean `6.65 Å`
  - `V106I`: mean `7.15 Å`
  - `SER105 OG`-to-`DOR`
  - `WT`: mean `7.44 Å`
  - `V106I`: mean `7.95 Å`
- In the precomputed feature table `results/analysis/ligand_pocket_features/tables/frame_features.csv`, `V106I` also remains close to the wild-type binding geometry:
  - `ligand_pose_rmsd_angstrom`
  - `WT`: `1.276 Å`
  - `V106I`: `1.345 Å`
  - `residue_min_distance_VAL106_angstrom`
  - `WT`: `3.389 Å`
  - `V106I`: `3.228 Å`

These data support a model in which `V106I` preserves the local packing role of residue `106` and does not induce the compensatory shift toward `SER105` that is seen in `V106A`.

The same local-packing interpretation is consistent with the lysine-side distances in `results/analysis/ligand_pocket_features/tables/frame_features.csv`:

- `residue_min_distance_LYS101_angstrom`
  - `WT`: `3.513 Å`
  - `V106I`: `3.504 Å`
  - `V106A`: `3.888 Å`

So `V106I` remains close to wild type at the `LYS101` side of the local wall, whereas `V106A` shifts away from it. By contrast, the deeper `LYS103` polar contacts are not lost in `V106A`, which indicates that the main effect is local wall repacking around residues `101/105/106` rather than complete disruption of the deeper `DOR` anchoring geometry.

The MM/GBSA analysis is consistent with the same interpretation. In `results/tables/holo/mmgbsa_replicate_metrics.csv`, `V106I` does not show a clean, stable energetic penalty comparable to `V106A`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106I`: mean `-127.88`, median `-141.15`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `V106I`: median `-265.86`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `V106I`: median `-23.62`

The `V106I` averages are strongly influenced by one outlying replicate, and the replicate-to-replicate variance is large. Taken together with the structural measurements, the energy analysis is more consistent with a mostly preserved binding mode than with a reproducible resistance-driving energetic defect.

## References
- Asante-Appiah E, et al. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- Giammarino F, et al. Prevalence and Phenotypic Susceptibility to Doravirine of the HIV-1 Reverse Transcriptase V106I Polymorphism in B and Non-B Subtypes. *J Infect Dis*. 2024. [PubMed](https://pubmed.ncbi.nlm.nih.gov/38206187/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
