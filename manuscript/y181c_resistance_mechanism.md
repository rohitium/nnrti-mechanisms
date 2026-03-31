# Y181C Doravirine Susceptibility Mechanism

## Summary
`Y181C` is a classical NNRTI resistance mutation, but doravirine (`DOR`) usually remains active against it. The structural explanation is that `DOR` does not rely on the `Y181` aromatic interaction that is important for rilpivirine (`RPV`) and several older NNRTIs. Instead, `DOR` places its chlorophenol group against `Y188`.

## Mechanistic Model
1. In the `RPV` complex, `Y181` forms part of the aromatic contact surface for the inhibitor.
2. In the `DOR` complex, the chlorophenol moiety instead stacks with `Y188`, and `Y181` is displaced toward the rim of the pocket.
3. The `Y181C` substitution removes an aromatic side chain at position `181`.
4. Because `DOR` does not depend strongly on `Y181`, the effect on `DOR` susceptibility is limited.

## Evidence Specific To Doravirine
- In serum-containing assays, `DOR` retained activity against `Y181C`, with an `IC50` of `31 nM` and inhibitory quotient `27` ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)).
- In `DOR` resistance-selection experiments at clinically relevant concentrations, `Y181C` did not produce breakthrough virus ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)).
- In the large clinical-isolate analysis, single unique `Y181C` had median fold change `1.6` (`1.2–1.8`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

## Structural Basis
Smith et al. explicitly compare the `RPV` and `DOR` binding modes and show that `RPV` stacks with `Y181`, whereas `DOR` stacks with `Y188`. In the `DOR` complex, `Y181` is shifted toward the rim of the pocket. They state that the lack of interaction with `Y181` explains why `Y181C` remains susceptible to `DOR` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).

The local `DOR`-bound MD trajectories are consistent with that geometry. In the aligned `WT` simulations, the `Y181` aromatic-ring centroid is closer to the pocket entrance centroid than the `Y188` ring centroid in all sampled frames, and `Y181` is farther from the chlorophenol-side ligand centroid than `Y188` in most frames. This is consistent with a rim-facing `Y181` and a more ligand-facing `Y188`.

The `Y181C` trajectories also preserve the main `DOR` binding geometry. In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, `Y181C` remains close to `WT` in ligand pose RMSD (`1.414` vs `1.276 A`), keeps the ligand on the palm side of the pocket (`ligand_palm_distance_angstrom` `2.709` vs `2.991 A`), and preserves the deeper `LYS103` anchor (`residue_min_distance_LYS103_angstrom` `2.907` vs `2.989 A`). The same table shows that `Y181C` does not weaken the `188`-side environment; `residue_min_distance_TYR188_angstrom` is `3.098 A` in `Y181C` versus `3.284 A` in `WT`.

This pattern is consistent with the structural interpretation from Smith et al.: `Y181` is not the dominant aromatic contact for `DOR`, so the `Y181C` substitution can be tolerated as long as the `Y188`-centered and deeper-pocket interactions remain intact.

## Energy Analysis
In `results/tables/holo/mmgbsa_replicate_metrics.csv`, `Y181C` also appears energetically weaker than `WT` by MM/GBSA:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `Y181C`: mean `-111.45`, median `-112.67`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `Y181C`: median `-264.32`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `Y181C`: median `172.70`

Because `Y181C` remains phenotypically susceptible to `DOR`, this is another case in which the MM/GBSA summary alone overstates the practical effect of the mutation. The structural interpretation remains that `DOR` is protected because it does not depend strongly on `Y181`.

## References
- Asante-Appiah E, et al. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- Feng M, Sachs NA, Xu M, Grobler J, Blair W, Hazuda DJ, Miller MD, Lai M-T. Doravirine suppresses common nonnucleoside reverse transcriptase inhibitor-associated mutants at clinically relevant concentrations. *Antimicrob Agents Chemother*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
