# V106A+F227L Resistance Mechanism

## Summary
`V106A/F227L` is a high-confidence doravirine (`DOR`) resistance pathway. It emerged in `DOR` resistance-selection experiments after the initial `V106A` step and is associated with very large reductions in `DOR` susceptibility. The current MD dataset supports a mechanism in which the `106-227` corridor opens strongly in the double mutant, consistent with upper-pocket repacking around the distal portion of `DOR`.

## Mechanistic Model
1. `V106A` removes a direct hydrophobic contact with the `DOR` pyridone core.
2. `F227L` alters the upper hydrophobic portion of the pocket that accommodates the distal end of `DOR`.
3. In the double mutant, the sidechain separation between positions `106` and `227` increases markedly relative to the matched wild-type panel.
4. This is consistent with coordinated opening and repacking of the upper `106-227` corridor used by `DOR`.

## Evidence Specific To Doravirine
- In `DOR` resistance-selection experiments, `V106A` was followed by `F227L` as `DOR` concentration increased ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)).
- Review summaries report `>150-fold` reduced susceptibility for `V106A/F227L` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).

## Structural Basis
The wild-type `DOR` structure shows direct contact with `V106` and involvement of `F227` in the upper hydrophobic region that surrounds the distal end of the inhibitor ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct `V106A/F227L-DOR` co-structure was not identified here.

## Evidence From The Current MD Dataset
The current trajectory-derived data support a strong corridor-opening signal for `V106A+F227L`:

- In `results/tables/holo/drm_sidechain_distance_timeseries_all_mutations.csv`, the sidechain-sidechain distance between positions `106` and `227` (`c1_to_c2`) is larger in the mutant than in the matched wild-type panel in all three replicates:
  - rep 1: `+1.15 Å`
  - rep 2: `+1.43 Å`
  - rep 3: `+1.42 Å`
- In the same file, the double mutant spends much more time in an open `c1_to_c2` state than the wild-type panel:
  - mutant `c1_to_c2 > 4.0 Å`: `95.5%`
  - wild-type `c1_to_c2 > 4.0 Å`: `28.1%`
  - mutant `c1_to_c2 > 5.0 Å`: `57.7%`
  - wild-type `c1_to_c2 > 5.0 Å`: `0.2%`
- The direct sidechain-to-ligand distances are less uniform than the corridor-opening coordinate. The `V106A` sidechain-to-`DOR` distance (`c1_to_dor`) is modestly larger than wild type, while the residue `227` sidechain-to-`DOR` distance (`c2_to_dor`) changes less consistently across replicates.

Taken together, these data support an opening of the `106-227` corridor as the clearest structural signature in the current MD dataset for `V106A+F227L`.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` is consistent with substantial weakening of the bound state:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106A+F227L`: mean `-121.85`, median `-123.36`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `V106A+F227L`: median `-267.88`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `V106A+F227L`: median `-17.57`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `V106A+F227L`: median `162.06`

This pattern supports a combined packing and polar-anchor failure mode, consistent with the strong `106-227` corridor opening seen in the trajectories.

## References
- de Béthune M-P, et al. Pharmaceutical, clinical, and resistance information on doravirine, a novel non-nucleoside reverse transcriptase inhibitor for the treatment of HIV-1 infection. 2020 review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)
- Feng M, Wang D, Grobler JA, Hazuda DJ, Miller MD, Lai M-T. In vitro resistance selection with doravirine (MK-1439), a novel nonnucleoside reverse transcriptase inhibitor with distinct mutation development pathways. *Antimicrob Agents Chemother*. 2015. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
