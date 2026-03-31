# Y318F Resistance Mechanism

## Summary
`Y318F` is one of the few single substitutions that crosses the doravirine (`DOR`) biological cutoff in published susceptibility datasets. The phenotypic effect is supported, but the structural mechanism is less direct than for `V106A` or `Y188L`, because a mutant `Y318F-DOR` structure was not identified here.

## Mechanistic Model
1. `Y318` is a distal NNRTI-binding-pocket residue long associated with NNRTI resistance.
2. `DOR` binds deeper in the pocket than `RPV`, but the distal pocket still contributes to ligand accommodation.
3. Replacing tyrosine with phenylalanine removes the phenolic hydroxyl and changes the local chemical environment without removing aromatic character.
4. That change is sufficient in published phenotyping studies to reduce `DOR` susceptibility into the resistant range.

## Evidence Specific To Doravirine
- In the large clinical-isolate analysis, single unique `Y318F` had median fold change `11.0` (`3.0–14.1`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).
- In clinical-development follow-up, site-directed `Y318F` conferred a `9.9-fold` reduction in `DOR` susceptibility ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

## Structural Basis
`Y318` has long been linked to NNRTI resistance more broadly ([Harrigan et al., 2002](https://pmc.ncbi.nlm.nih.gov/articles/PMC136283/)). In the `RPV` complex, `Y318` stacks with the benzonitrile group, whereas the published `DOR` structure emphasizes deeper interactions with `Y188`, `V106`, `L100`, and `P236` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct structural explanation for the `Y318F` effect on `DOR` was not identified here.

## Energy Analysis
The MM/GBSA signal for `Y318F` is modest compared with the stronger resistance mutations:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `Y318F`: mean `-147.34`, median `-150.00`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `Y318F`: median `-270.01`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `Y318F`: median `-29.44`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `Y318F`: median `148.30`

This near-WT energetic profile fits the idea that `Y318F` acts through a comparatively subtle distal-coupling mechanism rather than through gross loss of local pocket packing.

## References
- Asante-Appiah E, et al. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- Harrigan PR, Salim M, Stammers DK, Wynhoven B, Brumme ZL, McKenna P, Larder B, Kemp SD. A mutation in the 3' region of the human immunodeficiency virus type 1 reverse transcriptase (Y318F) associated with nonnucleoside reverse transcriptase inhibitor resistance. *J Virol*. 2002. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC136283/)
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
