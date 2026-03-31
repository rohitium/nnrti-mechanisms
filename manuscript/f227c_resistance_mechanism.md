# F227C Resistance Mechanism

## Summary
`F227C` is mechanistically plausible as a doravirine (`DOR`) resistance mutation because `F227` lies in the upper hydrophobic region of the NNRTI-binding pocket that accommodates the distal portion of `DOR`. Published `DOR` literature places `F227C/L` among the substitutions associated with the greatest reductions in susceptibility, although clinical failure patterns often involve `F227C` in combination with other mutations rather than as the dominant isolated pathway.

## Mechanistic Model
1. In the wild-type RT-`DOR` structure, the upper pocket containing `F227`, `Y188`, `W229`, and `L234` accommodates the distal end of `DOR`.
2. `F227` contributes aromatic bulk to this hydrophobic region.
3. `F227C` replaces phenylalanine with the smaller cysteine side chain.
4. This changes upper-pocket packing and reduces favorable interactions with the distal portion of `DOR`.

## Evidence Specific To Doravirine
- Reviews of `DOR` resistance patterns classify `F227C/L` among the canonical substitutions associated with the greatest reductions in susceptibility ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).
- Clinical-development analyses report substitutions at positions `106` and `227` as the most prevalent emergent `DOR` substitutions in virologic failure, with `F227C` commonly occurring alongside `V106I/A/M` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).
- The strongest published clinical-development examples are combination isolates such as `V106I/F227C` and `A98G/F227C/M184V` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

## Structural Basis
Smith et al. show that the upper hydrophobic pocket containing `F227` is displaced toward `DOR` relative to the `RPV` complex and contributes to the wild-type `DOR` binding environment ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct `F227C-DOR` mutant structure was not identified here.

## Evidence From The Current MD Dataset
The current MD dataset does not support a strong `106-227` corridor-opening effect in `F227C` alone.

- The direct `V106-F227` sidechain-sidechain minimum heavy-atom distance was computed from the raw `F227C` and matched `WT` trajectories for replicates `2` and `3`.
- In replicate `2`, the mean mutant-minus-WT shift was `-0.05 Å`.
- In replicate `3`, the mean mutant-minus-WT shift was `+0.45 Å`.
- Averaged across the two matched replicates, the mean shift was therefore only about `+0.20 Å`.

This is much smaller than the corresponding shifts in the double mutants `V106I+F227C` and `V106A+F227L`, where the `106-227` sidechain distance increases by approximately `+0.87 Å` and `+1.33 Å`, respectively, relative to the matched wild-type panels. In the current dataset, `F227C` alone therefore does not reproduce the strong corridor-opening signal seen in the `106/227` combination mutants.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` shows weaker binding than `WT`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `F227C`: mean `-118.07`, median `-121.36`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `F227C`: median `-265.11`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `F227C`: median `176.05`

This energetic pattern is directionally consistent with upper-pocket packing loss at residue `227`, even though the single-mutant structural signal is weaker than in the `106/227` combinations.

## References
- de Béthune M-P, et al. Pharmaceutical, clinical, and resistance information on doravirine, a novel non-nucleoside reverse transcriptase inhibitor for the treatment of HIV-1 infection. 2020 review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
