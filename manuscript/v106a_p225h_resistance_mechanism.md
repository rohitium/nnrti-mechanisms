# V106A+P225H Resistance Mechanism

## Summary
`V106A/P225H` is supported as a doravirine (`DOR`) resistance pathway in susceptibility summaries and clinical-development follow-up. The exact structural mechanism is less directly established than for `V106A/F227L` or `V106A/L234I`, but the available structural framework places `P225` in the same upper-pocket region involved in `DOR` binding.

## Mechanistic Model
1. `V106A` removes a direct contact between `DOR` and the `V106` side chain.
2. `P225H` alters a residue in the upper portion of the pocket near the distal `DOR` substituents.
3. The combination perturbs both the central contact region and a neighboring upper-pocket region used by `DOR`.
4. This produces a much larger reduction in susceptibility than expected from `V106A` alone.

## Evidence Specific To Doravirine
- In the transmitted-resistance review, `V106A/P225H` is listed with `>64-fold` reduced `DOR` susceptibility ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).
- Clinical-development follow-up described a resistant isolate containing `V106A/P225H/Y318F/K65R` with `>210-fold` reduced susceptibility to `DOR` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

## Structural Basis
The published wild-type `DOR` structure shows that the upper pocket containing `P225`, `F227`, and `L234` is shifted relative to the `RPV` complex and contributes to `DOR` accommodation ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct `V106A/P225H-DOR` mutant structure was not identified here.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` shows that `V106A/P225H` does not have a large net total-energy penalty despite its very high phenotype:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106A+P225H`: mean `-144.63`, median `-144.73`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `V106A+P225H`: median `-20.04`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `V106A+P225H`: median `145.45`

This is consistent with a mechanism dominated by pose rearrangement or gating rather than by one large static packing penalty.

## References
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
