# A98G+F227C Resistance Mechanism

## Summary
`A98G/F227C` is supported as a doravirine (`DOR`) resistance combination by clinical-development data from closely related isolates. The structural mechanism is less directly resolved than for `V106A`-based pathways, but published evidence supports a model in which `F227C` provides a canonical upper-pocket resistance component and `A98G` contributes an additional local conformational change.

## Mechanistic Model
1. `F227C` alters the upper hydrophobic region of the NNRTI pocket used by the distal end of `DOR`.
2. `A98G` by itself is not a canonical `DOR` resistance substitution.
3. In combination with a canonical resistance substitution such as `F227C`, `A98G` is associated with substantially reduced `DOR` susceptibility.
4. The likely effect is cooperative destabilization of the binding-pocket geometry that accommodates `DOR`.

## Evidence Specific To Doravirine
- Clinical-development data include an `A98G/F227C/M184V` isolate with `>93-fold` reduced susceptibility to `DOR` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).
- Reviews of transmitted resistance note that `A98G` can contribute to strong `DOR` resistance when present with canonical `DOR` resistance substitutions such as `F227C`, even though `A98G` alone is not usually associated with large `DOR` resistance ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).

## Structural Basis
The residue-specific structural mechanism for `A98G` in the `DOR` complex is not directly defined in the published wild-type structure. The strongest structural anchor in this combination is `F227`, which lies in the upper hydrophobic portion of the pocket that accommodates the distal end of `DOR` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct `A98G/F227C-DOR` co-structure was not identified here.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` is consistent with strong destabilization:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `A98G+F227C`: mean `-115.08`, median `-117.58`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `A98G+F227C`: median `-262.16`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `A98G+F227C`: median `172.09`

This pattern supports a mechanism in which the canonical `F227C` upper-pocket defect is amplified by an additional conformational contribution from `A98G`.

## References
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
