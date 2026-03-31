# K103N Doravirine Susceptibility Mechanism

## Summary
`K103N` is a major resistance mutation for older NNRTIs, but doravirine (`DOR`) usually retains activity against it. In published phenotyping studies, `K103N` alone causes little or no reduction in `DOR` susceptibility. The structural explanation is that `DOR` does not depend strongly on the `K103` side chain in the same way as efavirenz or nevirapine.

## Mechanistic Model
The proposed mechanism can be written as a short causal chain:

1. In the wild-type RT-`DOR` complex, the triazolone region of `DOR` contacts the main-chain atoms near `K103`.
2. `DOR` binds deeper in the NNRTI pocket than rilpivirine and uses a contact network centered more strongly on `Y188`, `V106`, `L100`, and `P236`.
3. Substituting lysine with asparagine at position `103` changes the side chain but leaves the local backbone geometry available.
4. The principal `DOR` contacts are therefore preserved well enough for antiviral activity to remain high.

## Evidence Specific To Doravirine
- In serum-containing antiviral assays, `DOR` remained potent against `K103N`, with an `IC50` of `21 nM` and inhibitory quotient `39` ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)).
- In a large clinical-isolate analysis, single unique `K103N` had median `DOR` fold change `1.0` (`0.7–1.3`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).
- Clinical-development reviews report maintained virologic suppression in participants with baseline `K103N` and no treatment-emergent `K103N` substitutions during `DOR` trials ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

## Structural Basis
The wild-type RT-`DOR` structure shows that `DOR` interacts with the main-chain atoms of `K103` rather than depending primarily on the `K103` side chain. The same structural analysis shows that `DOR` binds deeper in the pocket than `RPV` and uses a distinct interaction pattern centered on `Y188`, `V106`, `L100`, and `P236` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).

## Energy Analysis
In `results/tables/holo/mmgbsa_replicate_metrics.csv`, `K103N` shows weaker MM/GBSA binding than `WT` overall, but the component pattern is mixed:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `K103N`: mean `-116.21`, median `-119.15`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `K103N`: median `-39.63`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `K103N`: median `189.30`

So `K103N` shows more favorable direct electrostatics but a large opposing GB penalty. Because `K103N` remains phenotypically susceptible, the MM/GBSA total appears to overestimate the functional impact of this substitution in isolation.

## References
- Asante-Appiah E, et al. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- Feng M, Sachs NA, Xu M, Grobler J, Blair W, Hazuda DJ, Miller MD, Lai M-T. Doravirine suppresses common nonnucleoside reverse transcriptase inhibitor-associated mutants at clinically relevant concentrations. *Antimicrob Agents Chemother*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
