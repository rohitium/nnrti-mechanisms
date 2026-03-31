# V106A+L234I Resistance Mechanism

## Summary
`V106A/L234I` is a major doravirine (`DOR`) resistance pathway selected after the initial `V106A` step in in vitro resistance experiments. Published summaries report very large reductions in `DOR` susceptibility for this combination. The mechanism is a combined disruption of the `V106` core contact and the upper-pocket contour that accommodates the distal portion of `DOR`.

## Mechanistic Model
1. `V106A` removes a direct contact with the `DOR` pyridone core.
2. `L234I` changes the upper-pocket region near the distal end of `DOR`.
3. Together, these substitutions alter both the central anchoring interaction and the upper-pocket shape.
4. `DOR` binding is therefore strongly weakened.

## Evidence Specific To Doravirine
- In `DOR` resistance-selection experiments, `L234I` emerged after `V106A` as the alternative major pathway to `F227L` ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)).
- Review summaries report `>150-fold` reduced susceptibility for `V106A/L234I` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).

## Structural Basis
`L234` lies in the upper portion of the NNRTI pocket that is shifted toward `DOR` in the wild-type `DOR` complex and contributes to the environment around the distal end of the inhibitor ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct `V106A/L234I-DOR` co-structure was not identified here.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` shows a strong energetic penalty for `V106A/L234I`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106A+L234I`: mean `-139.34`, median `-146.80`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `V106A+L234I`: median `-264.32`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `V106A+L234I`: median `-19.07`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `V106A+L234I`: median `134.37`

The most prominent penalty is the loss of favorable van der Waals packing, which fits a primer-grip / upper-pocket reshaping mechanism.

## References
- de Béthune M-P, et al. Pharmaceutical, clinical, and resistance information on doravirine, a novel non-nucleoside reverse transcriptase inhibitor for the treatment of HIV-1 infection. 2020 review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)
- Feng M, Wang D, Grobler JA, Hazuda DJ, Miller MD, Lai M-T. In vitro resistance selection with doravirine (MK-1439), a novel nonnucleoside reverse transcriptase inhibitor with distinct mutation development pathways. *Antimicrob Agents Chemother*. 2015. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
