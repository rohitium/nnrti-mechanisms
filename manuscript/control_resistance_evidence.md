# Structural Mechanisms For Control Mutations In Doravirine Susceptibility

## Summary
Doravirine (`DOR`) binds the HIV-1 RT NNRTI-binding pocket (`NNIBP`) differently from rilpivirine (`RPV`). Structural analyses place `DOR` deeper in the pocket, with the chlorophenol moiety stacked against `Y188`, the pyridone core contacting `V106` and `L100`, and the triazolone group contacting `P236` and the main-chain atoms of `K103` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). In the same structural comparison, `Y181` is shifted toward the rim of the pocket in the `DOR` complex, whereas residues in the upper pocket, including `F227` and `L234`, are positioned closer to the bound inhibitor than in the `RPV` complex ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).

This binding mode explains the broad pattern seen across the control mutations:

- `K103N`, `Y181C`, and `G190A` are common NNRTI resistance mutations for older drugs, but `DOR` usually retains activity because these substitutions do not directly disrupt the main set of contacts that stabilize `DOR`.
- `V106A` and `Y188L` directly perturb key `DOR` contacts and therefore confer resistance as single mutations.
- `V106`-initiated combination pathways with `F227`, `L234`, or `P225` cause larger resistance because they alter multiple residues in the upper portion of the pocket that is used heavily by `DOR`.
- `V106I` behaves differently from `V106A`: it is observed as a polymorphism in untreated people, has little effect on `DOR` susceptibility by itself, and only becomes strongly associated with resistance in combination backgrounds.

## High-Confidence Negative Controls

### K103N
`K103N` is a major resistance mutation for efavirenz and nevirapine, but `DOR` was designed to retain activity against it. In serum-containing assays, `DOR` remained potent against `K103N`, with an `IC50` of `21 nM` and an inhibitory quotient of `39` ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)). In a large clinical-isolate analysis, single unique `K103N` had median `DOR` fold change `1.0` (`0.7–1.3`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

The structural explanation is that `DOR` does not depend on the same interaction network that makes older NNRTIs vulnerable to `K103N`. In the `DOR` complex, the triazolone group contacts the main-chain atoms of `K103`, and the drug does not rely on the `E138-K101` electrostatic arrangement that is important in the `RPV` complex ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). That geometry reduces the impact of changing the `K103` side chain from lysine to asparagine.

### Y181C
`Y181C` is another classical NNRTI resistance mutation, but `DOR` usually remains active against it. In serum-containing assays, `DOR` retained potency against `Y181C`, with an `IC50` of `31 nM` and inhibitory quotient `27`, and no viral breakthrough was observed in `DOR` selection experiments at clinically relevant concentrations ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)). In the large clinical-isolate dataset, single unique `Y181C` had median fold change `1.6` (`1.2–1.8`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

The structural basis is explicit in the `DOR` versus `RPV` comparison. `RPV` stacks with `Y181`, but `DOR` instead stacks with `Y188`, while `Y181` is shifted toward the rim of the pocket in the `DOR` complex ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). Smith et al. state that the lack of interaction with `Y181` explains why `Y181C` remains susceptible to `DOR`.

### G190A
`G190A` is a common transmitted NNRTI mutation, but `DOR` usually retains activity against it. `DOR` suppressed `G190A` in resistance-selection experiments at clinically relevant concentrations, and the DRIVE-BEYOND and DRIVE-SHIFT data supported continued antiviral activity in viruses containing `K103N`, `Y181C`, and/or `G190A` ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/); [Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)). In the large clinical-isolate analysis, single unique `G190A` had median fold change `1.2` (`1.0–1.4`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

No `G190A-DOR` co-structure was identified here. The available structural interpretation is indirect: `G190` is not among the residues highlighted as making the defining direct contacts with `DOR` in the wild-type structure, and `DOR` was specifically designed to retain activity against common transmitted mutations including `G190A` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/); [Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

### V106I
`V106I` is best supported as a polymorphism rather than a canonical `DOR` resistance mutation. In the 2024 JID study, site-directed `V106I` remained below the `DOR` biological cutoff, whereas `V106A`, `V106M`, and `Y188L` did not ([Giammarino et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38206187/)). In clinically derived `V106I` viruses, median fold change was `1.2` in subtype B and `1.8` in non-B viruses ([Giammarino et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38206187/)). In the large clinical-isolate dataset, single unique `V106I` had median fold change `0.8` (`0.6–1.3`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

Clinical-development analyses reached the same conclusion. Martin et al. state that available data support `V106I` as a polymorphism rather than a `DOR` resistance-associated substitution, noting that `V106I` does not reduce `DOR` potency in vitro and was present in baseline viruses from treatment-naive participants who responded to `DOR` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

The structural distinction from `V106A` is straightforward. `V106` lies in direct contact with the `DOR` pyridone core ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). Replacing valine with isoleucine preserves a branched hydrophobic side chain at that site. Replacing valine with alanine removes side-chain bulk and weakens that contact.

## High-Confidence Positive Controls

### V106A
`V106A` is a canonical `DOR` resistance mutation. In `DOR` resistance-selection experiments, `V106A` was the starting point for the major subtype B and subtype A pathways, followed by `F227L` or `L234I` as `DOR` concentration increased ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/); [Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)). The same review reports approximately `10-fold` reduced susceptibility for `V106A` alone and `>150-fold` for `V106A/L234I` and `V106A/F227L` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

The structural mechanism is direct. Smith et al. state that the branched hydrophobic side chain of `V106` interacts with the `DOR` pyridone core, and that this contact is lost in `V106A` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).

### Y188L
`Y188L` is a high-confidence single `DOR` resistance mutation. In the large clinical-isolate analysis, single unique `Y188L` had median fold change `41.0` (`25.0–250.0`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)). A subtype C phenotyping study also found high-level `DOR` resistance for `Y188L` ([Reddy et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11437401/)).

The structural mechanism is described in detail in [y188l_resistance_mechanism.md](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/manuscript/y188l_resistance_mechanism.md). Briefly, `DOR` stacks its chlorophenol ring against `Y188`; `Y188L` removes that aromatic surface and disrupts a key stabilizing interaction ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/); [Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)).

### Y318F
`Y318F` is one of the few single substitutions that crosses the `DOR` biological cutoff in the clinical-isolate dataset. Single unique `Y318F` had median fold change `11.0` (`3.0–14.1`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)). In clinical-development follow-up, site-directed `Y318F` conferred a `9.9-fold` reduction in `DOR` susceptibility ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

`Y318` lies in the distal `3'` region of RT and has long been linked to NNRTI resistance more broadly ([Harrigan et al., 2002](https://pmc.ncbi.nlm.nih.gov/articles/PMC136283/)). In the `RPV` complex, `Y318` stacks with the benzonitrile group, whereas `DOR` uses a different deeper binding geometry centered more strongly on `Y188`, `V106`, and `P236` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). A direct mutant co-structure explaining `Y318F-DOR` was not identified here, but the phenotypic evidence supports `Y318F` as a single-mutation positive control.

## Combination Pathways With Strong Doravirine Resistance

### V106A + F227L
`V106A/F227L` is a high-confidence `DOR` resistance pathway. `V106A` was the first-step mutation in selection experiments, and `F227L` emerged as `DOR` concentration increased ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)). Review summaries report `>150-fold` resistance for `V106A/F227L` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).

The mechanistic interpretation follows the wild-type structure. `V106` contacts the `DOR` pyridone core, while `F227` lies in the upper hydrophobic portion of the pocket that is displaced toward `DOR` in the `DOR` complex relative to the `RPV` complex ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). The combination therefore perturbs both the central core contact (`V106`) and the upper-pocket packing environment (`F227`).

### V106A + L234I
`V106A/L234I` is another major `DOR` resistance pathway selected in vitro after the initial `V106A` step ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)). Review summaries report `>150-fold` resistance for this combination ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).

`L234` is part of the same upper pocket region that is shifted toward `DOR` in the wild-type structure and that surrounds the distal portion of the inhibitor ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). As with `V106A/F227L`, the mechanistic pattern is cooperative disruption of both a direct core contact and the upper-pocket contour that helps accommodate the distal end of `DOR`.

### V106A + P225H
`V106A/P225H` is supported as a resistant combination in published susceptibility summaries. In the transmitted-resistance review, `V106A/P225H` is listed with `>64-fold` reduced `DOR` susceptibility ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)). Clinical-development data also reported a resistant isolate containing `V106A/P225H/Y318F/K65R` with `>210-fold` reduction in `DOR` susceptibility ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

The structural rationale is that `P225` lies in the upper portion of the pocket that shifts between the `RPV` and `DOR` complexes, near the region used by the distal substituents of `DOR` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). `V106A` removes a direct core contact, and `P225H` alters a nearby upper-pocket residue within the same `DOR`-dependent region.

### V106I + F227C
`V106I/F227C` is one of the clearest examples of a combination in which a weak or polymorphic single substitution becomes strongly resistant in a pairwise context. A `DOR`-resistant clinical isolate with `V106I/F227C` showed `>105-fold` reduced susceptibility ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)). Martin et al. note that `V106I` emerged in clinical failure often in combination with `F227C` rather than as an isolated `DOR`-selected substitution ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

The structural interpretation is again cooperative. `V106I` alone preserves the hydrophobic side chain at a direct-contact residue and usually has little phenotypic effect. `F227C` alters an upper-pocket residue that is part of the hydrophobic tunnel used by `DOR`. Together, these substitutions perturb both the `V106` contact region and the `F227`-containing upper-pocket environment.

### A98G + F227C
Clinical-development data include an `A98G/F227C/M184V` isolate with `>93-fold` reduced susceptibility to `DOR` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)). In the review of transmitted resistance, `A98G` is listed among substitutions that do not reduce `DOR` susceptibility by themselves but contribute to strong resistance in combination with canonical `DOR` resistance mutations such as `F227C` ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).

The residue-specific structural mechanism for `A98G` in the `DOR` complex is less directly defined than for `V106`, `Y188`, or `F227`. The available evidence supports a combination model in which `F227C` provides the canonical upper-pocket resistance component and `A98G` contributes an additional local change in pocket geometry or conformational preference.

## Special Case: F227C
`F227C` does not fit cleanly with the other high-confidence negative controls. Published reviews classify `F227C/L` among the canonical `DOR` resistance-associated substitutions that are associated with the greatest reductions in susceptibility ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)). Clinical-development analyses also report that substitutions at positions `106` and `227` were the most prevalent emergent `DOR` substitutions in virologic failure, with `F227C` commonly occurring alongside `V106I/A/M` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

At the same time, the clinical failure patterns usually involved `F227C` in combination with other substitutions rather than as the dominant isolated pathway. The strongest direct clinical-development examples are combination isolates such as `V106I/F227C` and `A98G/F227C/M184V` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).

The structural reason `F227C` is mechanistically plausible is that `F227` is part of the upper hydrophobic portion of the pocket that accommodates the distal end of `DOR` and shifts toward the inhibitor in the `DOR` complex ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). Replacing phenylalanine with cysteine reduces aromatic bulk and changes the local packing surface within a region that contributes directly to `DOR` binding. Direct `F227C-DOR` structural data were not identified here.

## References
- Asante-Appiah E, Lai J, Wan H, Yang D, Martin EA, Sklar P, Hazuda D, Petropoulos CJ, Walworth C, Grobler JA. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- de Béthune M-P, et al. Pharmaceutical, clinical, and resistance information on doravirine, a novel non-nucleoside reverse transcriptase inhibitor for the treatment of HIV-1 infection. 2020 review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)
- Feng M, Wang D, Grobler JA, Hazuda DJ, Miller MD, Lai M-T. In vitro resistance selection with doravirine (MK-1439), a novel nonnucleoside reverse transcriptase inhibitor with distinct mutation development pathways. *Antimicrob Agents Chemother*. 2015. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)
- Feng M, Sachs NA, Xu M, Grobler J, Blair W, Hazuda DJ, Miller MD, Lai M-T. Doravirine suppresses common nonnucleoside reverse transcriptase inhibitor-associated mutants at clinically relevant concentrations. *Antimicrob Agents Chemother*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)
- Giammarino F, et al. Prevalence and Phenotypic Susceptibility to Doravirine of the HIV-1 Reverse Transcriptase V106I Polymorphism in B and Non-B Subtypes. *J Infect Dis*. 2024. [PubMed](https://pubmed.ncbi.nlm.nih.gov/38206187/)
- Harrigan PR, Salim M, Stammers DK, Wynhoven B, Brumme ZL, McKenna P, Larder B, Kemp SD. A mutation in the 3' region of the human immunodeficiency virus type 1 reverse transcriptase (Y318F) associated with nonnucleoside reverse transcriptase inhibitor resistance. *J Virol*. 2002. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC136283/)
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Reddy N, Papathanasopoulos M, Steegen K, Basson AE. K103N, V106M and Y188L Significantly Reduce HIV-1 Subtype C Phenotypic Susceptibility to Doravirine. *Viruses*. 2024. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11437401/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
