# G190A Doravirine Susceptibility Mechanism

## Summary
`G190A` is a common NNRTI mutation, but doravirine (`DOR`) usually retains activity against it. Published phenotypic data place `G190A` among the mutations with little or no effect on `DOR` susceptibility. In the current MD dataset, `G190A` does not materially disrupt the deeper `DOR` anchoring network centered on `K103`, `V106`, and `Y188`, and it does not strongly displace the neighboring `β9-β10` hairpin residues such as `VAL179`.

## Mechanistic Model
1. The defining wild-type `DOR` contacts involve `Y188`, `V106`, `L100`, `P236`, and the main-chain region near `K103`.
2. `G190` is not highlighted as a defining direct `DOR` contact in the wild-type structure.
3. The `G190A` substitution lies in the `β9-β10` hairpin of the p66 palm subdomain, but the current MD data indicate that the key `DOR`-supporting contacts remain largely intact and that the neighboring `VAL179` side of this hairpin remains close to the ligand.
4. `DOR` susceptibility is therefore usually preserved.

## Evidence Specific To Doravirine
- `DOR` suppressed `G190A` in resistance-selection experiments at clinically relevant concentrations ([Feng et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)).
- Clinical-development reviews report maintained antiviral activity of `DOR` in viruses containing common transmitted mutations including `G190A` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).
- In the large clinical-isolate analysis, single unique `G190A` had median fold change `1.2` (`1.0–1.4`) ([Asante-Appiah et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)).

## Structural Basis
The published wild-type `DOR` structure emphasizes contacts with `Y188`, `V106`, `L100`, `P236`, and the `K103` main-chain region rather than with `G190` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). In HIV-1 RT structural nomenclature, residues `179`, `181`, `188`, and `190` are commonly described as part of the `β9-β10 hairpin` in the p66 palm subdomain ([Ren et al., 2004 review excerpt summarized in](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937); [Das et al., 2004/2005 review summary](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)). The current MD dataset provides mutant-specific structural evidence consistent with the same interpretation.

## Evidence From The Current MD Dataset
The local `G190A-DOR` trajectories support preservation of the main `DOR` binding geometry:

- In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, the ligand remains close to the wild-type pose:
  - `ligand_pose_rmsd_angstrom`
  - `WT`: `1.276 Å`
  - `G190A`: `1.418 Å`
- In the same file, the key contact-region distances remain similar to `WT`:
  - `residue_min_distance_VAL106_angstrom`
  - `WT`: `3.389 Å`
  - `G190A`: `3.238 Å`
  - `residue_min_distance_TYR188_angstrom`
  - `WT`: `3.284 Å`
  - `G190A`: `3.190 Å`
  - `residue_min_distance_LYS103_angstrom`
  - `WT`: `2.989 Å`
  - `G190A`: `3.018 Å`
- The local distance to residue `190` itself is also very similar:
  - `residue_min_distance_GLY190_angstrom`
  - `WT`: `3.323 Å`
  - `G190A`: `3.323 Å`
- The neighboring `VAL179` side of the same `β9-β10` hairpin also remains close to `DOR`:
  - `residue_min_distance_VAL179_angstrom`
  - `WT`: `3.514 Å`
  - `G190A`: `3.438 Å`
- In `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv`, the `VAL179` hydrophobic contacts remain much closer to the wild-type range in `G190A` than in the more resistant `G190E` and `G190S` mutants:
  - `hydrophobic: VAL179:CG1 - 2KW:F15` mean distance in `G190A`: `4.68 Å`
  - `hydrophobic: VAL179:CB - 2KW:F15` mean distance in `G190A`: `5.27 Å`
- In `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv`, the key polar `LYS103` contacts remain highly occupied in `G190A`:
  - `polar: LYS103:O - 2KW:N19` reference-like occupancy: `0.983`
  - `polar: LYS103:O - 2KW:N20` reference-like occupancy: `0.863`

These data support a model in which `G190A` causes, at most, a modest local perturbation of the p66 `β9-β10` hairpin while preserving both the `VAL179` side of that hairpin and the deeper `DOR` anchoring interactions that are most important for susceptibility.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` remains close to wild type:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `G190A`: mean `-150.23`, median `-144.11`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `G190A`: median `-267.14`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `G190A`: median `-26.09`

This is consistent with the structural interpretation that `G190A` causes only a limited perturbation of the p66 `β9-β10` hairpin and does not strongly destabilize `DOR` binding.

## References
- Asante-Appiah E, et al. Impact of HIV-1 Resistance-Associated Mutations on Susceptibility to Doravirine: Analysis of Real-World Clinical Isolates. *Antimicrob Agents Chemother*. 2021. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597775/)
- Feng M, Sachs NA, Xu M, Grobler J, Blair W, Hazuda DJ, Miller MD, Lai M-T. Doravirine suppresses common nonnucleoside reverse transcriptase inhibitor-associated mutants at clinically relevant concentrations. *Antimicrob Agents Chemother*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4808216/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Ren J, Esnouf R, Garman E, et al. Structural and biochemical effects of human immunodeficiency virus mutants resistant to non-nucleoside reverse transcriptase inhibitors. *Structure*. 2004 review excerpt. [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Das K, Clark AD Jr, Lewi PJ, et al. Conformational changes in HIV-1 reverse transcriptase induced by nonnucleoside reverse transcriptase inhibitor binding. review summary of NNIBP architecture. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)
