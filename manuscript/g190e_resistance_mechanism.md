# G190E Resistance Mechanism

## Summary
`G190E` is associated with substantially reduced doravirine (`DOR`) susceptibility. In the current dataset, `G190E` behaves differently from `G190A` because it strongly perturbs the p66 `β9-β10` hairpin and propagates that perturbation into the `VAL179` side of the NNRTI-binding pocket and into the deeper `LYS103`-linked contact network.

## Mechanistic Model
1. Residue `190` lies in the p66 `β9-β10` hairpin of HIV-1 reverse transcriptase.
2. `DOR` depends on a deeper binding geometry stabilized by contacts involving `Y188`, `V106`, `P236`, and the `K103` backbone region.
3. `G190E` introduces a larger, charged side chain into the `β9-β10` hairpin.
4. In the current MD dataset, this change is associated with strong displacement of the neighboring `VAL179` side of the hairpin and weakening of key `LYS103` polar contacts.
5. The result is a less favorable binding geometry and reduced `DOR` susceptibility.

## Evidence Specific To Doravirine
- In the susceptibility dataset used in this project, `G190E` has `18.0-fold` reduced `DOR` susceptibility.
- Published doravirine literature identifies `G190E` among the substitutions that can reduce `DOR` susceptibility, although the evidence base is more limited than for canonical mutations such as `V106A` or `Y188L` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).

## Structural Basis
The wild-type RT-`DOR` structure emphasizes contacts with `Y188`, `V106`, `L100`, `P236`, and the `K103` main-chain region rather than with `G190` directly ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). In HIV-1 RT structural nomenclature, residues `179`, `181`, `188`, and `190` are commonly described as part of the p66 `β9-β10 hairpin` ([Ren et al., 2004](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937); [Das et al., 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)).

## Evidence From The Current MD Dataset
The local `G190E-DOR` trajectories support a stronger perturbation than `G190A` or `G190S`:

- In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, the `VAL179` side of the p66 `β9-β10` hairpin moves much farther from `DOR`:
  - `residue_min_distance_VAL179_angstrom`
  - `WT`: `3.514 Å`
  - `G190E`: `5.737 Å`
  - `>4.0 Å` occupancy in `G190E`: `0.961`
- In the same file, the deeper contact network is also perturbed:
  - `residue_min_distance_LYS103_angstrom`
  - `WT`: `2.989 Å`
  - `G190E`: `3.163 Å`
  - `residue_min_distance_TYR181_angstrom`
  - `WT`: `3.328 Å`
  - `G190E`: `3.616 Å`
  - `residue_min_distance_PRO236_angstrom`
  - `WT`: `3.240 Å`
  - `G190E`: `3.460 Å`
- In `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv`, the `VAL179` hydrophobic contacts are strongly weakened:
  - `hydrophobic: VAL179:CG1 - 2KW:F15`
  - mean distance: `7.15 Å`
  - reference-like occupancy: `0.020`
  - `hydrophobic: VAL179:CB - 2KW:F15`
  - mean distance: `7.49 Å`
  - reference-like occupancy: `0.002`
- In the same file, the `LYS103` polar contacts are weakened relative to `G190A` and `G190S`:
  - `polar: LYS103:O - 2KW:N19` reference-like occupancy: `0.816`
  - `polar: LYS103:O - 2KW:N20` reference-like occupancy: `0.389`
  - `polar: LYS103:N - 2KW:N19` reference-like occupancy: `0.232`
- In `results/tables/holo/mmgbsa_replicate_metrics.csv`, `G190E` shows substantially weaker overall binding than `WT`:
  - `WT`: mean `binding_dg = -152.38`
  - `G190E`: mean `binding_dg = -122.80`

These data support a mechanism in which `G190E` destabilizes `DOR` binding by strongly displacing the p66 `β9-β10` hairpin, especially the `VAL179` side of that hairpin, and by weakening the deeper `LYS103`-linked anchoring interactions.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` is consistent with a large energetic penalty:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `G190E`: mean `-122.80`, median `-125.12`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `G190E`: median `-32.44`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `G190E`: median `182.86`

This pattern indicates that the net penalty is dominated less by van der Waals failure than by a much less favorable polar/solvation balance, which matches the strong `β9-β10` hairpin repatterning seen structurally.

## References
- Das K, Clark AD Jr, Lewi PJ, et al. Conformational changes in HIV-1 reverse transcriptase induced by nonnucleoside reverse transcriptase inhibitor binding. review summary of NNIBP architecture. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Ren J, Esnouf R, Garman E, et al. Structural and biochemical effects of human immunodeficiency virus mutants resistant to non-nucleoside reverse transcriptase inhibitors. *Structure*. 2004 review excerpt. [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
