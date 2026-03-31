# G190S Resistance Mechanism

## Summary
`G190S` is associated with modestly reduced doravirine (`DOR`) susceptibility. In the current dataset, `G190S` behaves as an intermediate case between `G190A` and `G190E`: it perturbs the p66 `β9-β10` hairpin and weakens the `VAL179`/aromatic-wall side of the pocket more than `G190A`, but it does not disrupt the deeper `LYS103` anchor as strongly as `G190E`.

## Mechanistic Model
1. Residue `190` lies in the p66 `β9-β10` hairpin.
2. `DOR` binds deeper in the NNRTI pocket and depends on a contact network involving `Y188`, `V106`, `P236`, and the `K103` backbone region.
3. `G190S` introduces a polar side chain at position `190`.
4. In the current MD dataset, this is associated with moderate displacement of the `VAL179`/`TYR181` side of the `β9-β10` hairpin and partial weakening of `LYS103`-linked contacts.
5. The resulting structural perturbation is larger than in `G190A` but smaller than in `G190E`, which is consistent with modest rather than high-level resistance.

## Evidence Specific To Doravirine
- In the susceptibility dataset used in this project, `G190S` has `5.2-fold` reduced `DOR` susceptibility.
- Published reviews place `G190S` among noncanonical substitutions that can reduce `DOR` susceptibility, although the evidence base is less extensive than for canonical mutations such as `V106A` or `Y188L` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).

## Structural Basis
The published wild-type `DOR` structure emphasizes deeper contacts with `Y188`, `V106`, `P236`, and the `K103` region rather than direct dependence on `G190` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). Residues `179`, `181`, `188`, and `190` are part of the p66 `β9-β10 hairpin` in standard HIV-1 RT structural nomenclature ([Ren et al., 2004](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937); [Das et al., 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)).

## Evidence From The Current MD Dataset
The local `G190S-DOR` trajectories support an intermediate level of structural perturbation:

- In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, the `VAL179` side of the p66 `β9-β10` hairpin is displaced farther from `DOR` than in `WT` or `G190A`, but less than in `G190E`:
  - `residue_min_distance_VAL179_angstrom`
  - `WT`: `3.514 Å`
  - `G190S`: `3.638 Å`
  - `G190E`: `5.737 Å`
  - `>4.0 Å` occupancy in `G190S`: `0.265`
- The aromatic-wall side is also perturbed:
  - `residue_min_distance_TYR181_angstrom`
  - `WT`: `3.328 Å`
  - `G190S`: `3.780 Å`
  - `residue_min_distance_TYR188_angstrom`
  - `WT`: `3.284 Å`
  - `G190S`: `3.425 Å`
- The deeper `LYS103` anchor remains much closer to wild type than in `G190E`:
  - `residue_min_distance_LYS103_angstrom`
  - `WT`: `2.989 Å`
  - `G190S`: `2.986 Å`
  - `G190E`: `3.163 Å`
- In `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv`, the `VAL179` hydrophobic contacts are weakened, but not to the extent seen in `G190E`:
  - `hydrophobic: VAL179:CG1 - 2KW:F15`
  - mean distance: `5.36 Å`
  - reference-like occupancy: `0.154`
  - `hydrophobic: VAL179:CB - 2KW:F15`
  - mean distance: `5.53 Å`
  - reference-like occupancy: `0.052`
- The `LYS103` polar contacts are partially weakened but remain stronger than in `G190E`:
  - `polar: LYS103:O - 2KW:N19` reference-like occupancy: `0.948`
  - `polar: LYS103:O - 2KW:N20` reference-like occupancy: `0.795`
  - `polar: LYS103:N - 2KW:N19` reference-like occupancy: `0.360`
- In `results/tables/holo/mmgbsa_replicate_metrics.csv`, overall binding is weaker than `WT` but much less impaired than in `G190E`:
  - `WT`: mean `binding_dg = -152.38`
  - `G190S`: mean `binding_dg = -140.01`
  - `G190E`: mean `binding_dg = -122.80`

These data support a mechanism in which `G190S` produces a moderate perturbation of the p66 `β9-β10` hairpin that weakens the `VAL179`/aromatic-wall side of the pocket without causing the stronger anchor disruption seen in `G190E`.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` is intermediate between `G190A` and `G190E`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `G190S`: mean `-140.01`, median `-130.77`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `G190S`: median `-264.66`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `G190S`: median `155.14`

This supports a moderate energetic penalty consistent with partial `β9-β10` hairpin distortion rather than the stronger polar/solvation disruption seen in `G190E`.

## References
- Das K, Clark AD Jr, Lewi PJ, et al. Conformational changes in HIV-1 reverse transcriptase induced by nonnucleoside reverse transcriptase inhibitor binding. review summary of NNIBP architecture. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1298242/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Ren J, Esnouf R, Garman E, et al. Structural and biochemical effects of human immunodeficiency virus mutants resistant to non-nucleoside reverse transcriptase inhibitors. *Structure*. 2004 review excerpt. [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S1357272504000937)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
