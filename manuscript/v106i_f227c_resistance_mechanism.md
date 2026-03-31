# V106I+F227C Resistance Mechanism

## Summary
`V106I/F227C` is a strong doravirine (`DOR`) resistance combination. Published clinical-development data report very large reductions in `DOR` susceptibility for this pair, even though `V106I` alone is usually a polymorphism with little effect on `DOR` potency. The current MD dataset supports a mechanism in which the `106-227` corridor opens and the residue `227` contact region becomes less favorable, while the residue `106` side chain remains close to the ligand.

## Mechanistic Model
1. `V106I` by itself largely preserves the hydrophobic side chain at a direct-contact residue and usually has little phenotypic effect.
2. `F227C` changes the upper hydrophobic portion of the pocket used by the distal end of `DOR`.
3. In the double mutant, the sidechain separation between positions `106` and `227` increases and the residue `227` side chain spends more time farther from `DOR` than in the corresponding wild-type background.
4. The net effect is upper-pocket repacking and weaker overall binding, rather than a simple loss of the residue `106` sidechain-to-ligand contact.

## Evidence Specific To Doravirine
- A clinical-development isolate containing `V106I/F227C` showed `>105-fold` reduced susceptibility to `DOR` ([Lai et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)).
- Martin et al. note that `V106I` emerged in virologic failure often in combination with `F227C`, whereas available data support `V106I` alone as a polymorphism rather than a `DOR` resistance mutation ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).

## Structural Basis
The wild-type `DOR` structure places `V106` in direct contact with the pyridone core and places `F227` in the upper hydrophobic pocket surrounding the distal end of the inhibitor ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)). No direct `V106I/F227C-DOR` co-structure was identified here.

## Evidence From The Current MD Dataset
The current trajectory-derived data support part of this mechanism directly and constrain the rest:

- In `results/tables/holo/drm_sidechain_distance_timeseries_all_mutations.csv`, the `V106I+F227C` panel separates the two mutated side chains into `c1` = `V106I` and `c2` = `F227C`. Relative to the matched wild-type panel, the mean `c1_to_c2` distance increases in all three replicates:
  - rep 1: `+1.15 Å`
  - rep 2: `+0.64 Å`
  - rep 3: `+0.83 Å`
- In the same table, the `F227C` sidechain-to-`DOR` distance (`c2_to_dor`) is larger in the double mutant than in the matched wild-type panel in all three replicates:
  - rep 1: `+0.48 Å`
  - rep 2: `+0.08 Å`
  - rep 3: `+0.29 Å`
- Also in the same table, the residue `227` side chain in the double mutant spends substantially more time farther than `4.0 Å` from `DOR` than the matched wild-type residue `227`:
  - mutant `c2_to_dor > 4.0 Å`: `35.8%`
  - wild-type `c2_to_dor > 4.0 Å`: `7.8%`
- By contrast, the residue `106` sidechain-to-`DOR` distance (`c1_to_dor`) does not show a direct weakening. In the double mutant, its mean value is slightly smaller than in the matched wild-type panel in all three replicates, and `c1_to_dor > 4.0 Å` is never observed in the double-mutant traces.
- In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, the minimum ligand distance to residue `227` is more often in a disengaged state in `V106I+F227C` than in `WT` or `V106I`:
  - `residue_min_distance_PHE227_angstrom > 4.0 Å`
  - `WT`: `8.6%`
  - `V106I`: `4.5%`
  - `V106I+F227C`: `35.5%`
- In the same file, the minimum ligand distance to residue `106` does not move outward in the double mutant. The mean `residue_min_distance_VAL106_angstrom` is `3.09 Å` in `V106I+F227C`, compared with `3.39 Å` in `WT` and `3.23 Å` in `V106I`.
- In `results/tables/holo/mmgbsa_replicate_metrics.csv`, the double mutant shows substantially weaker overall binding than `WT`, `V106I`, or `F227C`:
  - `WT`: mean `binding_dg = -152.38`
  - `V106I`: mean `binding_dg = -127.88`
  - `F227C`: mean `binding_dg = -118.07`
  - `V106I+F227C`: mean `binding_dg = -99.89`

Taken together, these data support an upper-pocket repacking model centered on positions `106` and `227`. They support widening of the `106-227` corridor and weaker residue `227` engagement with `DOR`. They do not support a simple interpretation in which the residue `106` sidechain itself loses contact with the ligand.

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` shows the strongest energetic destabilization in this set:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106I+F227C`: mean `-99.89`, median `-101.12`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `V106I+F227C`: median `-265.85`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `V106I+F227C`: median `160.84`

The replicate variance is high, but the net energetic penalty is extreme and directionally consistent with the structural picture of upper-pocket repacking and weakened residue-`227` engagement.

## References
- Lai M-T, Feng M, Xu M, Ngo W, Diamond TL, Hwang C, Grobler JA, Hazuda DJ, Asante-Appiah E. Doravirine and Islatravir Have Complementary Resistance Profiles and Create a Combination with a High Barrier to Resistance. *Antimicrob Agents Chemother*. 2022. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9112941/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
