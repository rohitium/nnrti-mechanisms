# V106A Resistance Mechanism

## Summary
`V106A` is a canonical doravirine (`DOR`) resistance mutation. It is a first-step mutation in in vitro `DOR` resistance selection and produces additional high-resistance pathways when combined with substitutions such as `F227L` or `L234I`. The structural mechanism includes loss of hydrophobic bulk at residue `106` and a local repacking of the `105-106` pocket wall that shifts `DOR` toward `SER105`.

## Mechanistic Model
1. In the wild-type RT-`DOR` structure, the pyridone core of `DOR` contacts `V106`.
2. `V106A` replaces valine with alanine and removes side-chain bulk at that contact site.
3. In the current MD dataset, this local cavity is accompanied by a shift of `DOR` toward `SER105`, indicating repacking of the `105-106` wall rather than simple loss of a single contact.
4. Additional upper-pocket substitutions can then amplify resistance further.

## Evidence Specific To Doravirine
- In `DOR` resistance-selection experiments, `V106A` was the starting point for major subtype B and subtype A pathways ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)).
- Review summaries report approximately `10-fold` reduced susceptibility for `V106A` alone and much larger reductions when paired with `F227L` or `L234I` ([Martin et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)).
- `V106A` is repeatedly listed among the canonical substitutions associated with the greatest reductions in `DOR` susceptibility ([Tang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)).

## Structural Basis
Smith et al. state that the branched hydrophobic side chain of `V106` interacts with the `DOR` pyridone core and that this contact is lost in `V106A` ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).

## Evidence From The Current MD Dataset
The local `V106A-DOR` trajectories indicate both weaker residue-`106` packing and a shift toward `SER105`:

- Direct residue-`106` sidechain-to-`DOR` distances increase relative to `V106I` and are slightly higher than wild type on average:
  - `WT`: mean `3.36 Å`
  - `V106A`: mean `3.54 Å`
  - `V106I`: mean `3.22 Å`
- `SER105` moves substantially closer to `DOR`:
  - `SER105` sidechain-to-`DOR`
  - `WT`: mean `6.65 Å`
  - `V106A`: mean `5.30 Å`
  - `V106I`: mean `7.15 Å`
  - `SER105 OG`-to-`DOR`
  - `WT`: mean `7.44 Å`
  - `V106A`: mean `5.98 Å`
  - `V106I`: mean `7.95 Å`
- In the precomputed feature table `results/analysis/ligand_pocket_features/tables/frame_features.csv`, `V106A` also shows a weaker residue-`106` contact than `V106I`:
  - `residue_min_distance_VAL106_angstrom`
  - `WT`: `3.389 Å`
  - `V106A`: `3.468 Å`
  - `V106I`: `3.228 Å`
- The same table shows that `V106A` places `DOR` closer to the palm side of the pocket:
  - `ligand_palm_distance_angstrom`
  - `WT`: `2.991 Å`
  - `V106A`: `2.637 Å`
  - `V106I`: `3.186 Å`

These data support a model in which `V106A` does not merely remove a single `V106-DOR` contact. Instead, it repacks the local `105-106` wall and induces a compensatory shift of `DOR` toward `SER105`, producing a new geometry that is less favorable than the wild-type or `V106I` arrangement.

The lysine-side distances support the same interpretation:

- In `results/analysis/ligand_pocket_features/tables/frame_features.csv`, `V106A` shifts `DOR` away from `LYS101`:
  - `residue_min_distance_LYS101_angstrom`
  - `WT`: `3.513 Å`
  - `V106A`: `3.888 Å`
  - `V106I`: `3.504 Å`
- `LYS101 > 4.0 Å` occupancy is much higher in `V106A` than in `WT` or `V106I`:
  - `WT`: `0.023`
  - `V106I`: `0.081`
  - `V106A`: `0.374`
- At the same time, the deeper `LYS103` polar contacts in `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv` are not lost in `V106A`:
  - `polar: LYS103:O - 2KW:N19` reference-like occupancy
  - `V106A`: `0.978`
  - `V106I`: `0.680`

This combination of observations indicates that `V106A` shifts the ligand away from the `LYS101` side of the local `101/105/106` wall while preserving a workable deeper `LYS103` anchor. The structural effect is therefore better described as local pose repacking than as uniform loss of all lysine-associated contacts.

The MM/GBSA analysis is consistent with a reproducibly less favorable bound state in `V106A`. In `results/tables/holo/mmgbsa_replicate_metrics.csv`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `V106A`: mean `-144.78`, median `-137.94`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `V106A`: median `-267.06`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `V106A`: median `-24.47`

So the energetic penalty in `V106A` is not limited to one term. The mutant shows weaker overall binding, with less favorable van der Waals packing and less favorable electrostatics than wild type. This matches the structural picture of a locally repacked `105-106` wall that produces a compensatory but suboptimal ligand geometry.

## References
- Feng M, Wang D, Grobler JA, Hazuda DJ, Miller MD, Lai M-T. In vitro resistance selection with doravirine (MK-1439), a novel nonnucleoside reverse transcriptase inhibitor with distinct mutation development pathways. *Antimicrob Agents Chemother*. 2015. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)
- Martin EA, Lai MT, Ngo W, et al. Review of Doravirine Resistance Patterns Identified in Participants During Clinical Development. *JAIDS*. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Rilpivirine and Doravirine have complementary efficacies against NNRTI-resistant HIV-1 mutants. *Retrovirology*. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Tang MW, et al. Potential role of doravirine for the treatment of HIV-1-infected persons with transmitted drug resistance. *AIDS Res Ther*. 2023. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903540/)
