# Y188L Resistance Mechanism

## Summary
`Y188L` is a mechanistically anchored doravirine resistance mutation. The core explanation is that residue 188 is part of the aromatic wall of the HIV-1 RT NNRTI-binding pocket, and doravirine relies on an aromatic interaction between its chlorophenol moiety and `Y188`. Replacing tyrosine with leucine removes that aromatic surface and can also alter the local pocket geometry. The result is a large loss of favorable binding interactions and a marked reduction in doravirine susceptibility.

## Mechanistic Model
The proposed mechanism can be written as a short causal chain:

1. In wild-type RT, `Y188` contributes an aromatic face in the NNRTI-binding pocket.
2. Doravirine positions its chlorophenol ring against this aromatic wall.
3. `Y188L` removes the aromatic ring at position 188, replacing it with an aliphatic leucine side chain.
4. This disrupts aromatic stacking and can locally reshape the pocket.
5. Doravirine binding is therefore weakened, leading to high-level phenotypic resistance.

## Evidence Specific To Doravirine
Several lines of evidence support this interpretation:

- The wild-type RT-doravirine crystal structure places doravirine deep in the NNRTI pocket and supports a direct role for `Y188` in stabilizing the ligand, including the chlorophenol-containing portion of the inhibitor ([Smith et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)).
- Merck’s doravirine resistance-selection paper explicitly states that `Y188L` removes the `pi-pi` interaction with doravirine and introduces unfavorable steric consequences, explaining the large resistance effect ([Feng et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)).
- Clinical and resistance reviews consistently list `Y188L` among the strongest single doravirine-resistance mutations, often with very large fold-change values ([Fashanu et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/); [de Béthune et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)).
- A recent subtype C phenotyping study reported very high doravirine resistance for `Y188L`, around `89-fold` in that assay system ([Reddy et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11437401/)).

## Local Trajectory Comparison With Y181C
The local `WT`, `Y181C`, and `Y188L` simulations support the same qualitative distinction reported by Smith et al. `Y181C` preserves a `WT`-like `DOR` pose, whereas `Y188L` shows a larger disruption on the `188` side of the pocket and a shift of the ligand toward the pocket entrance.

In `results/analysis/ligand_pocket_features/tables/frame_features.csv`:

- `ligand_pose_rmsd_angstrom` is larger in `Y188L` than in `Y181C` (`1.542` vs `1.414 A`).
- `ligand_palm_distance_angstrom` is larger in `Y188L` than in `Y181C` (`3.257` vs `2.709 A`).
- `ligand_entrance_distance_angstrom` is smaller in `Y188L` than in `Y181C` (`7.474` vs `7.997 A`).

Together, those values indicate that `Y188L` moves `DOR` closer to the entrance/rim and farther from the palm than `Y181C`.

The largest local disruption is on the `188` side of the pocket:

- `residue_min_distance_TYR188_angstrom` is `3.603 A` in `Y188L` versus `3.098 A` in `Y181C`.
- The fraction of frames with `residue_min_distance_TYR188_angstrom > 4 A` increases from `0.054` in `Y181C` to `0.153` in `Y188L`.

The atom-level contact traces in `results/tables/holo/dor_key_contacts_timeseries_all_mutations.csv` show the same pattern. Reference-like `TYR188` contact occupancies are lower in `Y188L` than in `Y181C`:

- `TYR188:CD2-C`: `0.131` in `Y188L` vs `0.620` in `Y181C`
- `TYR188:C-F`: `0.172` vs `0.328`
- `TYR188:CB-F14`: `0.117` vs `0.330`

Neighboring support contacts are also weaker in `Y188L`, including `VAL189:C-F` (`0.296` vs `0.396`) and `LYS103:O-N20` (`0.630` vs `0.874`).

These trajectory-derived differences are consistent with a `Y188L` mechanism in which loss of the aromatic side chain at position `188` weakens the chlorophenol-facing wall of the pocket, shifts `DOR` toward the entrance, and partially destabilizes the surrounding contact network. This pattern is not seen in `Y181C`.

## Broader NNRTI Context
This mechanism is also consistent with older NNRTI literature. `Y188` has long been recognized as part of the aromatic clamp that stabilizes NNRTI binding. Structural and biochemical studies with other NNRTIs showed that mutations at `Y181` and `Y188` often reduce susceptibility by disrupting aromatic stacking with inhibitor rings.

One mechanistic study is especially important because it directly tested aromaticity rather than just mutation identity. Substituting tyrosine with phenylalanine at `Y188` preserves aromatic character, and those variants retained or even improved NNRTI susceptibility. That result strongly supports the conclusion that aromatic `pi-pi` interactions at this site are a major contributor to NNRTI binding ([Maga et al., 2016](https://www.mdpi.com/1999-4915/8/10/263)).

## Energy Analysis
The MM/GBSA summary in `results/tables/holo/mmgbsa_replicate_metrics.csv` is consistent with marked energetic destabilization in `Y188L`:

- `binding_dg`
  - `WT`: mean `-152.38`, median `-148.67`
  - `Y188L`: mean `-108.83`, median `-110.57`
- `binding_dg_vdw`
  - `WT`: median `-272.24`
  - `Y188L`: median `-255.51`
- `binding_dg_electrostatic`
  - `WT`: median `-28.09`
  - `Y188L`: median `-18.53`
- `binding_dg_gb`
  - `WT`: median `148.21`
  - `Y188L`: median `161.10`

This pattern is consistent with loss of the aromatic `Y188` interaction and a less favorable bound state for the chlorophenol-containing end of `DOR`.

## References
- Feng M, Wang D, Grobler JA, Hazuda DJ, Miller MD, Lai M-T. In vitro resistance selection with doravirine (MK-1439), a novel nonnucleoside reverse transcriptase inhibitor with distinct mutation development pathways. *Antimicrob Agents Chemother*. 2015. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4291404/)
- Smith SJ, Zhao Y, Burke TR Jr, Hughes SH. Efficacies of doravirine and rilpivirine in vitro against wild-type and drug-resistant HIV-1. *Antimicrob Agents Chemother* / structural discussion in associated RT-DOR analyses. 2016. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4942337/)
- Maga G, et al. Mechanistic study of common non-nucleoside reverse transcriptase inhibitor-resistant mutations with K103N and Y181C substitutions. *Viruses*. 2016. [MDPI](https://www.mdpi.com/1999-4915/8/10/263)
- Fashanu OE, et al. Potential role of doravirine for treatment of HIV-1-infected persons with transmitted drug resistance. 2020 review of doravirine resistance data. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7655028/)
- de Béthune M-P, et al. Review of doravirine resistance pathways and NNRTI cross-resistance. 2020. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7055513/)
- Reddy T, et al. K103N, V106M and Y188L significantly reduce HIV-1 subtype C phenotypic susceptibility to doravirine. 2024. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11437401/)
