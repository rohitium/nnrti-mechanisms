# Is the 4NCG DOR pose specific to 4NCG? — evidence from 7Z2G

Drafted 25-Oct-2026. Addresses the reading of Reviewer 1's general comment in
which the objection is about the *structure*, not only about the sampling: that
4NCG may not represent the biologically relevant RT–DOR ensemble.

## What was done

A second experimental RT–DOR structure exists and was not used in the paper:

| | method | resolution | nucleic acid | ligand |
|---|---|---|---|---|
| **4NCG** | X-ray | 2.58 Å | none | 2KW (DOR) |
| **7Z2G** | cryo-EM, single particle | 3.65 Å | 38-mer DNA aptamer | 2KW (DOR) |

Different method, different laboratory, different year, and — importantly for
Reviewer 2's comment 2 — one of the two has nucleic acid bound. If the DOR pose
in 4NCG were a crystal-packing artefact, or an artefact of the binary complex,
7Z2G is where it would show.

Reproduce with:

```bash
PYTHONPATH=src python -m nnrti.cli.compare_dor_crystal_poses
PYTHONPATH=src python -m nnrti.cli.plot_crystal_pose_comparison
pymol -cq data/structures/4ncg-vs-7z2g-nnibp-overlay.pml
```

Outputs: `results/analysis/crystal_pose_comparison/` (summary.md, six tables,
the superposed 7Z2G coordinates, and three figures).

## What it shows

**The ligand.** Superposing on the 16 NNIBP-lining Cα atoms alone, DOR's 29
heavy atoms agree to **0.62 Å** RMSD, with a centroid displacement of **0.30 Å**.
Fitted ligand-on-ligand, the internal conformation agrees to **0.44 Å**. No
moiety is placed differently from the rest: chlorocyanophenyl 0.78 Å,
pyridinone 0.47 Å, triazolinone 0.54 Å, linker 0.87 Å. All four rotatable
dihedrals agree to within 15°, and all eight values (four per structure) fall
within 2.2 circular SD of the WT simulated mean — the simulated ligand
conformation matches both experiments, not just the one it started from.

**The pocket.** NNIBP Cα RMSD **0.49 Å**. Side-chain RMSD across the 15 p66
residues plus Glu138: median **0.80 Å**, and only Lys101 (2.34 Å, a
solvent-exposed rim lysine) and Phe227 (1.60 Å) exceed 1.2 Å. Two of 15
residues change χ1 by more than 30°, both β-branched valines (V106, V189)
rotating between staggered rotamers. Pocket volume measured by the routine used
on every trajectory: **198.3 Å³ (4NCG) vs 194.1 Å³ (7Z2G)**, a 2% difference.

**Without any superposition at all.** The minimum heavy-atom distance from each
NNIBP residue to DOR agrees between the two structures to a median of **0.16 Å**,
largest discrepancy **0.41 Å** (Lys103), Pearson *r* = **0.974**. These numbers
cannot be improved by a favourable choice of fitting atoms, because no fit
enters them. The four NNIBP gate distances agree to within 1.4 Å.

**Where the two structures *do* differ.** They differ, as they must — one is
bound to DNA. With the palm as the reference frame, Cα RMSD by subdomain:
p51 2.23 Å, p66 fingers 1.83 Å, thumb 1.62 Å, RNase H 1.37 Å, palm 1.21 Å,
connection 0.77 Å — and the **NNIBP lining 0.97 Å**, the second-lowest of all.
The differences sit in the subdomains that grip nucleic acid. The drug site is
the part that does not move.

## Draft text for the response

> We agree that the representativeness of the starting structure is a fair
> question, and we have now tested it directly rather than argued it. A second
> RT–DOR structure has been determined independently (PDB 7Z2G, 3.65 Å cryo-EM,
> HIV-1 RT with a 38-mer DNA aptamer and doravirine), by a different method and
> in a different laboratory, and with nucleic acid bound. We compared it with
> 4NCG.
>
> Superposed on the 16 NNIBP-lining Cα atoms, doravirine's 29 heavy atoms agree
> between the two structures to 0.62 Å RMSD, with a 0.30 Å centroid
> displacement; its four rotatable dihedrals agree to within 15°. The pocket
> that holds it agrees to 0.49 Å over the same Cα atoms and to a median
> side-chain RMSD of 0.80 Å, and the NNIBP volume differs by 2% (198 vs 194 Å³)
> by the same measurement used throughout this work. Without invoking any
> superposition, the minimum heavy-atom distance from each NNIBP residue to the
> drug agrees to a median of 0.16 Å (largest difference 0.41 Å, r = 0.974).
>
> The two structures do differ, in the subdomains that bind nucleic acid: with
> the palm as reference, p51, the fingers and the thumb move by 1.6–2.2 Å,
> while the NNIBP lining moves by 0.97 Å. We therefore do not claim that the
> binary and ternary complexes are the same enzyme; we claim, and now show,
> that the non-nucleoside inhibitor binding pocket and the pose of doravirine
> within it are the same in both, so the pose we sample around is a property of
> the drug–pocket complex rather than of one crystal form. The simulated WT
> ensemble is consistent with both experimental structures: all four dihedral
> values from each structure lie within 2.2 circular standard deviations of the
> WT simulated mean.
>
> This is added as Supplementary Figure SX and Supplementary Table SX, and
> noted in the Discussion. It does not address the separate question of
> sampling from independent starting poses, which we discuss honestly as a
> limitation.

## What this does *not* answer

- **R1.2a, independent starting poses.** Unchanged. Two experimental structures
  agreeing says the pose is right; it says nothing about whether 100 ns of MD
  seeded from it explores more than one basin. That limitation stands as
  written.
- **The ternary complex.** We have not simulated RT + DNA + DOR. The agreement
  shown here is between two *static* structures; it makes the binary-complex
  simulation defensible, it does not reproduce the ternary one.
- **7Z2G's resolution.** At 3.65 Å, side-chain positions carry real
  uncertainty, and side-chain-level agreement should not be over-read. The
  Cα-level, ligand-level and distance-level agreement is the robust part. Every
  NNIBP residue within 15 Å of the pocket centre is fully modelled in 7Z2G, so
  no atom was dropped from the comparison.

## Files

| | |
|---|---|
| `data/structures/apo_reference/7Z2G.{cif,pdb}` | the deposition, and a PDB conversion |
| `src/nnrti/cli/compare_dor_crystal_poses.py` | the comparison |
| `src/nnrti/cli/plot_crystal_pose_comparison.py` | the four-panel figure |
| `data/structures/4ncg-vs-7z2g-nnibp-overlay.pml` | the overlay render |
| `results/analysis/crystal_pose_comparison/summary.md` | every number above |

---

# Manuscript-ready text

Numbering follows what is already in use: Supplementary Figures run to 4 and
Supplementary Tables to 7, so this work becomes **SI section 8**,
**Supplementary Figure 5** and **Supplementary Table 8**. The workbook is
`revisions/Supplementary-Table-8.xlsx` (six sheets); the figure is
`results/analysis/crystal_pose_comparison/plots/crystal_pose_comparison.png`,
with the pocket overlay in `nnibp_overlay.png` as an optional panel E.

## SI section 8

> **8  An independent structure of the RT–DOR complex**
>
> Every simulation in this work starts from one structure, 4NCG, and a pose
> that is particular to that structure would propagate into every result. A
> second RT–DOR complex has since been determined by a different method, in a
> different laboratory, and in a different functional state: 7Z2G, a 3.65 Å
> single-particle cryo-EM structure of HIV-1 RT bound to a 38-mer DNA aptamer
> and doravirine. Comparing the two tests whether the pose we sample around is
> a property of the drug and its pocket or of one crystal form. Doravirine is
> the same chemical component in both depositions, so its 29 heavy atoms were
> matched by name, up to the six automorphisms of its bond graph that permute
> the trifluoromethyl fluorines. Side-chain comparisons likewise allow the
> equivalent terminal atoms of Phe, Tyr, Asp, Glu and Arg to exchange, since a
> ring flip is the same structure.
>
> Superposed on the 16 Cα atoms of the NNIBP-lining residues alone, doravirine
> agrees between the two structures to 0.62 Å heavy-atom RMSD, with a centroid
> displacement of 0.30 Å. Fitted ligand on ligand, its internal conformation
> agrees to 0.44 Å, and no ring system is placed differently from the rest
> (chlorocyanophenyl 0.78 Å, pyridinone 0.47 Å, triazolinone 0.54 Å). The four
> rotatable dihedrals agree to within 15°, and all eight values lie within 2.2
> circular standard deviations of the wild-type simulated mean, so the
> simulated ligand conformation is consistent with both experiments and not
> only with the structure it was seeded from. The pocket agrees to 0.49 Å over
> the same Cα atoms and to a median side-chain RMSD of 0.80 Å; only Lys101, a
> solvent-exposed residue at the rim, and Phe227 exceed 1.2 Å, and only the two
> β-branched valines V106 and V189 change χ1 by more than 30°. Measured by the
> routine used on every trajectory in this work, the NNIBP volume is 198 Å³ in
> 4NCG and 194 Å³ in 7Z2G.
>
> These comparisons all follow a superposition, and a superposition can be
> chosen to flatter. The agreement survives without one. The minimum
> heavy-atom distance from each NNIBP residue to doravirine agrees between the
> two structures to a median of 0.16 Å, with a largest discrepancy of 0.41 Å at
> Lys103 and a correlation of r = 0.974 across the sixteen residues; no fitting
> atoms enter that calculation.
>
> The two structures are not identical, and should not be: one has nucleic acid
> bound and the other does not. Taking the p66 palm as the reference frame, p51
> differs by 2.23 Å, the fingers by 1.83 Å and the thumb by 1.62 Å, while the
> NNIBP lining differs by 0.97 Å. The parts of the enzyme that grip nucleic
> acid move; the drug site does not. We therefore do not claim that the binary
> complex simulated here is interchangeable with the ternary complex. We claim
> the narrower thing the comparison supports: the pocket and the pose of
> doravirine within it are the same in both, so the starting geometry of this
> work is not an artefact of a single structure. At 3.65 Å the side-chain
> positions in 7Z2G carry real uncertainty, and the Cα-level, ligand-level and
> distance-level agreement should be weighed above the side-chain-level
> agreement; every NNIBP residue within 15 Å of the pocket centre is fully
> modelled in 7Z2G, so no atom was dropped from the comparison. This analysis
> does not address whether the simulations explore more than one basin around
> that pose, which is discussed separately as a limitation.

## Supplementary Figure 5 caption

> **Supplementary Figure 5. The doravirine pose is the same in two independent
> experimental structures.** 4NCG (X-ray, 2.58 Å, RT + DOR) compared with 7Z2G
> (single-particle cryo-EM, 3.65 Å, RT + DOR + 38-mer DNA aptamer). (A) Cα
> RMSD by subdomain with the p66 palm as the reference frame. The subdomains
> that bind nucleic acid differ most; the NNIBP lining (green) is among the
> least changed. (B) Minimum heavy-atom distance from each NNIBP residue to
> doravirine in each structure. No superposition enters this panel.
> (C) Per-residue side-chain RMSD after superposition on the 16 NNIBP Cα
> atoms, with terminal atoms of symmetric side chains allowed to exchange.
> (D) The four rotatable dihedrals of doravirine: the wild-type simulated
> distribution over 860 frames, with the value from each experimental
> structure. Per-residue and per-dihedral values are given in Supplementary
> Table 8.

Optional panel E, if a structural view is wanted:

> (E) The two ligands after superposition on the NNIBP Cα atoms alone (4NCG
> blue, 7Z2G red), with the NNIBP side chains of 4NCG in grey. Heavy-atom RMSD
> 0.62 Å over 29 atoms.

## Supplementary Table 8 caption

> **Supplementary Table 8. Agreement between the two experimental RT–DOR
> structures.** Summary quantities, per-residue comparison of the 16
> NNIBP-lining residues (Cα deviation, side-chain RMSD, χ1 and χ2, and minimum
> heavy-atom distance to doravirine in each structure), Cα RMSD by subdomain,
> ligand RMSD by moiety, the four rotatable dihedrals, and the four NNIBP gate
> distances. Superpositions use the 16 NNIBP Cα atoms unless stated; the
> residue-to-drug distances and the gate distances use none.

## Discussion sentence

One sentence, to sit with the existing caveat about the sampled pose:

> That pose is not peculiar to the structure we started from: doravirine and
> the residues lining its pocket agree to 0.62 and 0.49 Å respectively between
> 4NCG and 7Z2G, an independent cryo-EM structure of RT with doravirine and
> DNA bound (Supplementary Figure 5), although the sampling limitation above
> applies to the basin around that pose regardless.

If the Discussion also makes the nucleic-acid point (R2.2), a second sentence:

> The same comparison localises what nucleic acid changes: with the p66 palm as
> the reference frame, p51, the fingers and the thumb differ by 1.6–2.2 Å
> between the two structures while the NNIBP lining differs by 0.97 Å, so the
> mechanisms we cannot reach without nucleic acid are those that act through
> those subdomains rather than through the drug site itself.

## Methods addition

> Experimental structures were compared with gemmi; superpositions are
> Kabsch fits on the stated Cα sets with no outlier rejection. Ligand atoms
> were matched by chemical-component name, minimised over the automorphisms of
> the ligand bond graph, and side-chain RMSDs were minimised over the
> equivalent terminal atoms of symmetric residues. Pocket volumes for both
> structures were measured with the routine used on the trajectories; the value
> for 4NCG, 198.3 Å³, is the one already reported for that structure.

## Also worth updating

- Wherever 4NCG is cited as the starting structure, 7Z2G now needs a citation
  too (PDB 7Z2G; cryo-EM, 2022).
- Data availability / PDB accession list, if the manuscript carries one.

## Citation for 7Z2G

Singh, A. K.; De Wijngaert, B.; Bijnens, M.; Uyttersprot, K.; Nguyen, H.;
Martinez, S. E.; Schols, D.; Herdewijn, P.; Pannecouque, C.; Arnold, E.;
Das, K. Cryo-EM structures of wild-type and E138K/M184I mutant HIV-1 RT/DNA
complexed with inhibitors doravirine and rilpivirine. *Proc. Natl. Acad. Sci.
U.S.A.* **2022**, *119*, e2203660119. DOI: 10.1073/pnas.2203660119.
PDB 7Z2G, deposited 26-Feb-2022.

Note for the response: this paper also reports an E138K/M184I mutant RT–DOR
structure. We compared the wild-type entry only, which is the right comparator
for a wild-type starting structure.
Confirmed wild type at the two positions that paper mutates: 7Z2G has Glu138
and Met184.
