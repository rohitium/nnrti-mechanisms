# Handoff: insert the 7Z2G comparison into the .docx files

Written 25-Oct-2026 for whichever agent holds the .docx files. I did the
analysis and drafted the text but **did not touch any .docx**, because you are
working in them.

Everything you need is in this file. Background, full numbers and reproduction
commands: `revisions/R1-crystal-representativeness.md`. Raw outputs:
`results/analysis/crystal_pose_comparison/`.

## What was done and why

Reviewer 1's general comment can be read as doubting that the 4NCG pose
represents the biologically relevant RT–DOR ensemble. A second experimental
RT–DOR structure exists and was not used in the paper: **7Z2G**, a 3.65 Å
single-particle cryo-EM structure of HIV-1 RT with a 38-mer DNA aptamer and
doravirine, from a different group. The two were compared.

Headline results, all reproducible from
`PYTHONPATH=src python -m nnrti.cli.compare_dor_crystal_poses`:

| | |
|---|---|
| DOR heavy-atom RMSD, fit on the 16 NNIBP Cα only | **0.62 Å** (29 atoms) |
| DOR centroid displacement in that frame | 0.30 Å |
| DOR internal conformation (ligand fitted on ligand) | 0.44 Å |
| NNIBP Cα RMSD | 0.49 Å |
| NNIBP side-chain RMSD | median 0.80 Å, max 2.34 Å (Lys101) |
| Residue-to-DOR minimum distances, **no superposition** | median \|Δ\| 0.16 Å, max 0.41 Å, r = 0.974 |
| Largest difference in the four rotatable dihedrals | 14.8° (τ4) |
| NNIBP volume, same routine as the trajectories | 198.3 Å³ (4NCG) vs 194.1 Å³ (7Z2G) |
| Cα RMSD by subdomain, palm frame | p51 2.23, fingers 1.83, thumb 1.62, RNase H 1.37, palm 1.21, connection 0.77, **NNIBP lining 0.97** |

The last row is the rhetorical centre: the two structures differ in the
subdomains that grip nucleic acid, and agree at the drug site.

## Assets, already built

| what | where |
|---|---|
| Supplementary Figure 5 | `results/analysis/crystal_pose_comparison/plots/crystal_pose_comparison.png` (and `.pdf`, 300 dpi) |
| optional panel E | `results/analysis/crystal_pose_comparison/plots/nnibp_overlay.png` |
| Supplementary Table 8 | `revisions/Supplementary-Table-8.xlsx`, six sheets |
| CSV mirror, matching the other SI tables | `revisions/si_tables/SI-crystal_pose.csv` |

## Numbering

Checked against the current SI: Supplementary **Figures run to 4**,
Supplementary **Tables run to 7** (main-text tables take 1–4), SI **sections
run to 7**. So this becomes **SI section 8**, **Supplementary Figure 5**,
**Supplementary Table 8**. If you have added anything since, shift accordingly
and tell me so I can fix the audit.

---

# The insertions

Four files. Formatting convention already in use: additions in **red** in the
marked copies, plain text in the `-CLEAN` copies, identical wording in both.

## 1. `Supplementary-Text-Revision-1.docx` (and `-CLEAN`)

### 1a. New section at the very end

Section 7 ("Conformation of the unliganded binding pocket") is currently last.
Append after its Supplementary Figure 4 caption:

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

### 1b. Figure caption, under the inserted figure

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

If you also place the structural overlay as panel E, append:

> (E) The two ligands after superposition on the NNIBP Cα atoms alone (4NCG
> blue, 7Z2G red), with the NNIBP side chains of 4NCG in grey. Heavy-atom RMSD
> 0.62 Å over 29 atoms.

### 1c. Contents list, first page

Add a line matching the existing dot-leader format:

> 8  An independent structure of the RT–DOR complex .......... S10

Check the page number against the built document; S10 assumes section 7 ends on
S9, which is what the current contents list implies.

## 2. `Manuscript-Revision-1.docx` (and `-CLEAN`)

### 2a. Discussion

The sampling-limitation paragraph currently ends:

> "…so the complexes are stable but the ensembles do not visit alternative
> poses."

Append immediately after that sentence, in the same paragraph:

> That pose is not peculiar to the structure we started from: doravirine and
> the residues lining its pocket agree to 0.62 and 0.49 Å respectively between
> 4NCG and 7Z2G, an independent cryo-EM structure of RT with doravirine and
> DNA bound (Supplementary Figure 5), although the sampling limitation above
> applies to the basin around that pose regardless.

### 2b. Discussion, optional second sentence

Only if the Discussion makes the nucleic-acid point for Reviewer 2's comment 2.
It can follow the sentence above or sit in the nucleic-acid paragraph:

> The same comparison localises what nucleic acid changes: with the p66 palm as
> the reference frame, p51, the fingers and the thumb differ by 1.6–2.2 Å
> between the two structures while the NNIBP lining differs by 0.97 Å, so the
> mechanisms we cannot reach without nucleic acid are those that act through
> those subdomains rather than through the drug site itself.

### 2c. Methods

The structure-preparation paragraph begins "The starting holo structure was the
crystal structure of HIV-1 RT bound to doravirine (DOR; PDB ligand residue 2KW)
from PDB 4NCG…". Add at the end of that paragraph, or as a short paragraph of
its own:

> Experimental structures were compared with gemmi; superpositions are Kabsch
> fits on the stated Cα sets with no outlier rejection. Ligand atoms were
> matched by chemical-component name, minimised over the automorphisms of the
> ligand bond graph, and side-chain RMSDs were minimised over the equivalent
> terminal atoms of symmetric residues. Pocket volumes for both structures were
> measured with the routine used on the trajectories; the value for 4NCG,
> 198.3 Å³, is the one already reported for that structure.

### 2d. Supporting Information contents, near the end of the manuscript

Add a Supplementary Table 8 line in the style of the Supplementary Table 1–4
lines:

> Supplementary Table 8: Agreement between the two experimental RT–DOR
> structures: summary quantities, the per-residue comparison of the 16
> NNIBP-lining residues, Cα RMSD by subdomain, ligand RMSD by moiety, the
> rotatable dihedrals and the NNIBP gate distances (XLSX)

**Pre-existing problem, please check rather than copy me.** That list names
Supplementary Tables 1–4 and, in the Supplementary Text line, only
Supplementary Figures 1 and 2. Supplementary Tables 5–7 and Supplementary
Figures 3–4 exist but are not listed in either manuscript copy. Either the list
is out of date or the convention is deliberate. Whatever is decided for 5–7
applies to 8, and the Supplementary Text line needs Figures 3, 4 and 5 adding
if the list is being brought up to date.

## 3. `Response-to-Reviewers.docx`

Under Reviewer 1's general comment. Full drafted response is in
`revisions/R1-crystal-representativeness.md` under "Draft text for the
response" — it is a block quote of five paragraphs, written to be pasted as is.
It ends by stating plainly that this does **not** address R1.2a (independent
starting poses), which should stay in.

If Reviewer 2's comment 2 is answered separately, the subdomain result is worth
cross-referencing there too.

---

# The reference problem — decide before you insert

The reference list is **88 manually numbered entries**, "(1)" … "(88)", with
plain-text superscript citations. There is no field-based reference manager to
renumber for you.

The new citation is:

> Singh, A. K.; De Wijngaert, B.; Bijnens, M.; Uyttersprot, K.; Nguyen, H.;
> Martinez, S. E.; Schols, D.; Herdewijn, P.; Pannecouque, C.; Arnold, E.;
> Das, K. Cryo-EM structures of wild-type and E138K/M184I mutant HIV-1 RT/DNA
> complexed with inhibitors doravirine and rilpivirine. *Proc. Natl. Acad. Sci.
> U.S.A.* **2022**, *119*, e2203660119. DOI: 10.1073/pnas.2203660119.

ACS numbers references in order of first appearance. The first mention of 7Z2G
would be in Methods, immediately after the 4NCG citation, which is reference
14 — so strictly the new entry becomes **(15)** and references 15–88 shift to
16–89, along with every in-text superscript. Two options:

- **Renumber.** Correct, mechanical, and best done with a script over
  `word/document.xml` rather than by hand: highest-first substitution so (15)
  does not collide with (16). Verify the count is 89 afterwards and that no
  superscript was missed.
- **Append as (89).** One edit, no risk to the other 88, but out of sequence
  and a copy-editor will catch it.

**This is the author's call, not ours.** Ask before doing either.

## Two things worth saying in the response

- 7Z2G is the **wild-type** entry from that paper; the same paper also reports
  an E138K/M184I mutant RT–DOR structure. I verified it: 7Z2G has Glu138 and
  Met184. A reviewer who knows the paper will check, so say which one was used.
- The 198.3 Å³ measured here for 4NCG is **exactly** the value already in the
  SI for that structure (`apo_crystal_reference.csv`), so the new pocket volume
  contradicts nothing already written.

## After you insert

- Marked and CLEAN copies must match word for word; that is checked elsewhere
  in this revision and should stay true.
- `revisions/Reviewer-Coverage-Audit.md` already has the 7Z2G row under
  Reviewer 1 and a strengthened R2.2 row. Item 4 of "What remains" says this
  insertion is outstanding — **delete that item once it is done.**
- Supplementary Figure 5 should be placed at 300 dpi; the `.pdf` is vector and
  preferable if the SI accepts it.
