# Methods — Trajectory processing and structural features (extended)

The current section describes six analyses. The manuscript reports **eleven**
distinct structural quantities, and four of the ones it reports are not defined
anywhere: the DOR–RT centre-of-mass distance (Supplementary Figure 1B), the
NNIBP pocket volume (the G190E paragraph), the moiety-resolved contact counts
(the V106A paragraph) and the ring-centroid separation (the Tyr188 stacking
claim). Two definitions that *are* given are now wrong, and one describes an
analysis the paper never uses.

## What needs fixing in the existing text

| current text | problem |
|---|---|
| "within **4.5 Å**" | cutoff is **4.0 Å** — and Supplementary Figure 3A's caption already says 4.0 Å, so the two contradict each other |
| "the **number of protein heavy atoms** within…" | the quantity counts atom **pairs**; an RT atom close to three ring atoms contributes three |
| "The aryl ether torsion of DOR was measured as…" | this analysis appears nowhere in the Results — safe to delete, since you are trimming |

## Suggested replacement section

> **Trajectory processing and structural features.** Trajectory analyses were
> performed with MDAnalysis⁴³ and MDTraj.⁴⁴ Periodic-boundary artifacts were
> corrected by making protein chains whole and wrapping DOR to the nearest
> periodic image of the protein centre of mass. For structural comparisons,
> trajectories were aligned to the crystal structure reference using Cα atoms.
> All structural quantities were computed per frame and averaged within each
> replicate before being averaged across the three replicates, so that reported
> uncertainties reflect replicate-to-replicate variation rather than
> frame-to-frame fluctuation.
>
> DOR pose root mean squared deviation (RMSD) was computed after protein-pocket
> alignment from DOR heavy atoms relative to the crystallographic pose, and the
> DOR–RT centre-of-mass distance as the separation between the centre of mass of
> DOR and that of the protein. Residue–DOR proximity was computed as the minimum
> heavy-atom distance between DOR and each specified RT residue; named
> interactions between specific atoms, such as the hydrogen bond between the
> Lys103 main-chain carbonyl oxygen and the DOR triazolinone nitrogen, were
> measured as the distance between those two atoms directly.
>
> Contacts were counted as the number of protein–ligand heavy-atom **pairs**
> separated by less than 4.0 Å, so an RT atom lying close to several ligand atoms
> contributes once per pair. Burial of a ring was computed over the atoms of that
> ring, and total DOR burial over all DOR heavy atoms, in the same way. For the
> moiety-resolved analysis, DOR heavy atoms were partitioned by bond
> connectivity into its three ring systems — chlorocyanophenyl, pyridinone and
> triazolinone — each taking its own exocyclic substituents, with the shared
> ether oxygen and the methylene assigned to a linker group; contacts were then
> counted separately against each moiety.
>
> Ring planes were obtained by singular value decomposition of the ring atom
> coordinates, the interplanar angle between two rings was taken as the angle
> between their normals, and the ring-centroid separation as the distance between
> their mean atomic positions. The NNIBP volume was computed per frame on a
> 0.75 Å cubic grid spanning a 10 Å sphere centred on the Cα centroid of the
> sixteen pocket-lining residues (p66 L100, K101, K103, V106, T107, V108, V179,
> Y181, Y188, V189, G190, F227, W229, L234 and Y318; p51 E138⁶), counting grid
> points further from every protein heavy atom than that atom's van der Waals
> radius plus a 1.4 Å solvent probe. This definition requires no bound ligand and
> is therefore directly comparable between holo and apo simulations.

---

## Notes on choices

- **Replicate-before-frame averaging** is stated explicitly because it is what
  makes the reported ± values standard errors over three independent
  trajectories rather than over ~180 correlated frames, and because frame counts
  are not uniform across replicates (139–540 per replicate).
- **The pocket-volume residue list** cites ref 6, which is where it comes from
  (`src/analysis/metrics.py`); worth keeping the citation so the definition is
  not mistaken for an arbitrary choice.
- **The ligand-free pocket definition** matters more than it looks: it is what
  permits the apo/holo comparison, and it is why the G190E volume expansion can
  be attributed to the pocket rather than to the drug's displacement.
- The aryl ether torsion paragraph is dropped. The analysis exists
  (`compute_dor_torsions.py`) but no torsion is reported anywhere in the
  manuscript, so describing it in Methods invites a reviewer to ask where the
  result is.
