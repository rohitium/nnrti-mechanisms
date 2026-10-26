# Reviewer 1, comment 3 — draft response and manuscript edits

**Revised after re-reading the manuscript.** The cancellation argument the
reviewer is asking about is already in the Discussion in three places, and
Supplementary Table 7 already reports the holo-minus-apo opening under the name
`induced_fit_opening`. Nothing new needs to be built. Four things need doing:

1. extend Supplementary Table 7 from 8 genotypes to the 13 that now have apo data;
2. attach the literature to claims the paper already makes;
3. move the apo numbers out of the Discussion into Results, where they belong;
4. correct one statement the data no longer supports.

Numbers will be re-stamped when the apo runs reach 100 ns (4–6 days).

---

## What the manuscript already says

| Where | What it already argues |
|---|---|
| Discussion ¶593 | "G190E and A98G+F227C both distort the DOR–RT interface … but comparable distortion is present in the apo simulations, so the penalty largely cancels in the thermodynamic cycle and the residual favors binding" |
| Discussion ¶593 | K103N+M230L: interface unperturbed, Met230 in the primer grip, binary complex has no template/primer |
| Discussion ¶594 | The unliganded NNIBP is not an open cavity; a mutation affecting both states in the same direction gives ΔΔGbind ≈ 0 |
| Discussion ¶596 | Same cancellation argument, stated for the induced-fit pocket |
| Suppl. Table 7 | `induced_fit_opening` = holo − apo, per genotype, plus a crystal-reference sheet (1DLO/1HMV apo vs 4NCG holo) and the caveat that the apo cavity is under-solvated so openings are upper bounds |

So the response should point at this, not re-derive it.

---

## The response

**Response.** The thermodynamic-cycle cancellation the reviewer describes is
the explanation we had reached, and the apo simulations the reviewer asks for
were already the basis for it. We have extended that analysis to every genotype
for which we now hold apo trajectories, moved the numbers into the Results
where they belong, and added the literature that bears on it. The comparison
supports the reviewer's hypothesis, with one qualification we now state
explicitly: remodeling of the unliganded pocket is larger than wild type for
all three mismatched genotypes, but it does not by itself separate them from
the rest of the panel.

To address this comment, we have added this subsection to the **Results**
section, after "Structural consequences of DOR-associated substitutions":

> **The unliganded pocket and the cost of opening it**
>
> The NNIBP is not an open cavity in unliganded RT, so DOR binding must pay a
> conformational cost. We measured that cost by simulating the unliganded
> enzyme and subtracting the apo NNIBP volume from the holo value for the same
> genotype (Supplementary Table 7).
>
> The unliganded pocket is smaller than the bound pocket in every genotype, by
> 117 to 242 Å³, and the simulated wild-type apo volume (91 Å³) sits between
> the two unliganded crystal structures 1DLO (102 Å³) and 1HMV (104 Å³), while
> the DOR-bound 4NCG value (198 Å³) lies well above both. Wild type pays
> 137 Å³ to open the pocket. All three genotypes whose ΔΔGbind misses the
> phenotype pay more: G190E by 105 Å³, A98G+F227C by 38 Å³ and K103N+M230L by
> 20 Å³. The cost is not specific to them, however: V106I, Y188L and
> V106I+F227C also pay more than wild type, and for these three ΔΔGbind tracks
> the phenotype. What distinguishes the mismatched genotypes is that they pay
> this cost without losing contact with DOR once bound — G190E keeps
> chlorocyanophenyl burial above wild type (22.7 ± 2.0 against 19.5 ± 0.5) and
> a stacking angle closer to coplanar (12.0 ± 0.9°), whereas Y188L loses about
> a third of that burial (13.7 ± 1.5; Table 3).

To address this comment, we have added this paragraph to the **Methods**
section, after "Trajectory processing and structural features":

> **Unliganded RT simulations.** Unliganded (apo) systems were built by removing
> DOR from the corresponding prepared holo structure and were solvated,
> minimized, heated and run under the same protocol. The NNIBP volume was
> computed per frame as above, and the pocket-opening cost reported as the
> holo minus the apo volume for the same genotype. Apo replicates differ in
> length, so volumes were compared on a window common to all of them. Frames in
> which the pocket-lining Cα atoms were no longer a compact cluster, defined as
> a radius of gyration above 10 Å, were excluded; this affected one G190E
> replicate and one K103N+M230L replicate, in both cases because the subunit
> unfolded locally. The apo cavity was not re-solvated with bulk water before
> measurement, which biases apo volumes low, so the reported openings are upper
> bounds.

To address this comment, we have **modified the existing paragraph** in the
**Discussion** (¶593), keeping its argument and removing the numbers now given
in the Results:

> Three genotypes show a mismatch between ΔΔGbind and the phenotype, and they
> do so for different reasons. G190E and A98G+F227C distort the DOR–RT
> interface (ΔΔETotal = 1.94 ± 0.40 and 2.47 ± 0.24 kcal/mol), but they distort
> the unliganded pocket as well, and by more than wild type. A mutation that
> perturbs both states in the same direction perturbs both legs of the
> thermodynamic cycle, so the penalty largely cancels and the residual favors
> binding (ΔΔGbind = −1.02 ± 0.38 and −1.53 ± 0.87 kcal/mol). G190E is the
> clearest case: the largest pocket expansion in the panel and the largest
> opening cost, yet no loss of aromatic packing, so little is left in the bound
> state for ΔΔGbind to register. A mutation can raise the free-energy cost of
> forming the pocket without raising the free energy of binding to it, and the
> cycle is by construction least able to see that.
>
> For G190E a second factor concerns the phenotype rather than the calculation.
> Substitutions at position 190 impair replication severely: a site-directed
> G190E virus replicates too poorly for its NNRTI susceptibility to be
> determined at all, and patient-derived G190E viruses replicate at 20–60% of
> reference only because they carry compensatory substitutions elsewhere in RT
> and in Gag-PR, most often L74V or V75.⁷⁶ The clinical fold-change we compare
> against is therefore measured on viruses that are not G190E alone, while our
> calculation changes one residue in an otherwise wild-type enzyme. Part of
> this mismatch is a difference in what is being measured.
>
> K103N+M230L is different again. Its opening cost is close to wild type and the
> DOR interface is unperturbed: chlorocyanophenyl burial is the second highest
> in the panel (22.8 ± 2.4) and the Lys103 anchor is intact (3.05 ± 0.04 Å).
> Met230 lies in the p66 primer grip, the β12–β13 hairpin that positions the
> primer terminus at the polymerase active site; substitutions there alter
> template-primer binding and RNase H cleavage.⁷⁹ M230L reduces NNRTI
> susceptibility and also impairs RT function and replication capacity.⁸⁰ Our
> system is a binary RT–DOR complex with no template/primer, so this mechanism
> is absent from the simulation rather than poorly sampled by it.
>
> A large opening cost is therefore necessary but not sufficient to produce a
> mismatch. ΔΔGbind misses a genotype when the pocket is remodeled and the
> bound-state contacts survive.

---

## Notes for the author

1. **No new table.** Supplementary Table 7 already holds this. It needs
   extending from 8 genotypes to the 13 that now have apo trajectories, and
   the K103N+M230L row needs re-computing — its apo mean changes once the nine
   frames in which the subunit unfolds are removed.

2. **One existing claim should be dropped.** The Discussion currently says the
   G190E apo expansion was traced to "reconstruction of periodic images for a
   chain-broken topology rather than a conformational event". The data does not
   support that: on the correctly imaged trajectory the pocket-lining Cα radius
   of gyration still reaches 42 Å in that replicate, and an independent test
   found exactly two genuinely unfolded apo replicates in the whole set — this
   one and K103N+M230L rep 1. It is local unfolding. That replicate is being
   re-run.

3. **The G190E replication-capacity point is the one genuinely new thing here,
   and I think it is the strongest part of the answer.** Huang et al. (ref 76,
   already cited) report that a site-directed G190E virus replicates too poorly
   to phenotype (<0.01% of reference), that its susceptibility could not be
   determined for that reason, and that the patient-derived G190E viruses which
   can be phenotyped replicate at 20–60% because of compensatory changes. The
   current Discussion gestures at this ("may not be a binding effect") without
   making the argument.

4. **Supplementary Table 7's note about the matched window needs updating.** It
   says "The wild-type apo trajectories are 10 ns and set the comparison
   window". Wild type now has three replicates at 83–100 ns, and every apo run
   is heading to 100 ns, so that sentence expires in about a week.
