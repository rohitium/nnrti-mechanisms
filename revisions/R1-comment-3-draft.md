# Reviewer 1, comment 3 — draft response and manuscript edits

Numbers below are from the apo panel as it stands today. They will be
re-stamped once all apo replicates reach 100 ns (4–6 days); the structure of
the argument does not depend on that.

---

## The response

**Reviewer 1, comment 3.** *The mismatches suggest different resistance
mechanisms. [...] K103N+M230L shows only a small calculated binding penalty
despite 36-fold resistance. A98G+F227C and G190E show substantial resistance
despite favorable calculated ΔΔGbind. G190E is particularly notable because it
produces the largest pocket expansion, yet no binding penalty. The authors
should examine whether mutation-induced changes occur in both apo and holo RT
and whether such remodeling could therefore be poorly reflected in ΔΔGbind.*

**Response.** We agree, and we have done the apo/holo comparison the reviewer
asks for. It supports the reviewer's hypothesis, with one qualification that
we think is important and that we now state in the paper: remodeling of the
unliganded pocket is larger than wild type for all three mismatched genotypes,
but it is not by itself what separates them from the rest of the panel.

To address this comment, we have added this paragraph to the **Methods**
section:

> **Unliganded RT simulations.** We simulated unliganded (apo) RT for wild type
> and for a subset of genotypes using the same protocol as the holo systems.
> Each apo system was built by removing DOR from the corresponding prepared
> holo structure, then solvated, minimized, heated and run under identical
> conditions. The NNIBP volume was computed per frame as described above. Apo
> replicates differ in length, so volumes were compared on a window common to
> all of them. Frames in which the pocket-lining Cα atoms were no longer a
> compact cluster, defined as a radius of gyration above 10 Å, were excluded
> before averaging; this removed 158 of 200 frames in one G190E replicate and
> 9 of 209 frames in one K103N+M230L replicate, in both cases because the
> subunit unfolded locally rather than because the pocket changed shape. We
> report the pocket-opening cost as the difference between the holo and apo
> NNIBP volumes of the same genotype.

To address this comment, we have added this subsection to the **Results**
section, after "Structural consequences of DOR-associated substitutions":

> **The unliganded pocket and the cost of opening it**
>
> The NNIBP does not exist as an open cavity in unliganded RT, so every
> genotype must pay some conformational cost to accommodate DOR. We measured
> that cost directly by running apo simulations and comparing the NNIBP volume
> with the holo value for the same genotype (Table 4).
>
> The unliganded pocket is smaller than the bound pocket in every genotype we
> simulated, by 117 to 242 Å³. Wild type pays 140 Å³. Three genotypes pay
> substantially more: G190E (242 Å³), V106I (196 Å³) and Y188L (195 Å³). Two
> pay less: V106A+P225H (117 Å³) and Y181C (120 Å³).
>
> All three genotypes whose ΔΔGbind misses the phenotype pay more than wild
> type to open the pocket: G190E by 103 Å³, A98G+F227C by 36 Å³ and
> K103N+M230L by 14 Å³. The cost is not specific to them, however. V106I
> (+56 Å³), Y188L (+55 Å³) and V106I+F227C (+39 Å³) also pay more than wild
> type, and for these three ΔΔGbind tracks the phenotype. What distinguishes
> the mismatched genotypes is that they pay this cost without losing contact
> with DOR once bound: G190E retains chlorocyanophenyl burial above wild type
> (22.7 ± 2.0 against 19.5 ± 0.5) and a stacking angle closer to coplanar
> (12.0 ± 0.9° against 13.4 ± 0.4°), whereas Y188L loses about a third of that
> burial (13.7 ± 1.5; Table 3).

> **Table 4.** NNIBP volume in the unliganded and DOR-bound states, and the
> difference between them. Apo values are means over three replicates on a
> matched window; holo values are reproduced from Table 3.
>
> | Genotype | Apo (Å³) | Holo (Å³) | Opening cost (Å³) | vs WT (Å³) |
> |---|---|---|---|---|
> | G190E | 42.0 | 284.3 | 242.3 | +102.7 |
> | V106I | 51.8 | 247.3 | 195.5 | +55.9 |
> | Y188L | 47.0 | 241.6 | 194.6 | +55.0 |
> | V106I+F227C | 54.4 | 232.8 | 178.4 | +38.8 |
> | A98G+F227C | 62.6 | 238.1 | 175.5 | +35.9 |
> | K103N+M230L | 97.0 | 250.5 | 153.5 | +13.9 |
> | V106A | 103.5 | 250.5 | 147.0 | +7.4 |
> | K103N+P225H | 100.0 | 239.7 | 139.7 | +0.1 |
> | WT | 88.8 | 228.4 | 139.6 | 0.0 |
> | V106A+L234I | 100.0 | 237.3 | 137.3 | −2.3 |
> | Y181C | 137.5 | 257.3 | 119.8 | −19.8 |
> | V106A+P225H | 106.7 | 223.3 | 116.6 | −23.0 |

To address this comment, we have modified this paragraph in the
**Discussion** section (replacing the text that previously reported these
numbers for the first time):

> Three genotypes show a mismatch between ΔΔGbind and the phenotype, and they
> do so for three different reasons.
>
> G190E and A98G+F227C distort the DOR–RT interface (ΔΔETotal = 1.94 ± 0.40 and
> 2.47 ± 0.24 kcal/mol), but they distort the unliganded pocket as well, and by
> more than wild type (Table 4). A mutation that perturbs both states in the
> same direction perturbs both legs of our thermodynamic cycle, so the penalty
> largely cancels and the residual favors binding (ΔΔGbind = −1.02 ± 0.38 and
> −1.53 ± 0.87 kcal/mol). G190E is the clearest case: it has the largest pocket
> expansion in the panel and the largest opening cost, yet it loses no aromatic
> packing, so little is left in the bound state for ΔΔGbind to register. A
> mutation can therefore raise the free-energy cost of forming the pocket
> without raising the free energy of binding to it, and our cycle is by
> construction least able to see exactly that.
>
> For G190E there is a second reason, which concerns the phenotype rather than
> the calculation. Substitutions at position 190 impair replication severely:
> the replication capacity of a site-directed G190E virus is too low to measure,
> and its NNRTI susceptibility could not be determined for that reason, while
> patient-derived G190E viruses replicate at 20–60% of reference because they
> carry compensatory substitutions elsewhere in RT and in Gag-PR, most often
> L74V or V75.76 The clinical fold-change we compare against is therefore
> measured on viruses that are not G190E alone, whereas our calculation changes
> one residue in an otherwise wild-type enzyme. Some of this mismatch is a
> difference in what is being measured, not an error in the calculation.
>
> K103N+M230L is different again. Its opening cost is only slightly above wild
> type, and the DOR interface is unperturbed: chlorocyanophenyl burial is the
> second highest in the panel (22.8 ± 2.4) and the Lys103 anchor is intact
> (3.05 ± 0.04 Å). Met230 lies in the p66 primer grip, the β12–β13 hairpin that
> positions the primer terminus at the polymerase active site; substitutions
> there alter template-primer binding and RNase H cleavage.79 M230L reduces
> NNRTI susceptibility and also impairs RT enzymatic function and replication
> capacity.80 Our system is a binary RT–DOR complex with no template/primer, so
> this mechanism is absent from the simulation rather than poorly sampled by it.
>
> A large opening cost is therefore necessary but not sufficient to produce a
> mismatch. V106I, Y188L and V106I+F227C also pay more than wild type to open
> the pocket, and for these ΔΔGbind tracks the phenotype, because each also
> loses contact with DOR in the bound state — most clearly Y188L, which loses
> about a third of its chlorocyanophenyl burial. ΔΔGbind misses a genotype when
> the pocket is remodeled and the bound-state contacts survive.

---

## Notes for the author

1. **This adds a Table 4** to the main text. The apo data currently sits only
   in Supplementary Table 7. Given that two reviewers ask about apo/holo
   remodeling, a main-text table seems right, but it is easy to move back.

2. **The structural fix matters more than the wording.** The apo numbers
   currently appear for the first time in the Discussion, with no Methods
   paragraph describing the apo simulations at all and no Results subsection.
   The edits above put the method in Methods and the numbers in Results.

3. **One existing claim should be dropped.** The current Discussion says the
   G190E apo expansion was traced to "reconstruction of periodic images for a
   chain-broken topology rather than a conformational event". That is not what
   the data shows: on the correctly imaged trajectory the pocket-lining Cα
   radius of gyration still reaches 42 Å in that replicate, and an independent
   test identified exactly two genuinely unfolded apo replicates in the whole
   set — this one and K103N+M230L rep 1. It is a local unfolding event, not an
   imaging artifact. The Methods text above states it that way, and that
   replicate is being re-run.

4. **Numbers to re-stamp when the apo runs finish** (4–6 days): every value in
   Table 4 and the three figures quoted in the Results text.

5. **The G190E replication-capacity point is new and I think it is the
   strongest part of the answer.** Huang et al. (ref 76, already cited) report
   that a site-directed G190E virus replicates too poorly to measure at all
   (<0.01% of reference) and that its NNRTI susceptibility could not be
   determined for that reason; the patient-derived G190E viruses that *can* be
   phenotyped replicate at 20–60% only because they carry compensatory changes
   such as L74V or V75, and in Gag-PR. So the 18-fold DOR reduction in Table 1
   was measured on viruses carrying more than G190E, while our calculation
   changes one residue in an otherwise wild-type background. This is a real
   limitation of the comparison, not of the simulation, and saying so directly
   is more convincing than attributing the whole gap to cycle cancellation.
   The current Discussion gestures at this ("may not be a binding effect") but
   does not make the compensatory-mutation argument.

6. **Ref 79 (Ghosh) supports the primer-grip claim as written** — the p66
   primer grip is the β12–β13 hairpin and alterations there change
   template-primer positioning and RNase H cleavage. Ref 80 (Xu) supports the
   M230L replication-capacity claim. Both were already cited; I have only made
   what they say more specific.
