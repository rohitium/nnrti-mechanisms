# Introduction mechanism sentences — drugs verified against the sources

Checked 2026-09-01 against the published abstracts and PDB depositions, not
recall. Two of my earlier guesses were wrong, in opposite directions.

| ref | paper | NNRTIs actually studied |
|---|---|---|
| 7 | Tantillo 1994, *J Mol Biol* 243:369–387 | survey; **nevirapine, TIBO, α-APA** |
| 8 | Schauer 2014, *NAR* 42:11687–11696 | **efavirenz, nevirapine, rilpivirine** |
| 9 | Hsiou 2001, *J Mol Biol* 309:437–445 | **loviride, HBY 097** |
| 10 | Lai 2016, *Viruses* 8(10):263 | **efavirenz, nevirapine, delavirdine** |
| 11 | Ren 2001, *J Mol Biol* 312:795–805 | **nevirapine, efavirenz, UC-781, TNK-651, PETT-2** |

## Findings that change the text

**Ref 8 is correct, and names efavirenz.** I had flagged this as a possible
miscitation; it is not. The abstract states that K103N "does not prevent binding
between **efavirenz** and RT–T/P but instead allows formation of a stable and
productive RT–T/P–dNTP complex, possibly through disruption of the E138–K101
salt bridge." The draft's claim is exactly supported, and the drug is efavirenz.

**Ref 9 did not use nevirapine or efavirenz.** Hsiou 2001 solved the unliganded
K103N RT structure and complexes with **loviride** and **HBY 097**, and proposed
that a hydrogen-bond network between Asn103 and Tyr188 stabilises the closed
(unliganded) pocket, slowing inhibitor entry.

**Ref 10 makes the Y181C stacking claim drug-dependent, which the draft states
as universal.** Lai 2016 reports Y181C resistance of **78-fold to nevirapine**
and **16-fold to delavirdine**, but only **2.7-fold to efavirenz**. The π–π
stacking loss therefore matters for the first-generation compounds and much less
for efavirenz. This is worth stating precisely, because it is the same argument
the paper later makes for DOR: an inhibitor that does not depend on Tyr181 is
largely unaffected by its loss.

**Ref 11 studied Y188C, not Y188L**, alongside Y181C — PDB 1JKH (Y181C·efavirenz),
1JLA (Y181C·TNK-651), 1JLB (Y181C·nevirapine), 1JLC (Y181C·PETT-2),
1JLF (Y188C·nevirapine), 1JLG (Y188C·UC-781). Its "second generation" means
efavirenz and UC-781, not the ETR/RPV/DOR generation the paper discusses later;
worth avoiding that collision.

**Ref 10's author list in the draft is correct** (Lai, Munshi, Lu, Feng,
Hrin-Solt, McKenna, Hazuda, Miller). Journal is *Viruses* 2016, 8(10), 263.

## Suggested rewrite

> NNRTI resistance mutations can decrease HIV-1 susceptibility by affecting
> interactions between the inhibitor and the NNIBP.⁷ K103N alters allosteric RT
> dynamics and permits formation of a stable, productive RT–template/primer–dNTP
> complex even while **efavirenz** remains bound,⁸ and independently stabilises
> the closed, unliganded form of the pocket through a hydrogen-bond network
> between Asn103 and Tyr188, slowing inhibitor entry — as shown structurally with
> **loviride and HBY 097**.⁹ Y181C removes an aromatic stacking interaction
> between the inhibitor and the NNIBP, conferring high-level resistance to
> **nevirapine** and **delavirdine** (78- and 16-fold) but little to
> **efavirenz** (2.7-fold), consistent with the differing dependence of these
> inhibitors on Tyr181.¹⁰,¹¹ G190A adds a bulge in a compact binding region,
> causing steric conflict with the inhibitor.⁶

The Y181C clause is the substantive improvement: as written the draft implies a
uniform mechanism, when the source shows a 30-fold spread across drugs — and
that spread is precisely the phenomenon the paper goes on to exploit for DOR.

## Sources

- Schauer 2014 — https://pmc.ncbi.nlm.nih.gov/articles/PMC4191400/
- Hsiou 2001 — https://pubmed.ncbi.nlm.nih.gov/11371163/
- Lai 2016 — https://pmc.ncbi.nlm.nih.gov/articles/PMC5086599/
- Ren 2001 — https://pubmed.ncbi.nlm.nih.gov/11575933/ ; PDB 1JKH, 1JLB, 1JLG
- Tantillo 1994 — https://pubmed.ncbi.nlm.nih.gov/7525966/
