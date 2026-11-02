#!/usr/bin/env python3
"""Move the structural results into Results, and say what the metrics mean.

Three things were wrong and are fixed here.

The Results subsection on the unliganded pocket spent a paragraph on pocket
volume, which is the least discriminating of the five metrics in Figure 4 --
three genotypes separate from WT and the errors are wide -- and said
nothing about the metrics that do separate. Figure 4's two most striking cells
had no text at all: A98G+F227C, whose pocket lining doravirine fails to
rigidify, and K103N+M230L, whose primer grip is the one place the panel shows a
mutation-induced change in the unliganded enzyme for that genotype. The
subsection is rewritten around what the figure actually resolves.

The V106 confidence intervals were introduced in the Discussion, which is for
interpreting results, not presenting them. They move to the Results paragraph
that reports the displacement.

Methods described the volume and the contacts but never defined the Calpha
RMSF, its superposition frame, the regions it is averaged over, the chi1
convention, or the induced-fit difference. A reader could not reproduce three
of the five columns of Figure 4.

Usage
-----
    PYTHONPATH=src python ops/docx/revise_results_section.py [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import docx_surgery as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MS = ROOT / "revisions/Manuscript-Revision-1.docx"
RESP = ROOT / "revisions/Response to Reviewers.docx"

HEADING = "Mutation-induced changes in RT conformation, bound and unliganded"

PARA_VOLUME = (
    "Simulating the unliganded enzyme for three 100 ns replicates per genotype lets "
    "every structural observable be referenced twice: against WT within one "
    "state (Figure 3), and as the change binding DOR produces in it, (bound − "
    "unliganded), referenced in turn to WT (Figure 4). The second is the "
    "quantity a thermodynamic cycle cannot see, because a mutation that moves both "
    "states by the same amount cancels in ∆∆Gbind by construction. The NNIBP is not "
    "an open cavity in unliganded RT, and WT pays 152 ± 15 Å³ of pocket volume "
    "to open it (Supplementary Table 7). Pocket volume is, however, the least "
    "discriminating of the five metrics in Figure 4: only three genotypes separate "
    "from WT, G190E and G190S opening the pocket 71 ± 26 Å³ and 65 ± 27 Å³ "
    "more on binding (95% intervals [10, 131] and [3, 127] Å³) and V106A+P225H 48 ± "
    "18 Å³ less ([−90, −6] Å³), while the other fifteen are within error. The V106A "
    "class pays least because its unliganded pocket is already wide, exceeding the "
    "wild-type apo volume by 40 to 52 Å³."
)

PARA_CHI1 = (
    "The side-chain dihedrals separate more genotypes. In WT, binding DOR "
    "rotates Tyr181 by only −18° and Tyr188 by +53°, from unliganded values of −68° "
    "and −97°. A98G+F227C and G190S depart sharply from that pattern: their "
    "unliganded Tyr181 sits close to WT's, at −68° and −52°, but binding "
    "drives it to −162° and −153°, a rotamer the wild-type complex does not adopt, "
    "giving induced-fit shifts of −76 ± 11° ([−101, −51]°) and −83 ± 28° ([−148, "
    "−18]°). The A98G+F227C value is the largest effect in Figure 4 relative to its "
    "error, at 6.9 standard errors. Y318F does the converse at Tyr188: its unliganded "
    "Tyr188 starts at +174°, a rotamer no other genotype visits, and DOR returns it "
    "to the common bound value of −38°, an induced-fit shift of +95 ± 32° ([22, "
    "169]°). Residue 318 does not line the pocket, so this is a distal substitution "
    "reorganising the unliganded tyrosine rather than the bound complex."
)

PARA_RMSF = (
    "The Cα root-mean-square fluctuation measures how much a region moves about its "
    "own mean position once rigid-body motion of the pocket has been removed, and it "
    "is the metric on which the two genotypes whose ∆∆Gbind misses the phenotype "
    "most clearly differ from WT. In WT, binding DOR damps both regions "
    "we track: the pocket lining falls from 0.452 ± 0.011 Å unliganded to 0.361 ± "
    "0.023 Å bound, and the primer grip from 0.570 ± 0.043 Å to 0.429 ± 0.021 Å. In "
    "A98G+F227C the pocket lining is not damped at all, moving from 0.394 ± 0.027 Å "
    "to 0.415 ± 0.009 Å, an induced-fit shift of +0.11 ± 0.04 Å from WT ([0.02, "
    "0.20] Å) and the second largest effect in Figure 4. Its bound interface is "
    "otherwise close to WT's, the only resolved contact metric being a −0.24 ± "
    "0.06 Å change in the Val179 distance, so DOR occupies this pocket without "
    "rigidifying its walls — a difference the bound-state interface energy does not "
    "register."
)

PARA_PRIMER = (
    "In K103N+M230L the effect is in the primer grip and runs the other way. Its "
    "unliganded primer grip is markedly more mobile than WT's, 0.772 ± 0.041 Å "
    "against 0.570 ± 0.043 Å, and binding DOR suppresses that excess, to 0.474 ± "
    "0.012 Å. The induced fit is therefore twice WT's, −0.30 Å against −0.14 Å, "
    "an induced-fit shift of −0.16 ± 0.06 Å ([−0.31, −0.01] Å). Residues 225–235 "
    "contact the primer terminus and lie outside the drug site, so this is an "
    "allosteric consequence of filling the pocket rather than a change in the "
    "interface: K103N+M230L has no computed binding penalty (∆∆Gbind = 0.59 ± 0.41 "
    "kcal/mol) and an unperturbed DOR interface, yet its unliganded enzyme differs "
    "from WT in this region more than any other genotype in the panel. V106I "
    "shows the same metric in the opposite direction, +0.14 ± 0.06 Å ([0.01, 0.27] Å), "
    "its primer grip being equally mobile with and without the drug. Across both "
    "states, 66 of 1,287 genotype-by-feature comparisons have a 95% interval "
    "excluding zero, close to the number expected if none differed, so the "
    "comparisons above are reported because they clear that threshold by a wide "
    "margin and not because they were selected after the fact (Supplementary Table 4)."
)

V106_STATS = (
    " Treating the replicate as the unit of replication, the displacement toward "
    "Ser105 is well resolved: −1.35 ± 0.09 Å for V106A against WT (95% "
    "interval [−1.61, −1.09] Å) and +1.85 ± 0.14 Å for V106I against V106A ([1.25, "
    "2.45] Å). The accompanying change in pyridinone contact count is not resolved "
    "(−3.2 ± 2.2 contacts for V106A against WT, [−9.3, +2.9]), and V106A "
    "against V106M is within error on every metric examined (Supplementary Table 4)."
)

METHODS_RMSF = (
    "Per-residue Cα root-mean-square fluctuation was computed after superposing every "
    "frame on the Cα atoms of the sixteen pocket-lining residues, so that the reported "
    "fluctuation is motion about each residue's own mean position with rigid-body "
    "translation and rotation of the pocket removed, and was then averaged over the "
    "residues of a region. Regions were the pocket lining defined above; the primer "
    "grip, residues 225–235; strand β6, residues 100–110; strands β9–β10, residues "
    "179–190; the YMDD motif, residues 183–186; and the p66 fingers (1–85, 118–155), "
    "palm (86–117, 156–236), thumb (237–318) and connection (319–426) subdomains. "
    "Side-chain χ1 for Tyr181 and Tyr188 was taken as the N–Cα–Cβ–Cγ dihedral, and "
    "both its mean and any difference between genotypes were computed on the circle, "
    "since a linear mean of angles lying either side of the branch cut is meaningless. "
    "Secondary structure was assigned per residue per frame by DSSP and reported as "
    "the fraction of a region's residues in helix or in strand. Every observable is "
    "reported in two comparisons: the shift of a genotype from WT within one "
    "state, and the induced fit, meaning the change that binding DOR produces in the "
    "observable, (bound − unliganded), itself referenced to WT. The second is a "
    "difference of four replicate means and is zero when a mutation moves the bound "
    "and unliganded states alike."
)

A_RESULTS_OLD = "The NNIBP is not an open cavity in unliganded RT, so DOR binding must pay a"
A_HEADING_OLD = "The unliganded pocket and the cost of opening it"
A_V106_RESULTS = "One of the most dramatic observations is in simulations of V106 substitutions"
A_V106_DISC = "The displacement toward Ser105 is well resolved"
A_METHODS_DUP = "Three unliganded replicates of 100 ns were run for each genotype"
A_METHODS_RING = "Ring planes were obtained by singular value decomposition"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ms, _ = S.load(MS)
    before = len(S.paragraphs(ms))

    # --- the Results subsection, rewritten around what Figure 4 resolves ----
    # The author removed the subsection heading in Word; the paragraphs stay
    # where they are rather than regaining one.
    body = S.find_para(ms, A_RESULTS_OLD)
    S.replace_span(body, S.para_text_raw(body), PARA_VOLUME, red=True)
    S.insert_after(S.parent_of(ms, body), body,
                   [S.build_paragraph([("text", t)], body)
                    for t in (PARA_CHI1, PARA_RMSF, PARA_PRIMER)])

    # --- the V106 intervals belong with the V106 result --------------------
    S.append_run(S.find_para(ms, A_V106_RESULTS), V106_STATS, red=True)
    disc = S.find_para(ms, A_V106_DISC)
    raw = S.para_text_raw(disc)
    S.delete_span(disc, raw[raw.index(" " + A_V106_DISC):])

    # --- Methods: define the metrics the figures are built from ------------
    ring = S.find_para(ms, A_METHODS_RING)
    S.insert_after(S.parent_of(ms, ring), ring,
                   [S.build_paragraph([("text", METHODS_RMSF)], ring)])
    # The apo volume protocol is stated in the structural-metrics paragraph
    # already; drop the copy that landed in the system-preparation paragraph.
    dup = S.find_para(ms, A_METHODS_DUP)
    raw = S.para_text_raw(dup)
    start = raw.index(" " + A_METHODS_DUP)
    keep = " Three unliganded replicates of 100 ns were run for each genotype, matching the bound-state protocol."
    S.replace_span(dup, raw[start:], keep, red=True)

    # --- Supplementary Table 10 is now tabs of Table 4 ---------------------
    for p in S.paragraphs(ms):
        if "Supplementary Table 10" in S.para_text(p):
            S.replace_span(p, "Supplementary Table 10", "Supplementary Table 4", red=True)
    si10 = [p for p in S.paragraphs(ms)
            if S.para_text(p).startswith("Supplementary Table 10:")]
    for p in si10:
        S.parent_of(ms, p).remove(p)

    after = len(S.paragraphs(ms))
    print(f"manuscript paragraphs {before} -> {after}")

    live = {
        "results": "“" + S.para_text(S.find_para(ms, "Simulating the unliganded enzyme for three 100 ns")) + "”",
        "v106": "“" + S.para_text(S.find_para(ms, A_V106_RESULTS)) + "”",
        "v106_disc": "“" + S.para_text(S.find_para(ms, "In all four V106A-containing genotypes, DOR slides")) + "”",
    }

    # --- Response: re-sync the quotes this touched -------------------------
    rp, _ = S.load(RESP)
    for anchor, text in (
        ("The NNIBP is not an open cavity in unliganded RT, so DOR binding must pay a", live["results"]),
        ("In all four V106A-containing genotypes, DOR slides", live["v106_disc"]),
    ):
        p = S.find_para(rp, anchor)
        S.replace_span(p, S.para_text_raw(p), text, red=True)
    for p in S.paragraphs(rp):
        if "Supplementary Table 10" in S.para_text(p):
            S.replace_span(p, "Supplementary Table 10", "Supplementary Table 4", red=True)

    if a.dry_run:
        print("dry run: nothing written")
        return 0
    S.save(MS, MS, ms)
    S.save(RESP, RESP, rp)
    print("wrote manuscript and response")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
