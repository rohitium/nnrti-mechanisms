#!/usr/bin/env python3
"""Insert the promised Methods and Results text, then re-sync the Response quotes.

The Response document follows one rule: it quotes the manuscript verbatim. An
audit of every quoted block against the manuscript and the Supplementary Text
found four blocks containing sentences that appear in neither -- text drafted in
the Response as the answer to a reviewer and never inserted into the paper. A
reviewer checking the quotes against the manuscript would have found them
missing.

Two of those blocks are fixed by writing the paper: Methods gains the rest of
the unliganded-simulation protocol, and Results gains the subsection on the
cost of opening the pocket that the Response says was added. The other blocks
are fixed by regenerating the quote from the live manuscript, which is also
what keeps them correct after the Discussion reorganisation: a quote derived
from the source cannot drift from it.

Usage
-----
    PYTHONPATH=src python ops/docx/sync_response_quotes.py [--dry-run]
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
SI = ROOT / "revisions/Supplementary-Text-Revision-1.docx"

METHODS_APO = (
    " Three unliganded replicates of 100 ns were run for each genotype, matching the "
    "bound-state protocol. The NNIBP volume was computed per frame as above, and the "
    "pocket-opening cost reported as the holo minus the apo volume for the same "
    "genotype, with the standard errors of the two means combined in quadrature. Frames "
    "in which the pocket-lining Cα atoms were no longer a compact cluster, defined as a "
    "radius of gyration above 10 Å, were excluded; this affected one G190E replicate, in "
    "which the subunit unfolded locally, and is the only such exclusion in either state. "
    "The apo cavity was not re-solvated with bulk water before measurement, which biases "
    "apo volumes low, so the reported openings are upper bounds."
)

RESULTS_HEADING = "The unliganded pocket and the cost of opening it"
RESULTS_BODY = (
    "The NNIBP is not an open cavity in unliganded RT, so DOR binding must pay a "
    "conformational cost. We measured that cost by simulating the unliganded enzyme for "
    "three 100 ns replicates per genotype and subtracting the apo NNIBP volume from the "
    "holo value for the same genotype (Supplementary Table 7). The unliganded pocket is "
    "smaller than the bound pocket in every genotype, by 104 to 223 Å³, and the simulated "
    "wild-type apo volume, 76 ± 9 Å³, lies below the unliganded crystal structures 1DLO "
    "(102 Å³) and 1HMV (104 Å³), while the DOR-bound 4NCG value of 198 Å³ lies well above "
    "all three. Wild type pays 152 ± 15 Å³ to open the pocket. Of the three genotypes "
    "whose ∆∆Gbind misses the phenotype, only G190E pays significantly more, by 71 ± 26 "
    "Å³; A98G+F227C (22 ± 19 Å³) and K103N+M230L (2 ± 28 Å³) are within error of wild "
    "type. Nor is a large opening cost specific to the mismatched genotypes: G190S pays "
    "65 ± 27 Å³ more than wild type, second only to G190E and likewise resolved at three "
    "replicates, yet ∆∆Gbind tracks its phenotype, while V106I, Y188L and V106I+F227C pay "
    "40 ± 17, 35 ± 23 and 19 ± 25 Å³ more without the difference being resolved. In the "
    "opposite direction, every V106A-containing genotype pays less than wild type because "
    "its unliganded pocket is already wide: V106A exceeds the wild-type apo volume by "
    "52 ± 11 Å³, and V106A+L234I, V106A+P225H and V106A+F227L by 40 to 50 Å³. The V106A "
    "class therefore pre-organises the unliganded pocket toward the bound geometry rather "
    "than reshaping the interface only after DOR arrives."
)

G190E_PHENOTYPE = (
    " The comparison is further complicated by what can be measured: the clinical fold "
    "change for G190E is determined on patient-derived viruses carrying additional "
    "substitutions elsewhere in the genome, whereas our calculation changes one residue "
    "in an otherwise wild-type enzyme, so part of this mismatch is a difference in what "
    "is being compared."
)

A_METHODS_APO = "Apo systems were built by removing DOR from the corresponding prepared holo"
A_LAST_RESULT = "One of the most dramatic observations is in simulations of V106 substitutions"
A_CLOSING = "More generally, resistance can act through large-scale allosteric change"
A_SI9 = "The NNIBP is not an open cavity in the unliganded RT enzyme."


def quote(text: str) -> str:
    return "“" + text + "”"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    # ---------------------------------------------------------------- manuscript
    ms, _ = S.load(MS)
    before = len(S.paragraphs(ms))

    S.append_run(S.find_para(ms, A_METHODS_APO), METHODS_APO, red=True)

    last = S.find_para(ms, A_LAST_RESULT)
    heading_src = S.find_para(ms, "Structural consequences of DOR-associated substitutions")
    body_src = S.find_para(ms, A_CLOSING)
    S.insert_after(S.parent_of(ms, last), last, [
        S.build_paragraph([("text", RESULTS_HEADING)], heading_src),
        S.build_paragraph([("text", RESULTS_BODY)], body_src),
    ])

    closing = S.find_para(ms, A_CLOSING)
    raw = S.para_text_raw(closing)
    tail = "may not be a binding effect at all."
    S.replace_span(closing, tail, tail + G190E_PHENOTYPE, red=True)

    after = len(S.paragraphs(ms))
    print(f"manuscript paragraphs {before} -> {after}")

    # Pull the text the Response must quote, from the manuscript as it now stands.
    live = {
        "dde": S.para_text(S.find_para(ms, "Several genotypes show inconsistencies between")),
        "m230": S.para_text(S.find_para(ms, "Similarly, for K103N+M230L, our calculations")),
        "system": S.para_text(S.find_para(ms, "The second limitation is the model system")),
        "v106": S.para_text(S.find_para(ms, "In all four V106A-containing genotypes, DOR slides")),
        "closing": S.para_text(S.find_para(ms, A_CLOSING)),
        "methods_apo": S.para_text(S.find_para(ms, A_METHODS_APO)),
        "results_apo": S.para_text(S.find_para(ms, "The NNIBP is not an open cavity in unliganded RT")),
    }

    # ------------------------------------------------------- Supplementary Text
    si, _ = S.load(SI)
    par = S.find_para(si, A_SI9)
    S.replace_span(par, "V106I+F227C (171 ± 20 Å³) and A98G+F227C (175 ± 12 Å³)",
                   "K103N+M230L (154 ± 24 Å³) and A98G+F227C (175 ± 12 Å³)", red=True)

    # ------------------------------------------------------------------ Response
    rp, _ = S.load(RESP)

    def requote_para(p, text: str) -> None:
        S.replace_span(p, S.para_text_raw(p), quote(text), red=True)

    def requote(anchor: str, text: str) -> None:
        requote_para(S.find_para(rp, anchor), text)

    # The same Discussion paragraph is quoted under R1.3 and again under R2.2.
    # The two copies are byte-identical, so they are addressed by position, and
    # they now diverge: R1.3 keeps the mechanism argument, R2.2 keeps only the
    # point about which mechanisms a binary complex forecloses.
    dupes = S.find_paras(rp, "but comparable distortion is present in the apo simulations")
    if len(dupes) != 2:
        raise RuntimeError(f"expected that quote twice, found {len(dupes)}")
    r13, r22 = dupes

    requote_para(r13, live["dde"])
    S.insert_after(S.parent_of(rp, r13), r13,
                   [S.build_paragraph([("text", quote(live["m230"]))], r13),
                    S.build_paragraph([("text", quote(live["closing"]))], r13)])
    requote_para(r22, live["system"])
    # R2.3: the V106 paragraph as it now reads.
    requote("We quantified the uncertainty on these contrasts", live["v106"])

    p = S.find_para(rp, "Rather than treat the comparisons the reviewer names as special cases")
    S.replace_span(p, "in both the bound and unliganded states, using Welch intervals on the "
                      "difference of replicate means with the replicate as the unit of "
                      "replication (Supplementary Table 4)",
                   "in both the bound and unliganded states, using Welch intervals on the "
                   "difference of replicate means with the replicate as the unit of "
                   "replication (Supplementary Table 10)", red=True)

    requote("Apo systems were built by removing DOR from the corresponding prepared holo",
            live["methods_apo"])
    requote("The NNIBP is not an open cavity in unliganded RT, so DOR binding must pay a "
            "conformational cost. We measured that cost by simulating the unliganded enzyme and",
            live["results_apo"])

    if a.dry_run:
        print("dry run: nothing written")
        return 0
    S.save(MS, MS, ms)
    S.save(SI, SI, si)
    S.save(RESP, RESP, rp)
    print("wrote manuscript, supplementary text and response")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
