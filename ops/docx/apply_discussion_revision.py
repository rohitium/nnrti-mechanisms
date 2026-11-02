#!/usr/bin/env python3
"""Reorganise the Discussion and finish the figure renumbering.

What this does, and why each piece is here:

Discussion. The Discussion mixed mechanism with limitation, so the same caveat
had to be restated in three paragraphs. Each mechanism paragraph keeps its
mechanism and loses its caveat; the caveats are gathered into three paragraphs
that each state one thing once -- sampling, the model system, and the energy
functions. The third of those is new: nothing in the paper interpreted the
force field or the charge model as a limitation, although the charge-model
sensitivity had been measured.

The cancellation sentence is corrected rather than removed. It was right, but
it said "the apo simulations", which before this revision could only mean the
alchemical apo leg and now also names 3 x 100 ns equilibrium ensembles and a
figure built on holo - apo. It now names the leg and quantifies it.

Figures. The mechanism figure was renumbered to Figure 5 in Results but its
caption still read Figure 3, which the two new heatmaps also claimed, and the
Discussion still pointed at Figure 3A/3B for panels that are now Figure 5.
Both heatmap captions were one-line stubs.

Methods and the Supporting Information list. The only statement of how the
confidence intervals were computed was inside a Discussion passage being cut,
so it moves to the statistics paragraph of Methods. The Supporting Information
list still described the September submission: one Supplementary Text with two
figures and four tables, against nine sections, six figures and ten tables.

Every inserted run is red, per the revision convention. Citations are moved as
field elements, so a relocated sentence keeps its citation bound to the same
Zotero item; the two clones of existing citations get fresh citationIDs.

Usage
-----
    PYTHONPATH=src python ops/docx/apply_discussion_revision.py [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import docx_surgery as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MS = ROOT / "revisions/Manuscript-Revision-1.docx"

# --------------------------------------------------------------- anchors
A_G190 = "G190A resists first-generation NNRTIs"
A_V106 = "We quantified the uncertainty on these contrasts"
A_DDE = "Several genotypes show inconsistencies between"
A_M230 = "Similarly, for K103N+M230L, our calculations suggest negligible"
A_GENERAL = "More generally, resistance can act through large-scale allosteric change"
A_TYR = "Loss of aromatic stacking at Tyr181 and Tyr188"
A_FIG5 = "Trajectory-level evidence for two distinct resistance mechanisms"
A_FIG3 = "WT-referenced structural features related to RT-DOR contacts"
A_FIG4 = "RT Structural features extracted from equilibrium MD simulations"
A_STATS = "Associations between the computed energies and phenotype were assessed"
A_SITEXT = "Supplementary Text: Convergence of the DOR pose"
A_SIT4 = "Supplementary Table 4: Per-replicate values of the structural observables"
A_FF = "Protein parameters were assigned from the Amber14 OpenMM force-field files"
A_GB = "Here, ∆EvdW and ∆Eelec are the van der Waals and electrostatic"

# --------------------------------------------------------------- new prose
V106_TAIL = (
    "The displacement toward Ser105 is well resolved (−1.35 Å [−1.73, −0.97] for V106A "
    "against wild type, and +1.85 Å [1.25, 2.45] for V106I against V106A), whereas the "
    "accompanying change in pyridinone contact count is not (−3.2 contacts [−9.3, +2.9]), "
    "and the V106A against V106M contrast is within error on every metric we examined "
    "(Supplementary Table 10). We therefore rest the V106 argument on the displacement."
)

CANCEL_NEW = (
    "but the same perturbation is present in the unliganded leg of the thermodynamic "
    "cycle, so most of it cancels in the difference: the alchemical transformation to "
    "G190E costs −53.8 kcal/mol in the bound leg and −52.8 kcal/mol in the unliganded "
    "leg, leaving a residual of −1.02 kcal/mol (Supplementary Table 3). What ∆∆Gbind "
    "reports for these genotypes is a small remainder of a much larger perturbation, "
    "and here that remainder favors binding"
)

L_SAMPLING = (
    "Three limitations bound every interpretation above, and each applies to all of it "
    "rather than to one genotype. The first is sampling. Each system was simulated in "
    "triplicate from a single seeded pose, so the replicates explore one basin of the "
    "bound ensemble and their agreement measures convergence within that basin, not "
    "evidence that it is the only basin a resistant genotype occupies. The unliganded "
    "ensembles are the harder case, because the unoccupied pocket relaxes on timescales "
    "that 100 ns does not resolve. The consequence is visible across the panel as a "
    "whole: of the 1,287 genotype-by-feature structural comparisons we computed in the "
    "two states, 66 have a 95% confidence interval that excludes zero, close to what "
    "would be expected if none of them differed (Supplementary Table 10). Most "
    "structural features in this panel are therefore unresolved at three replicates, "
    "and we report them as unresolved rather than as unchanged. The comparisons we "
    "interpret above are those that clear that threshold by a wide margin."
)

L_SYSTEM_1 = (
    "The second limitation is the model system. We simulate a binary RT–DOR complex "
    "with no template/primer and no incoming dNTP, so a resistance mechanism acting "
    "through the nucleic-acid interface is structurally absent from our calculations. "
    "For this panel that is not a minor exclusion. Met230 and Pro225 both lie in the "
    "p66 primer grip, and a change in how the template is held, or in how often RT "
    "adopts a catalytically productive conformation, would not move ∆∆Gbind or "
    "∆∆ETotal as we measure them, yet could still lower susceptibility in cells. The "
    "same applies to the distortions we attribute to G190E, which might instead act on "
    "the primer grip,"
)
L_SYSTEM_2 = (
    " a region our system does not contain. What makes the binary model defensible for "
    "the drug site is that the pocket appears insensitive to the omission: in a cryo-EM "
    "structure of RT bound to doravirine and duplex DNA, the NNIBP Cα atoms and the DOR "
    "pose agree with 4NCG to 0.49 Å and 0.62 Å, while the fingers, the thumb and the "
    "p51 subunit move by 1.6–2.2 Å (Supplementary Section 1). The mechanisms we cannot "
    "reach lie outside the pocket rather than inside it."
)

L_FF_1 = (
    "The third limitation is in the energy functions themselves. Every number we report "
    "rests on fixed-charge, pairwise-additive force fields — ff14SB for the protein,"
)
L_FF_2 = (
    " OpenFF with Gasteiger partial charges for doravirine"
)
L_FF_3 = (
    " — and the MM/GBSA energies rest additionally on a generalised-Born treatment of "
    "solvent. Three consequences deserve stating plainly. First, the two energetic "
    "methods in this study do not share a ligand charge model: the equilibrium "
    "simulations, and the MM/GBSA energies taken from them, use Gasteiger charges, "
    "whereas the alchemical legs use the AM1-BCC charges OpenFF assigns by default. "
    "Re-scoring the same snapshots under AM1-BCC leaves every sign of ∆∆ETotal unchanged "
    "but moves its value by up to 1.33 kcal/mol (Supplementary Table 6), which is larger "
    "than every gap between adjacent genotypes in the ∆∆ETotal ranking, whose median is "
    "0.16 kcal/mol. The sign of an interface penalty is robust to the charge model; the "
    "rank of an individual genotype is not. Second, the shift falls almost entirely on "
    "the polar terms, and is largest where a substitution changes the formal charge. "
    "G190E is both cases at once: its ∆∆ETotal is dominated by the polar solvation term "
    "(∆∆EGB = 2.50 ± 0.32 kcal/mol), which an implicit solvent may overestimate for a "
    "glutamate buried in a hydrophobic pocket, and that same term moves by 1.91 kcal/mol "
    "under AM1-BCC charges. G190E is therefore the genotype whose energetics these force "
    "fields are least equipped to describe, and also one of the three whose ∆∆Gbind does "
    "not match its phenotype; we cannot separate those two facts with the methods used "
    "here. Third, these are systematic errors rather than statistical ones. They do not "
    "diminish with more replicates, and the confidence intervals quoted throughout this "
    "work do not contain them; MM/GBSA estimates are known to depend on the force field "
    "and the simulation protocol as much as on the structure being scored."
)
L_FF_4 = (
    " We have therefore drawn conclusions where a sign or a rank order carries them, and "
    "have avoided resting an argument on a difference of about one kcal/mol."
)

GENERAL_1 = (
    " G190 substitutions are a concrete instance: they increase susceptibility to DLV "
    "while severely impairing replication capacity, G190E most of all,"
)
GENERAL_2 = (
    " so the phenotype against which we score G190E may not be a binding effect at all."
)

METHODS_STATS = (
    " Differences between genotypes in a structural observable were quantified with the "
    "replicate as the unit of replication: the difference of replicate means, its "
    "standard error, a two-sided Welch interval on that difference, and Hedges' g with "
    "the small-sample correction. At three replicates per genotype a 95% interval "
    "requires 2.78 standard errors for a mutant-against-wild-type difference, and 2.31 "
    "for a difference of bound-minus-unliganded changes, which has more degrees of "
    "freedom; intervals are reported on that basis throughout (Supplementary Table 10)."
)

FIG3_CAPTION = (
    "Changes in RT–DOR contacts across the genotype panel, referenced to wild type. "
    "Each cell gives the shift of a structural metric from its wild-type value, in that "
    "metric's own units, above the standard error of the shift. Errors combine the "
    "genotype's and wild type's standard errors over three independent 100 ns replicates "
    "in quadrature, since wild type is itself a three-replicate mean. Colour is scaled "
    "within each column, so the darkest cell in a column is the largest shift for that "
    "metric; it carries no meaning across columns, which are in ångström, degrees and "
    "counts. Genotypes are grouped by the resistance category of Table 1 and ordered by "
    "fold reduction in doravirine susceptibility within each group. Metrics for which no "
    "genotype's 95% confidence interval excluded zero are not shown; on this criterion "
    "the Tyr188 / chlorocyanophenyl interplanar angle was removed. Of the 126 "
    "genotype–metric comparisons plotted, 21 have a 95% interval excluding zero. "
    "Per-comparison effect sizes, confidence intervals and p values for every metric, "
    "including those not plotted, are in Supplementary Table 10."
)

FIG4_CAPTION = (
    "Mutation-induced changes in RT's own conformation, as the difference between the "
    "bound and unliganded states, referenced to wild type. Each cell gives (holo − apo) "
    "for a genotype minus (holo − apo) for wild type: the extent to which binding "
    "doravirine reshapes the enzyme, and how far that reshaping differs from wild "
    "type's. A value near zero means the mutation moves the bound and unliganded states "
    "by the same amount, concordant remodelling, which cancels in a thermodynamic cycle "
    "and is therefore invisible to ∆∆Gbind by construction. Errors propagate the "
    "standard errors of all four three-replicate means in quadrature. Colour, grouping "
    "and ordering are as in Figure 3. Metrics for which no genotype's 95% confidence "
    "interval excluded zero are not shown; on this criterion the Cα RMSF of the YMDD "
    "motif, the p66 thumb and the p66 connection subdomains, and the β9–β10 strand "
    "content were all removed, and we report these as unresolved at three replicates "
    "rather than as unchanged. Of the 88 non-blank genotype–metric comparisons plotted, "
    "11 have a 95% interval excluding zero."
)

SI_TEXT = (
    "Supplementary Text: representativeness of the RT–DOR crystal structure and recovery "
    "of the crystallographic pose from perturbed starting poses (Supplementary Figures 1 "
    "and 2), convergence of the molecular dynamics ensembles (Supplementary Figure 3), "
    "non-equilibrium switching work distributions for every alchemical transformation "
    "(Supplementary Figure 4), finite-size corrections for the charge-changing alchemical "
    "legs, sensitivity of the MM/GBSA energies to the ligand charge model and to the "
    "snapshot relaxation protocol, protein-level conformational analysis (Supplementary "
    "Figure 5), and conformation of the unliganded binding pocket (Supplementary Figure "
    "6) (PDF)"
)

SI_TABLES = [
    "Supplementary Table 5: WT-referenced MM/GBSA energies under two snapshot relaxation "
    "protocols, with their agreement and their correlation with phenotype (XLSX)",
    "Supplementary Table 6: WT-referenced MM/GBSA energies under the Gasteiger and "
    "AM1-BCC ligand charge models, and the per-atom partial charges of each (XLSX)",
    "Supplementary Table 7: Unliganded and bound NNIBP volume for every genotype, with "
    "the induced-fit opening and crystal-structure reference values (XLSX)",
    "Supplementary Table 8: Comparison of the doravirine-bound NNIBP between PDB 4NCG "
    "and the independently determined structure 7Z2G (XLSX)",
    "Supplementary Table 9: Placement of Tyr181 and Tyr188 in RT bound to four NNRTIs "
    "(XLSX)",
    "Supplementary Table 10: Effect size and 95% confidence interval for every genotype "
    "and structural feature in both states, with the Welch statistics behind each (XLSX)",
]


def grab_citation(tree, anchor: str, shown: str, new_id: str):
    """Clone the field in `anchor` whose superscript reads `shown`."""
    par = S.find_para(tree, anchor)
    runs = [r for r in par if r.tag == S.W + "r"]
    groups, group, depth = [], [], 0
    import xml.etree.ElementTree as ET
    for run in runs:
        xml = ET.tostring(run, encoding="unicode")
        if 'w:fldCharType="begin"' in xml:
            depth += 1
        if depth:
            group.append(run)
        if 'w:fldCharType="end"' in xml:
            depth -= 1
            if depth == 0:
                groups.append(group)
                group = []
    for g in groups:
        if any(S.run_text(r) == shown for r in g):
            return S.new_citation(g, new_id, shown)
    raise LookupError(f"no citation reading {shown!r} in {anchor[:40]!r}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--doc", type=Path, default=MS)
    a = ap.parse_args()

    tree, _ = S.load(a.doc)
    before_paras = len(S.paragraphs(tree))
    before_fields = sum(S.field_count(p)[0] for p in S.paragraphs(tree))
    # The bibliography is one ZOTERO_BIBL field spanning every reference
    # paragraph, so two paragraphs are legitimately unbalanced on their own.
    # What matters is that the set does not grow.
    before_unbalanced = sum(1 for p in S.paragraphs(tree)
                            if S.field_count(p)[0] != S.field_count(p)[1])

    # ---- figure numbering, finished -------------------------------------
    S.replace_span(S.find_para(tree, A_FIG5),
                   "Figure 3: Trajectory-level", "Figure 5: Trajectory-level", red=True)
    S.replace_span(S.find_para(tree, A_TYR), "(Figure 3A)", "(Figure 5A)", red=True)
    S.replace_span(S.find_para(tree, "In all four V106A-containing genotypes, DOR slides"),
                   "(Figure 3B, Table 3)", "(Figure 5B, Table 3)", red=True)
    S.replace_span(S.find_para(tree, A_FIG3),
                   "WT-referenced structural features related to RT-DOR contacts extracted "
                   "from equilibrium MD simulations.", FIG3_CAPTION, red=True)
    S.replace_span(S.find_para(tree, A_FIG4),
                   "RT Structural features extracted from equilibrium MD simulations of the "
                   "RT-DOR system in holo and apo configurations.", FIG4_CAPTION, red=True)

    # ---- citations that have to travel ---------------------------------
    g190 = S.find_para(tree, A_G190)
    tail = S.para_text_raw(g190)
    cut = tail.index(" Notably, G190 substitutions increase susceptibility to DLV")
    moved = S.delete_span(g190, tail[cut:])
    if len(moved) != 2:
        raise RuntimeError(f"expected to recover 2 citations from the G190 paragraph, got {len(moved)}")
    cite_dlv, cite_grip = (S.redden_field(g) for g in moved)      # 76, 77

    cite_ff14sb = grab_citation(tree, A_FF, "37", "rsFFlim1")
    cite_gast = grab_citation(tree, A_FF, "40", "rsFFlim2")
    cite_mmgbsa = grab_citation(tree, A_GB, "34,48", "rsFFlim3")

    # ---- mechanism paragraphs lose their caveats -----------------------
    v106 = S.find_para(tree, A_V106)
    t = S.para_text_raw(v106)
    S.replace_span(v106, t[t.index("We quantified the uncertainty"):], V106_TAIL, red=True)

    dde = S.find_para(tree, A_DDE)
    S.replace_span(dde,
                   "but comparable distortion is present in the apo simulations, so the "
                   "penalty largely cancels in the thermodynamic cycle and the residual "
                   "favors binding", CANCEL_NEW, red=True)
    t = S.para_text_raw(dde)
    S.delete_span(dde, t[t.index(" In case of G190E,"):])

    m230 = S.find_para(tree, A_M230)
    t = S.para_text_raw(m230)
    S.delete_span(m230, t[t.index(" Our system is a binary RT–DOR complex containing no"):])

    # ---- the three limitation paragraphs -------------------------------
    new_paras = [
        S.build_paragraph([("text", L_SAMPLING)], g190),
        S.build_paragraph([("text", L_SYSTEM_1), ("field", cite_grip),
                           ("text", L_SYSTEM_2)], g190),
        S.build_paragraph([("text", L_FF_1), ("field", cite_ff14sb),
                           ("text", L_FF_2), ("field", cite_gast),
                           ("text", L_FF_3), ("field", cite_mmgbsa),
                           ("text", L_FF_4)], g190),
    ]
    S.insert_after(S.parent_of(tree, m230), m230, new_paras)

    # ---- the closing paragraph takes the replication-capacity point ----
    general = S.find_para(tree, A_GENERAL)
    t = S.para_text_raw(general)
    marker = "none of which are captured in these calculations."
    S.replace_span(general, marker, marker + "\u0000", red=False)
    # split the sentinel into the moved sentence plus its citation
    runs = [r for r in general if r.tag == S.W + "r"]
    host = next(r for r in runs if "\u0000" in S.run_text(r))
    S._set_text(host, S.run_text(host).replace("\u0000", ""))
    at = list(general).index(host)
    general.insert(at + 1, S.make_run(GENERAL_1, None, red=True))
    for offset, run in enumerate(cite_dlv, start=2):
        general.insert(at + offset, run)
    general.insert(at + 2 + len(cite_dlv), S.make_run(GENERAL_2, None, red=True))

    # ---- Methods: how the intervals were computed ----------------------
    S.append_run(S.find_para(tree, A_STATS), METHODS_STATS, red=True)

    # ---- Supporting Information list -----------------------------------
    sitext = S.find_para(tree, A_SITEXT)
    S.replace_span(sitext, S.para_text_raw(sitext), SI_TEXT, red=True)
    si_t4 = S.find_para(tree, A_SIT4)
    S.replace_span(si_t4, "summarized in Table 3 (XLSX)",
                   "summarized in Table 3, in both the bound and the unliganded state (XLSX)",
                   red=True)
    S.insert_after(S.parent_of(tree, si_t4), si_t4,
                   [S.build_paragraph([("text", line)], si_t4) for line in SI_TABLES])

    # ---- verification ---------------------------------------------------
    paras = S.paragraphs(tree)
    after_fields = sum(S.field_count(p)[0] for p in paras)
    unbalanced = sum(1 for p in paras if S.field_count(p)[0] != S.field_count(p)[1])
    print(f"paragraphs {before_paras} -> {len(paras)} (+{len(paras)-before_paras})")
    print(f"citation fields {before_fields} -> {after_fields} (+{after_fields-before_fields})")
    if unbalanced != before_unbalanced:
        raise RuntimeError(f"paragraphs with a field split across them: "
                           f"{before_unbalanced} -> {unbalanced}")
    print(f"no new split fields ({unbalanced} paragraphs share the bibliography field, "
          f"as before)")
    red = sum(ET.tostring(p, encoding="unicode").count('w:val="FF0000"') for p in paras)
    print(f"red runs: {red}")

    if a.dry_run:
        print("dry run: nothing written")
        return 0
    S.save(a.doc, a.doc, tree)
    print(f"wrote {a.doc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
