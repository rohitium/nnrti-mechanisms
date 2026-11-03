# JCIM ci-2026-030865 — Revision 1 status

"Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase".
Satija, Tao, Shafer. **Due 10-Nov-2026.** Last updated 02-Nov-2026, 16:00.

Written as a handoff. The author is away for several days; this is the only
live note in `revisions/`, and §1 is the part to read before touching anything.

---

## 1. Read this before editing any document

**The author edits the .docx files in Word, in parallel with scripted edits.**
This has already caused one near-miss: between two of my passes they rewrote
both figure captions, changed "wild type" to "WT" throughout, and deleted a
Results subsection heading. Before writing to any document:

```
git diff --quiet HEAD -- revisions/Manuscript-Revision-1.docx && echo clean || echo AUTHOR EDITED
```

If it reports an edit, read the file and build on it. Never regenerate a
`.docx`; edit in place. A `~$` file in `revisions/` means Word has that
document open right now — do not write to it.

**Tooling.** `ops/docx/docx_surgery.py` does run-aware surgery on a paragraph
and, unlike a naive text replacer, *moves Zotero citation fields* rather than
refusing to go near them. Its helpers: `find_para` / `find_paras` (address by
unique text, never by index — indices shift), `para_text_raw` (slice spans from
this, not `para_text`, which normalises the stray double spaces Word leaves),
`replace_span`, `delete_span` (returns the citation fields it removed, so a
relocated sentence keeps them), `build_paragraph`, `insert_after`. Three
applied scripts use it and are worth reading as examples:
`apply_discussion_revision.py`, `sync_response_quotes.py`,
`revise_results_section.py`.

**House rules.**
- Every change goes in **red** (`<w:color w:val="FF0000"/>`); the helpers do
  this by default. Supplementary Table .xlsx files are exempt.
- **Citations are live Zotero fields.** Never renumber by hand. The
  bibliography is a `ZOTERO_BIBL` field, so **the author must refresh Zotero in
  Word** — moved citations are out of numerical order until then.
- The Response quotes the manuscript **verbatim**. After any manuscript edit,
  re-run the quote audit in §3 and re-sync whatever drifted.
- Back up before a pass: `revisions/_backups/<date-time-label>/`.

## 2. Where everything stands

| document | state |
|---|---|
| `Manuscript-Revision-1.docx` | matches HEAD. Discussion reorganised, Results carry the structural findings, Methods complete |
| `Supplementary-Text-Revision-1.docx` | matches HEAD. Nine sections, six figures |
| `Response to Reviewers.docx` | matches HEAD. Every quoted sentence verified against the paper |
| `Supplementary-Table-1..9.xlsx` | complete set, all rebuilt or verified 02-Nov |
| `Supporting-Information-Cover-Page.docx`, `TOC-Graphic.*` | done |

Reviewers 2 and 3 are answered in full. Reviewer 1 is answered except as noted
in §6. The one scientific request not acted on is **R1.2a, independent starting
poses**: replicates differ only by 0.1 Å jitter, and the response says so
plainly. If there is a second round, that is where it comes from.

## 3. The two checks to run after any change

```bash
# 1. every Response quote must still appear in the manuscript or the SI
PYTHONPATH=src ~/miniconda3/envs/nnrti-prep/bin/python - <<'EOF'
import sys, re; sys.path.insert(0,'ops/docx')
import docx_surgery as S
from pathlib import Path
def alltext(p):
    x,_=S.load(Path(p)); return " ".join(S.para_text(e) for e in S.paragraphs(x))
MS=alltext('revisions/Manuscript-Revision-1.docx')
SI=alltext('revisions/Supplementary-Text-Revision-1.docx')
rp,_=S.load(Path('revisions/Response to Reviewers.docx'))
tot=miss=0
for p in S.paragraphs(rp):
    tx=S.para_text(p)
    if not (tx.startswith("\u201c") and tx.rstrip().endswith("\u201d")): continue
    q=tx.strip("\u201c\u201d").lstrip("\u2026").strip()
    s=[x for x in re.split(r'(?<=[.!?]) +', q) if len(x)>40]
    m=[x for x in s if x[:55] not in MS and x[:55] not in SI]
    tot+=len(s); miss+=len(m)
    for x in m: print("  MISSING:", x[:90])
print(f"{miss} of {tot} quoted sentences not in the paper")
EOF
```

It must print **0**. It was 19 when I first ran it, from four blocks drafted in
the Response and never inserted into the paper.

```bash
# 2. the figures and Supplementary Table 4 must agree
PYTHONPATH=src ~/miniconda3/envs/nnrti-prep/bin/python -m nnrti.cli.build_si_structural_tables
```

It prints a warning for any row where the p column and the panel criterion
disagree. It must print none.

## 4. The statistics — settled, do not reopen without reading this

**Pooled-variance (Student's) two-sided t test**, fixed degrees of freedom:
**4** for a genotype against WT, **8** for the induced fit, which is the
interaction term of a genotype × ligand-state design. Resolved means the shift
exceeds **2.776** or **2.306** standard errors respectively.

Welch was tried and rejected. With three replicates per group the standard
error is *identical* either way — pooling variances and propagating standard
errors are the same arithmetic at equal n — so Welch changes only the degrees
of freedom, and it estimates those from variances carrying 2 df each. The
estimate ranged 1.0–4.0 against a nominal 4, which cost real power (Figure 4
fell from 12 resolved cells to 5) and bought nothing: Welch earns its keep for
*unequal* group sizes, and these are equal by design. Hedges' g was removed at
the author's request; the effect size is the shift in the observable's own
units.

Counts under the adopted test: Figure 3, 56 of 527 resolved; Figure 4, 12 of
188; **68 of 715 overall against 36 expected by chance.**

**A claim may be discussed without being resolved** — the author's position,
and the right one — but only with a reason from outside this dataset: a
reviewer's question, prior literature, or a mechanism the model structurally
excludes. K103N+M230L's primer grip (p = 0.04, at the edge) is the worked
example; the text rests it on M230's location and the literature, not on the
margin. Do not go trawling the ~650 unresolved comparisons for stories.

## 5. Simulations, and the one job in flight

**I could not reach Sherlock at the end of this session — SSH returned
"Permission denied", so the credential has expired.** Re-authenticate before
assuming anything below is current.

Last known state (02-Nov 12:23): all 60 apo replicates are at 100 ns **except
v106a rep_02**, job **50786462**, which was at 92.6 ns and running on
sh02-13n13, having resumed correctly into the canonical tree
(`$SCRATCH/nnrti-mechanisms-git`). It has either finished or hit its wall clock
by now. Locally that replicate still reads 34.3 ns, so it is *excluded* from
every current number and V106A apo is n = 2.

### When it lands, regenerate in this order

```bash
SHERLOCK_USER=rsatija bash ops/sync/rsync_apo.sh pull      # trajectory + state
PYTHONPATH=src python -m nnrti.cli.fix_pbc_trajectories    # SEE WARNING BELOW
PYTHONPATH=src python -m nnrti.cli.compare_apo_ensembles
PYTHONPATH=src python -m nnrti.cli.compute_pocket_opening_cost
PYTHONPATH=src python -m nnrti.cli.compare_apo_to_crystal
PYTHONPATH=src python -m nnrti.cli.build_structural_metrics_table
PYTHONPATH=src python -m nnrti.cli.compute_structural_effect_sizes
PYTHONPATH=src python -m nnrti.cli.plot_structural_heatmap
PYTHONPATH=src python -m nnrti.cli.build_si_structural_tables
PYTHONPATH=src python ops/docx/build_si_workbooks.py
```

Use `~/miniconda3/envs/nnrti-prep/bin/python` — bare `python` lacks MDAnalysis
and MDTraj. **Analysis runs locally; Sherlock is for MD only.**

**`fix_pbc_trajectories` skips existing `*_analysis_pbcfix.dcd` files**, so the
extended v106a rep_02 will keep a stale 34 ns sidecar and nothing will warn
you. Compare its frame count against the raw `*_analysis.dcd` and pass
`--force` for that replicate only. DCDs are gitignored — back up before
forcing.

### Text that changes when V106A reaches n = 3

- Results, the pocket-volume paragraph: "V106A exceeds the WT apo volume by
  52 ± 11 Å³" and the V106A-class sentence.
- Supplementary Text section 9: "One replicate, V106A replicate 2, reached only
  34 ns and is excluded from the matched comparison" — delete it.
- Supplementary Table 7's V106A row (n goes 2 → 3) and Table 4 throughout.
- Re-run both checks in §3 afterwards.

The V106A-class conclusion itself does not depend on this replicate: three
combination genotypes carry it at n = 3 each.

## 6. Open work

- **Refresh Zotero in Word** and check the reference list. Nothing was
  renumbered by hand, but moved citations are out of order until a refresh.
- **Supplementary Text section 9 duplicates the Results.** It repeats the
  opening-cost numbers (152 ± 15, +71 ± 26, 104–223 Å³) the main text now
  states. Cut it back to what Results does not say.
- **"Supplementary Section 1" in the Discussion is the only section-number
  citation in the paper.** Everything else cites figures and tables. Change it
  to Supplementary Figure 1 or Table 8.
- **Methods lists three RMSF regions that appear in no figure or table** — β6,
  β9–β10 and the p66 palm — and of 18 helix/strand fractions only β9–β10 strand
  content is tabulated. Either trim the Methods list to what is reported, or
  add the three to `RT_METRICS` in `plot_structural_heatmap.py`, which puts
  them into Table 4 automatically. One line either way; the author had not
  chosen when the session ended.
- **The G190E replication-capacity specifics** the Response once asserted —
  20–60% of reference, compensatory L74V or V75 — are not in the manuscript. I
  could not reach Huang 2003 to verify them. If confirmed, they strengthen the
  G190E paragraph.
- **The primer-grip finding is not in the main text.** Primer-grip RMSF is
  resolved in 9 of 19 genotypes in the bound state, always in the same
  direction. It is the most systematic signal in the panel and currently
  appears only in the Supplementary Figure 5 caption.
- `mutation_component_summary.csv` is stale; use
  `binding_energy/tables/ddg_full.csv`.
- **Author actions:** ORCID iDs, the funding field (NIAID R24AI136618) on the
  submission form, and an English-language pass over the pre-existing text.

## 7. The numbers that matter

**Effect sizes, one method for every feature.**
The statistics live in the two tabs of Supplementary Table 4, built by
`nnrti.cli.build_si_structural_tables` from the same per-replicate data the
figures use, so a tab cannot drift from the figure it backs. See §4 for the
test. At n = 3, twice the standard error is **not** a 95% criterion — it is
p = 0.12.

- V106A toward Ser105: **−1.35 ± 0.09 Å [−1.61, −1.09]**, p = 0.004.
- K103N side chain to DOR: −3.31 ± 0.15 Å [−3.73, −2.89], p < 0.001.
- G190E clears 95% on **eight** holo features.
- V106M clears it on **none**, so "within error on every metric we examined" stands.
- 68 of 715 reach 95%, against 36 expected by chance.

**Apo pocket and opening cost** (3 × 100 ns, all 20 genotypes;
`results/analysis/apo_ensembles/tables/pocket_opening_cost.csv`). The holo arm
reproduces the published value to 0.4 Å³ (228.4 against 228), so the recompute
is trustworthy. Wild-type apo is **76.4 ± 9.4 Å³** and opening costs
**152.1 ± 15.0 Å³**, not 137. With errors propagated,
only **G190E (+70.7 ± 26.2 Å³)**, G190S (+65.1 ± 27.0) and V106A+P225H
(−47.9 ± 18.3) differ from wild type. A98G+F227C is +22.5 ± 18.9 and
K103N+M230L is +1.8 ± 28.4, so the claim that all three mismatched genotypes pay
more does not survive, and a large nominal cost is not peculiar to them: G190S
pays the second most and ΔΔGbind tracks its phenotype.

Three genotypes differ from wild type in apo volume, all V106A-containing and
all **more open**: V106A +52 ± 11 Å³ (p = 0.02), V106A+L234I, V106A+P225H. V106A
pre-organises the unliganded pocket toward the bound geometry, which is a
measured instance of what R1.2 asks about. Y181C is a milder second case
(apo 103.5 Å³, χ1 at +40° against WT's −68°).

One replicate is still excluded, G190E apo rep_02, which loses 158 of 200 frames
to the pocket-Cα Rg > 10 Å filter (local unfolding, not a pocket). That is now
the **only** exclusion; Methods still describes two. The apo cavity is not
re-solvated before measurement, so apo volumes are biased low and openings are
upper bounds. This is the one Methods caveat worth keeping. Delete the other,
"apo replicates differ in length", which is no longer true.

**The primer grip loosens, consistently.** Of the 171 genotype-by-region
fluctuation comparisons in the bound state, 17 are resolved at 95%, and they are
not scattered: the primer grip (225–235) becomes more mobile in **nine of the
nineteen genotypes**, by 0.08 to 0.16 Å, always in the same direction. G190E is
more mobile than wild type in four regions at once — pocket lining, primer grip,
β6 and the p66 fingers, the last by 0.56 Å [0.12, 0.99]. Secondary structure
moves by at most 0.11 of strand character and only 4 of 342 comparisons are
resolved. This is the most direct thing the data says about Reviewer 1's
question on the P225 and M230 genotypes, and a binary complex cannot connect it
to a mechanism. It is currently only in the Supplementary Figure 5 caption.

**DOR in water: binding stops the clock rather than narrowing the
distribution.** Across 60 bound replicates × 4 torsions × 100 ns, 24 µs of
cumulative torsion-time, DOR leaves its bound rotamer **three times**. Free in
water each torsion leaves several times per 100 ns at matched frame spacing, and
3–50 times at full resolution. Circular SD is 8–12° bound against 68–116° free.
The largest mutational effect on ligand rigidity is a quarter of the
bound-to-free gap, so **every genotype holds DOR in essentially the same
internal conformation** and resistance does not work by letting the ligand relax.
ΔG_preorg is 3–4 kcal/mol; do not quote a decimal, because τ3 crosses its
barrier only 3–6 times per 100 ns and the joint population multiplies through
it. ΔG_preorg is a free-ligand property, identical for every genotype, so it
cancels exactly in a mutant-minus-WT difference and cannot explain why ΔΔGbind
misses G190E.

**Charge model.** AM1-BCC now covers WT, V106A, Y188L, G190E and K103N+M230L at
3 replicates each (`charge_model_comparison.csv`, finished 02-Nov 02:59). The
largest ΔΔETotal shift is **1.33 kcal/mol (G190E)**, and all four signs are
preserved. The old "≤0.81 kcal/mol" in the coverage audit is superseded; the
documents already say 1.33.

## 8. Where the detail went

`_archive/` holds the eleven superseded memos. Two are still worth opening:

- `R1-crystal-representativeness.md` — the full 4NCG vs 7Z2G evidence and the
  manuscript-ready SI text. The text is already in the documents; the file is
  the record of how the numbers were obtained.
- `deleted-table-notes.md` — verbatim text deleted from the supplementary
  workbook Notes tabs, including the method and limitation lines that have
  nowhere to live except the Supplementary Text.

Technical traps that cost time live in `CLAUDE.md`, not here.
