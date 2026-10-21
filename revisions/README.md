# JCIM ci-2026-030865 — Revision 1

Deadline **10-Nov-2026**. Everything for the resubmission lives here.

## Files to submit

| File | Upload as |
|---|---|
| `Manuscript-Revision-1-CLEAN.docx` | the manuscript (clean, unmarked) |
| `Manuscript-Revision-1.docx` | Supporting Information for Review Only (changes in red) |
| `Response-to-Reviewers.docx` | point-by-point response |
| `Supplementary-Text-Revision-1-CLEAN.docx` | Supporting Information for publication |
| `Supporting-Information-Cover-Page.docx` | SI cover page |
| `Supplementary-Table-5.xlsx` | new SI table (MM/GBSA relaxation sensitivity) |
| `TOC-Graphic.pdf` or `.tif` | TOC graphic (also already on the manuscript's last page) |

Supplementary Tables 1–4 are unchanged; take them from `paper/submission/`.

## Read these first

- **`Reviewer-Coverage-Audit.md`** — every reviewer point against what exists.
- **`Task-List.md`** — effort-ordered tracker with progress.
- `Revision-Plan.md` — the original triage, kept for the reasoning.

## Status

All 12 reviewer points are addressed and every `[PENDING]` marker in the
response is resolved. One number may still move slightly: the AM1-BCC
comparison currently rests on 8 of 9 replicates (Y188L has 2 of 3; the third is
computing). The conclusion does not depend on it.

## What needs you

1. **ORCID iDs** — each author validates their own in Paragon Plus.
2. **Funding** — NIAID R24AI136618 is in the manuscript; confirm it on the form.
3. ~~Zotero citations~~ — **done.** You inserted them in Word on 20-Oct; the
   manuscript now carries 96 Zotero fields (was 93) and they renumbered into ACS
   order as 46 (Kuhn 2005), 47 (Rastelli 2010) and 34,48 (Genheden & Ryde with
   Weis 2006). Nothing further needed.
4. **English-language pass** — the ~1,100 words added here were audited; the
   pre-existing text has not been. The editor raised it, so consider a service.
5. **Decide on R1.2a** — independent starting poses is the one scientific request
   not acted on. See the audit for the argument we make instead.

## What changed in the science

Most of the revision is new analysis on existing data. Five results changed what
the manuscript claims rather than only adding to it:

- The classification rule is `(estimate − SEM) > 0.5 kcal/mol`, not `> 0.5`.
  The Methods said the latter, which was wrong and looser than what was run.
- MM/GBSA ΔΔE is **materially sensitive** to the snapshot relaxation protocol
  (Supplementary Table 5). The manuscript now says so.
- **V106A vs V106M is within error on every metric**, and the V106A pyridinone
  contact change is not resolved at three replicates. The Discussion now rests
  the V106 argument on the Ser105 displacement, which is resolved (g ≈ 9).
- All three genotypes where ΔΔG_bind misses the phenotype have a larger
  induced-fit gap than WT, and M230L and A98G are the only substituted positions
  that cannot contact DOR. Reviewer 1's pre-binding hypothesis now has positive
  evidence, not just a plausible story.
- DOR stays within 0.64–0.93 Å of the crystallographic pose in every genotype.
  That is reported because it cuts against us: the complexes are stable, but the
  ensembles never visit an alternative pose, which is Reviewer 1's objection
  stated numerically.

## Reproducing

```bash
PYTHONPATH=src python ops/qa/verify_revision_numbers.py
```

Re-derives every number quoted in the documents from its source file. Run it
before submitting — the data files get regenerated and the prose does not.
