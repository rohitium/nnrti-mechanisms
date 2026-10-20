# JCIM ci-2026-030865 — Revision 1 task list

Deadline **10-Nov-2026**. Tasks are ordered by increasing effort within each tier.
Effort is focused working time, not wall clock. `R1.x` = Reviewer 1 point x.

**Progress: 6 / 23 complete** — last updated 20-Oct-2026

---

## Tier 0 — Editorial, no computation

- [ ] **0.1** Funding (NIAID R24AI136618) present in both the manuscript file and the Paragon submission form — already in the manuscript as of commit `3eed59d`, needs confirming on the form. _(Editor · 10 min)_
- [ ] **0.2** ORCID iDs validated and attached to each author's Paragon profile — each author must do their own. _(Editor · 15 min + co-authors)_
- [ ] **0.3** SI cover page: "Supporting Information" header, title, full author list and affiliations exactly as in the manuscript, then a contents list of the SI files. _(Editor #2 · 45 min)_

> **Dropped:** the "NRNTI" / "APPROVAL" / "hinges on" typos (R3.2). None of these strings
> exists in `DorDRM-MD-09-02-26-ACS.docx` or the SI, and the ACS-generated reviewer PDF is
> not retrievable from the portal. Answer R3.2 by stating that the manuscript has been
> proofread throughout and these corrections made.

## Tier 1 — Text-only scientific revisions

- [x] **1.1** Define "correctly predicted" quantitatively: state the ΔΔG_bind cutoff in Abstract, Results and Methods and give the explicit confusion matrix. Numbers exist in `results/analysis/classification_performance/classification_metrics.json` (TP 7, TN 3, FP 1, FN 2; MCC 0.50; AUC 0.78, p = 0.074). _(R3.1 · 2–3 h)_
- [ ] **1.2** Full proofread / English-language pass; state in the response that this was done. _(Editor, R3.2 · 1 day or outsourced)_
- [x] **1.3** Gasteiger discussion: DOR is neutral, so the charge model perturbs ΔE_elec and ΔE_GB far more than ΔΔG_bind, which largely cancels in the thermodynamic cycle. Write regardless of whether 4.1 is done — it frames that result. _(R2.5 part · 2 h)_
- [x] **1.4** FEP Methods expansion: titratable NNIBP residues assigned identically in WT and mutant; pmx hybrid-residue construction; Gapsys soft-core parameters (α, σ, λ-schedule) as actually run. _(R2.4 part · 3–4 h)_
- [x] **1.5** Limitations paragraph on the binary RT–DOR system: no nucleic acid, no incoming dNTP; name the foreclosed mechanisms (primer grip for M230L, the K103N closed-pocket/entry-rate argument) and connect them to the specific false negatives. _(R2.2, R1.1 · 3–4 h)_
- [ ] **1.6** Discussion rewrite on the mismatches (K103N+M230L, A98G+F227C, G190E). Draft as text first; strengthen once 2.7 lands. _(R1.3, R1.1 · 4–6 h)_

## Tier 2 — Reanalysis of data already on disk, no new MD

- [ ] **2.1** Mutation-map figure: all 18 genotypes on the RT–DOR structure, positions relative to DOR and the NNIBP. Built from 4NCG. _(R1.4 · 4–6 h)_
- [ ] **2.2** WT-referenced heatmap of the Table 3 structural metrics, from `paper/tables/Table-3-structural.csv` and `results/analysis/mechanisms/`. _(R1.4 · 3–4 h)_
- [ ] **2.3** Bootstrap confidence intervals on the key pairwise contrasts (V106A vs V106I/M, G190A vs G190E) from the per-replicate tables and Supplementary Table 4. Effect sizes with CIs, not means ± SEM. _(R2.3 part · 1 day)_
- [ ] **2.4** Autocorrelation times and effective sample size for the structural observables, and the same for the MM/GBSA and FEP estimates. Block averaging already exists in `results/analysis/md_convergence/tables/block_averages.csv`; this extends `compute_md_convergence`. _(R2.1 · 1–1.5 days)_
- [ ] **2.5** Restore the DOR torsion CLI (τ1–τ4) from git history and rerun on current holo trajectories. No PBC correction needed (intramolecular). _(R2.3 part · 1 day)_
- [ ] **2.6** Restore RMSF and DCCM code from history (`ce2dccb`, `967a600`, `bf3eb23`); rerun on local holo DCDs; add secondary-structure analysis and pocket conformational distributions. _(R1.4 · 1.5–2 days)_
- [ ] **2.7** Matched-window apo comparison, WT-referenced. Truncate all mutant apo trajectories to their first 10 ns to match the WT apo reference, then compare pocket volume, RMSF and loop conformations across G190E, A98G+F227C, F227C, K103N+M230L, V106A, V106A+P225H, V106I+F227C. Then ask separately whether the mutant observables drift between 10 ns and 100 ns — that determines whether a WT extension is needed and how long it must be. _(R1.2(b), R1.3 · 1.5 days)_

## Tier 3 — Reruns of existing pipelines with changed settings

- [ ] **3.1** MM/GBSA relaxation sensitivity. All 60 published rows used `snapshot_relaxation=unrestrained`, 2000 iterations — exactly what R3.3 flags. Rerun a subset (WT + ~5 genotypes = 18 replicates) under `none` and `h_relax`, both already implemented in `src/nnrti/md/openmm/mmgbsa.py`, and report ΔΔE_Total and components side by side. CPU only (OpenCL is unsafe here). _(R3.3 · 1 day setup + 2–3 days compute)_ — part (a), the written justification, is already in the revised Methods; this covers part (b), the numbers.
- [ ] **3.2** Complete the finite-size correction. `src/nnrti/fep/charge_correction.py` implements only the periodicity self-energy term and its docstring states the higher-order Rocklin terms are omitted. Either add net-charge undersolvation and discrete-solvent terms for the two charge-changing legs (K103N, G190E), or compute and report their magnitude as a bounded error. _(R2.4 substantive · 1.5–2 days)_ — the Methods now state plainly which terms were and were not applied; this task supplies the magnitude estimate.

## Tier 4 — New computation on existing trajectories

- [ ] **4.1** Regenerate DOR parameters with AM1-BCC (OpenFF default) and rerun MM/GBSA on existing trajectories for a subset. No new MD. _(R2.5 · 2 days)_
- [ ] **4.2** *Optional:* rerun 2–3 FEP legs with AM1-BCC charges to show ΔΔG_bind is insensitive to the charge model. Stronger answer, costs Sherlock time. _(R2.5 · 2–3 days + queue)_

## Tier 5 — New simulations — **deferred pending 2.7**

Decide only after the matched-window analysis says whether a WT extension is needed.

- [ ] **5.1** *Conditional:* extend WT apo from 10 ns to 100 ns, 3 replicates, restarting from `*_md_final.pdb`. ~0.7 GPU-day per replicate; ~1 day wall clock. Unlocks a full-length WT-referenced apo comparison against G190E, A98G+F227C, F227C and K103N+M230L. _(R1.2(b) · 1 day + analysis)_
- [ ] **5.2** *Conditional:* independent starting poses. Current replicates differ only by 0.1 Å jitter, so the reviewer's objection stands as written. Re-dock or substantially perturb the DOR pose for 2–3 representative systems (WT, Y188L, G190E), 3 × 100 ns each, and test whether the crystallographic pose is recovered. _(R1.2(a) · 1 day setup + ~3 days GPU + 1 day analysis)_

> **Not planned:** V106A+P225H apo extension. Dropped — it was carried from the reviewer's
> wording ("P225- and M230-containing genotypes") rather than from a gap in the evidence.
> K103N+M230L apo already exists at 37/39/75 ns. Revisit only if 2.7 shows the P225 class
> behaving differently.

## Tier 6 — Assembly and submission (last)

- [ ] **6.1** TOC graphic: .tif or .pdf, last page, labelled "TOC Graphic", no caption, legible text, original artwork. _(Editor · 4–6 h)_
- [ ] **6.2** Point-by-point response covering all 12 scientific points. _(All · 1–1.5 days)_
- [ ] **6.3** Clean manuscript + tracked/highlighted copy as "Supporting Information for Review Only" + clean SI for publication. Regenerate Tables 2–4 and the SI tables from whatever reruns land. _(Editor · 1 day)_
- [ ] **6.4** *Optional:* cover art + 80-word caption. _(Editor · 2–3 h)_

---

## Data preservation — done, outside the reviewer list

- [x] **D.1** Audit Sherlock for data not present on the Mac. Found: **G190E apo, 100 ns × 3** (`results/md_runs/apo/g190e/`), absent locally; FEP end-state equilibrium trajectories for the five uncleaned legs (5 ns each, 20 ns for the G190E leg). The Mac held 13 GB of a ~555 GB dataset.
- [x] **D.2** Pull G190E apo to the Mac — 3 × 200 frames, solvent-stripped, 112 MB, verified.
- [ ] **D.3** Tier-3 FEP archive to `$GROUP_HOME` before the `$SCRATCH` purge. FEP leg mtime is 2026-08-13, so purge falls around **11 Nov 2026** — one day after the revision deadline. Job `48676833` submitted 20-Oct; covers all 19 legs with endpoint ensembles and switch seeds but not the redundant `equil.trr`. Check with `bash ops/slurm/cluster/check_backup_status.sh`.

> A tier-1 backup from 2026-08-28 already covers four legs. It does **not** cover
> `wt_to_G190E` (158 GB, the largest, excluded at the time because it was still running)
> or the 14 small legs.
