# JCIM ci-2026-030865 — Revision 1 triage plan

> **Status:** this is the original triage of the reviewer points, kept for the reasoning.
> The live task tracker is [Task-List.md](Task-List.md), which supersedes the effort
> estimates and the R1.2 scope below after the Sherlock data audit of 20-Oct-2026.

Deadline: **10-Nov-2026** (21 days from 20-Oct). Three reviewers, 12 scientific points
plus editorial/formatting items.

Categories: **[T]** text/Methods only · **[R]** reanalysis of data already on disk ·
**[S]** new simulation required.

---

## Reviewer 1

### R1.1 — Weak correlation may reflect mechanism, not just FEP error [T]
Discussion rewrite. Ties to R1.3. No new computation needed beyond what R1.2/R1.3 produce.

### R1.2 — Convergence within a basin ≠ adequate conformational sampling [S] — *the schedule driver*
Two requests:
- **(a) Independent starting poses.** Current replicates differ only by 0.1 Å coordinate jitter
  (Methods, p. 5) — reviewer is right that this samples one basin. Need genuinely
  independent holo starts for 2–3 representative systems (suggest WT, Y188L, G190E):
  re-dock / perturb the DOR pose, 3 × 100 ns each, then check the crystallographic
  pose is recovered (DOR RMSD, contacts).
  ≈ 9 new holo runs.
- **(b) Longer apo simulations for P225/M230 genotypes vs WT.** *Revised 20-Oct after
  auditing Sherlock.* Far more apo data exists than this plan originally assumed, split
  across two machines. At 100 ns × 3: A98G+F227C and F227C (Mac), **G190E (Sherlock only,
  since pulled)**. K103N+M230L is at 37/39/75 ns. WT, V106A, V106A+P225H and V106I+F227C
  are at 10 ns × 3. All three mismatch genotypes therefore have substantial apo data and
  **WT is the only real gap** — it is the reference, so its 10 ns caps every comparison.
  The twelve 300 ns runs this plan first proposed are not needed. See Task-List 2.7: do the
  matched 10 ns comparison first, and let it decide whether WT needs extending at all.

### R1.3 — Discuss the mismatches mechanistically [R]
K103N+M230L (small penalty, 36× resistant), A98G+F227C and G190E (favorable ΔΔG, resistant).
Feeds directly off R1.2(b): compare apo and holo remodeling for exactly these genotypes.
Trajectory data for holo is complete (3 × 100 ns for all 21 systems).

### R1.4 — Structural analysis needs a clearer summary [R]
- Figure mapping all mutations onto the RT–DOR structure relative to DOR/NNIBP — new figure,
  built from 4NCG; no simulation needed.
- WT-referenced heatmap of Table 3 metrics — `results/analysis/mechanisms/` +
  `paper/tables/Table-3-structural.csv` already hold everything; plotting work only.
- RMSF, secondary-structure changes, pocket conformational distributions — **the RMSF and
  DCCM code existed before the Sep-2026 repo trim** (deleted in `ce2dccb`, `967a600`,
  `bf3eb23`); restore from git history rather than rewriting. Runs on the 13 GB of local
  holo DCDs.

## Reviewer 2

### R2.1 — Convergence diagnostics (block averaging, autocorrelation) [R]
Block averaging already exists: `results/analysis/md_convergence/tables/block_averages.csv`
(`compute_md_convergence --n-blocks`). Missing: autocorrelation times / effective sample size,
and the same treatment applied to the MM/GBSA and FEP estimates. Extension of an existing CLI.

### R2.2 — Binary RT–DOR system omits nucleic acid and dNTP [T]
Honest limitations paragraph: primer grip (M230L) and the K103N allosteric/entry-gate argument
cannot be captured in a binary complex. Pairs naturally with R1.3 — it is the leading
explanation for the K103N+M230L false negative.

### R2.3 — Effect sizes with CIs; per-atom / torsion-level support [R]
- Bootstrap CIs on the key pairwise contrasts (V106A vs V106I/M, G190A vs G190E) from the
  per-replicate tables already in `results/analysis/mechanisms/` and Supplementary Table 4.
- Torsion analysis: the DOR torsion CLI (τ1–τ4) also predates the trim — restore from history
  and rerun on current trajectories. No PBC correction needed for intramolecular torsions.

### R2.4 — Charge-changing legs: protonation states, hybrid topologies, soft-core, finite-size [T, small R]
Mostly a Methods expansion, but with one substantive gap: `src/nnrti/fep/charge_correction.py`
implements only the periodicity self-energy term and its own docstring states it omits the
higher-order Rocklin terms (net-charge undersolvation, discrete-solvent). Either add those
terms or report their estimated magnitude explicitly. Also needs: explicit statement that
titratable NNIBP residues were assigned identically in WT and mutant, the pmx hybrid-residue
construction, and the Gapsys soft-core parameters actually used (α, σ, λ-schedule).

### R2.5 — Gasteiger vs AM1-BCC/RESP [S, small]
Gasteiger is the weakest methodological choice in the paper and the one a reviewer can most
easily insist on. Recommend doing the comparison rather than arguing: regenerate DOR parameters
with AM1-BCC (OpenFF default, `openmmforcefields` supports it directly), rerun MM/GBSA on
existing trajectories for a subset (cheap — no new MD), and ideally rerun FEP for 2–3 legs.
Note DOR is neutral, so the charge model affects ΔE_elec/ΔE_GB far more than ΔΔG_bind, which
is largely cancelled in the thermodynamic cycle — that argument plus a subset recomputation
should satisfy the point.

## Reviewer 3

### R3.1 — Define "correctly predicted" quantitatively [T]
A threshold is already implicit in `results/analysis/classification_performance/`
(TP/TN/FP/FN, MCC, AUC). State the ΔΔG_bind cutoff in the Abstract and Methods and make the
confusion matrix explicit.

### R3.2 — Typos [T]
"NRNTI" → "NNRTI" (p. 16); p. 30 wording; delete the redundant "APPROVAL" header (p. 35).
Also: the editor flagged English-language editing overall.

### R3.3 — Full minimization in MM/GBSA [R] — *confirmed, must act*
The published numbers used `snapshot_relaxation=unrestrained`, 2000 iterations
(`results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv`, all 60 rows), which is
exactly what the reviewer describes. **The code already supports `h_relax` and `none`**
(`src/nnrti/md/openmm/mmgbsa.py`, `--snapshot-relaxation`), so the sensitivity test is a flag
change and a CPU rerun, no new MD. Plan: rerun a subset (WT + ~5 genotypes, 18 replicates)
under all three settings, report ΔΔE_Total and components side by side, and justify the choice
on the documented H-atom/SHAKE clash artifact.

---

## Editorial / submission items
- TOC graphic (.tif or .pdf, last page, labelled "TOC Graphic", no caption).
- Supporting Information cover page (title, full author list and affiliations, contents list).
- Point-by-point response; clean manuscript; tracked/highlighted copy as "Supporting
  Information for Review Only"; clean SI for publication.
- ORCID iDs validated for all authors; funding (NIAID R24AI136618) in both the form and the file.
- Optional cover art + 80-word caption.

---

## Critical path and the one decision to make
Everything except **R1.2** is reanalysis, text, or a cheap rerun and fits comfortably inside
three weeks. R1.2 needs ~21 new MD runs on Sherlock plus their analysis, which is tight but
not impossible if the queue cooperates.

Options:
1. **Full scope** — do R1.2(a) and (b) as above; request a short deadline extension from
   Kirchmair now (extensions of 2–4 weeks are routinely granted and cost nothing to ask for).
2. **Reduced scope** — do R1.2(b) only for the four apo systems (the specific request), and
   answer R1.2(a) with a smaller 2-system pose test. Likely fits 10-Nov.
3. **Argue it** — decline new simulation and address R1.2 in Discussion only. Not recommended;
   this is Reviewer 1's central objection and the one they framed as conceptual.
