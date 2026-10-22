#!/usr/bin/env python3
"""Build the Supporting Information table workbooks.

Supplementary Tables 1-4 are XLSX, so Tables 5-7 are too rather than being
embedded in the Supplementary Text. Each workbook carries the summary the
text refers to, the per-genotype values behind it, and a Notes sheet stating
how the numbers were produced.

Usage
-----
    PYTHONPATH=src python ops/docx/build_si_workbooks.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revisions"
REL = ROOT / "results/analysis/relaxation_sensitivity/tables"


def _write(path: Path, sheets: dict[str, pd.DataFrame], notes: list[str]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name[:31], index=False)
        pd.DataFrame({"Notes": notes}).to_excel(
            xl, sheet_name="Notes", index=False, header=False)
    print(f"wrote {path.relative_to(ROOT)}  sheets: "
          f"{list(sheets) + ['Notes']}")


def table5() -> None:
    panel = pd.read_csv(REL / "relaxation_panel_long.csv")
    agree = pd.read_csv(REL / "relaxation_protocol_agreement.csv")
    pheno = pd.read_csv(REL / "relaxation_phenotype_correlation.csv")
    per_rep = pd.read_csv(
        ROOT / "results/.checkpoints/.checkpoint_mmgbsa_relaxsens_h_relax.csv")
    keep = ["structure", "mutation", "safe_label", "replicate", "fold_reduction",
            "binding_dg", "binding_dg_sem", "binding_dg_vdw", "binding_dg_vdw_sem",
            "binding_dg_electrostatic", "binding_dg_electrostatic_sem",
            "binding_dg_gb", "binding_dg_gb_sem", "mmgbsa_snapshot_relaxation"]
    _write(OUT / "Supplementary-Table-5.xlsx",
           {"Table 5": panel, "Table 5b agreement": agree,
            "Table 5b phenotype": pheno,
            "Per replicate h_relax": per_rep[[c for c in keep if c in per_rep]]},
           [
               "Supplementary Table 5. Sensitivity of MM/GBSA energies to the snapshot "
               "relaxation protocol.",
               "",
               "Sheet 'Table 5': WT-referenced ddE in kcal/mol for all 19 mutant "
               "genotypes, mean +/- SEM over three replicates, under two protocols.",
               "  h_relax    hydrogens relaxed, heavy atoms restrained (k = 10000 kJ/mol/nm^2)",
               "  full_min   all atoms free; the protocol used for Table 2 of the main text",
               "",
               "Sheet 'Table 5b agreement': Spearman correlation between the two protocols "
               "across the 19 genotypes, per energy component.",
               "Sheet 'Table 5b phenotype': each protocol against log10 fold reduction in "
               "susceptibility, over the 13 genotypes with an established phenotype.",
               "",
               "Sampling is identical across protocols and matches the main-text run: 100 "
               "snapshots spread evenly over the post-equilibration region of each 100 ns "
               "replicate, 0.25 discard fraction, 2000 iteration cap, CPU platform, double "
               "precision.",
               "",
               "Single-point evaluation, with no relaxation, was tested on six genotypes and "
               "is described in the Supplementary Text. It is not tabulated here: its mean is "
               "set by a small number of snapshots carrying steric contacts rather than by the "
               "protocol, so the values are not comparable with the two columns above.",
           ])


def table6() -> None:
    """Charge model: Gasteiger (used throughout) against AM1-BCC."""
    comp = pd.read_csv(REL / "charge_model_comparison.csv")
    charges = pd.read_csv(REL / "dor_charge_models.csv")
    per_rep = pd.read_csv(ROOT / "results/.checkpoints/.checkpoint_mmgbsa_am1bcc.csv")
    keep = ["structure", "mutation", "safe_label", "replicate", "fold_reduction",
            "binding_dg", "binding_dg_sem", "binding_dg_vdw", "binding_dg_vdw_sem",
            "binding_dg_electrostatic", "binding_dg_electrostatic_sem",
            "binding_dg_gb", "binding_dg_gb_sem", "mmgbsa_charge_method"]
    done = sorted(comp.genotype.unique())
    _write(OUT / "Supplementary-Table-6.xlsx",
           {"Table 6": comp,
            "DOR partial charges": charges,
            "Per replicate AM1-BCC": per_rep[[c for c in keep if c in per_rep]]},
           [
               "Supplementary Table 6. Sensitivity of the MM/GBSA energies to the ligand "
               "partial charge model.",
               "",
               "Sheet 'Table 6': WT-referenced ddE in kcal/mol, mean +/- SEM over three "
               "replicates, for each energy component under both charge models.",
               "  gast     Gasteiger charges, as used throughout the manuscript",
               "  am1bcc   AM1-BCC charges assigned by openff-toolkit",
               "The same trajectories and the same snapshots are scored under both models, "
               "so the difference isolates the charge assignment.",
               "",
               "Sheet 'DOR partial charges': per-atom charges under both models. AM1-BCC is "
               "the more polarised assignment (mean |q| 0.228 e against 0.138 e; RMS "
               "difference 0.165 e, maximum 0.390 e), with the largest differences on the "
               "triazolinone and pyridinone carbonyl and amide atoms.",
               "",
               "Sheet 'Per replicate AM1-BCC': the individual replicate values behind the "
               "means, alongside the wild-type reference.",
               "",
               f"Genotypes re-scored under AM1-BCC so far: {', '.join(done)}, each against "
               "the wild-type reference. Re-scoring is a per-replicate MM/GBSA calculation "
               "on the existing trajectories and does not require new molecular dynamics.",
           ])


def table7() -> None:
    """Unliganded pocket volume beside the bound-state volume."""
    apo = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_matched_window_summary.csv")
    holo = pd.read_csv(ROOT / "results/analysis/pocket_volume/tables/pocket_volume_genotype.csv")
    crystal = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_crystal_reference.csv")
    frames = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_frame_metrics.csv")

    apo = apo.assign(mutation=apo.safe_label.str.upper().str.replace("_", "+", regex=False))
    t = apo.merge(holo, on="mutation", how="left").rename(columns={
        "pocket_volume": "apo_volume", "pocket_volume_sd": "apo_volume_sd",
        "volume_mean": "holo_volume", "volume_sem": "holo_volume_sem"})
    t["induced_fit_opening"] = (t.holo_volume - t.apo_volume).round(1)
    t = t[["mutation", "n_rep", "apo_volume", "apo_volume_sd", "holo_volume",
           "holo_volume_sem", "induced_fit_opening", "chi1_Y181", "chi1_Y188",
           "dVolume_vs_WT"]].sort_values("mutation").reset_index(drop=True)
    t[["holo_volume", "holo_volume_sem"]] = t[["holo_volume", "holo_volume_sem"]].round(1)

    per_rep = (frames[frames.window == "matched_10ns"]
               .groupby(["safe_label", "replicate"], as_index=False)
               .agg(apo_volume_mean=("pocket_volume", "mean"),
                    apo_volume_std=("pocket_volume", "std"),
                    n_frames=("pocket_volume", "size"),
                    chi1_Y181=("chi1_Y181", "mean"),
                    chi1_Y188=("chi1_Y188", "mean")))

    _write(OUT / "Supplementary-Table-7.xlsx",
           {"Table 7": t, "Per replicate apo": per_rep, "Crystal reference": crystal},
           [
               "Supplementary Table 7. Unliganded (apo) and bound (holo) NNIBP volume.",
               "",
               "Volumes are the convex-hull proxy over the 16 NNIBP-lining Calpha atoms, "
               "in cubic angstroms.",
               "",
               "Sheet 'Table 7': per genotype, the apo volume over the matched window beside "
               "the bound-state volume of Table 3 of the main text.",
               "  apo_volume           mean over three replicates; sd is across frames",
               "  holo_volume          bound state, over full production",
               "  induced_fit_opening  holo minus apo: the opening DOR has to pay for",
               "  chi1_Y181, chi1_Y188 mean side-chain chi1 in the apo ensemble, degrees",
               "  dVolume_vs_WT        apo volume relative to wild-type apo",
               "",
               "The wild-type apo trajectories are 10 ns and set the comparison window, so "
               "every apo trajectory is cut to a matched 10 ns window. The bound-state "
               "volumes come from full 100 ns production, so the two columns are not matched "
               "in sampling; the apo column is the shorter of the two.",
               "",
               "Sheet 'Crystal reference': the same proxy evaluated on crystal structures. "
               "The unliganded crystals 1DLO and 1HMV bracket the simulated apo ensemble, "
               "and the DOR-bound 4NCG sits well above both.",
               "",
               "The apo pockets reported here were equilibrated without relocating bulk water "
               "into the cavity, which leaves the cavity under-solvated and biases the volume "
               "low. The induced-fit openings are therefore upper bounds.",
           ])


def main() -> int:
    table5()
    table6()
    table7()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
