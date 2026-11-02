#!/usr/bin/env python3
"""Build the Supporting Information table workbooks.

Supplementary Tables 1-4 are XLSX, so Tables 5-7 are too rather than being
embedded in the Supplementary Text. Each workbook carries the summary the text
refers to and the per-genotype values behind it.

There is deliberately no Notes sheet. Every column that has a unit says so in
its name (``_kcal``, ``_A3``, ``_deg``, ``_e``), so a number can be read where
it sits. What a Notes sheet cannot express -- why a protocol was chosen, what a
measurement does not account for -- is prose and belongs in the Supplementary
Text, not in a tab beside the data.

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

#: Columns holding a physical quantity must end in one of these, so the unit
#: travels with the number instead of living in a note.
_UNIT_SUFFIXES = ("_kcal", "_A3", "_A", "_deg", "_e", "_p", "_ns")

#: Counts, labels and dimensionless ratios, which need no unit.
_UNITLESS = {
    "genotype", "energy_component", "relaxation_protocol", "replicate",
    "n_replicates", "n_genotypes", "n_frames", "n_ca_atoms", "element",
    "atom_index", "source", "state", "sign_agreement", "dor_fold_reduction",
    "R2_vs_log10_fold", "spearman_h_relax_vs_full_min", "spearman_vs_log10_fold",
}

#: Energy columns in the MM/GBSA checkpoints and the component stem each maps
#: to. The stems are the same on every sheet, so moving between a summary and
#: its per-replicate data needs no translation.
_ENERGY_COLUMNS = {
    "binding_dg": "total",
    "binding_dg_vdw": "vdw",
    "binding_dg_electrostatic": "elec",
    "binding_dg_gb": "gb",
}

#: Per-replicate numbers for each MM/GBSA variant.
_CHECKPOINTS = {
    "full_min": "results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv",
    "h_relax": "results/.checkpoints/.checkpoint_mmgbsa_relaxsens_h_relax.csv",
    "gasteiger": "results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv",
    "am1bcc": "results/.checkpoints/.checkpoint_mmgbsa_am1bcc_extended.csv",
}

#: kJ/mol in the checkpoints; every sheet reports kcal/mol.
_KJ_PER_KCAL = 4.184


def _check_columns_self_describing(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    """Fail the build rather than ship a column whose unit is unstated.

    With no Notes sheet there is nowhere else a unit could be written down, so
    this is the only thing between a reader and a number they cannot interpret.
    """
    bad = [f"{sheet}:{c}" for sheet, df in sheets.items() for c in map(str, df.columns)
           if c not in _UNITLESS and not c.endswith(_UNIT_SUFFIXES)]
    if bad:
        raise AssertionError(
            f"{path.name}: these columns state neither a unit nor a known unitless "
            f"quantity, and there is no Notes sheet to explain them: {', '.join(bad)}")


def _write(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    _check_columns_self_describing(path, sheets)
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name[:31], index=False)
    print(f"wrote {path.relative_to(ROOT)}  sheets: {list(sheets)}")


def _display(label: str) -> str:
    """Genotype label as the tables print it: WT, V106A, K103N+M230L.

    A rule rather than a lookup, so a genotype added later is formatted without
    anyone remembering to extend a dictionary.
    """
    if label.lower() == "wt":
        return "WT"
    return label.upper().replace("_", "+")


def _stem(component: pd.Series) -> pd.Series:
    """ddE_Total -> Total. The column name already says these are ddE values."""
    return component.astype(str).str.replace("ddE_", "", regex=False)


def _per_replicate_pair(sources: dict[str, str], genotypes: list[str]) -> pd.DataFrame:
    """One row per genotype and replicate, two variants side by side.

    Checkpoint energies are kJ/mol and every sheet here is kcal/mol, so they
    are converted. The summary sheets were always kcal/mol while the
    per-replicate sheets used to be dumped raw, leaving two sheets of one
    workbook differing by a factor of 4.184 with nothing saying so.
    """
    names = list(sources)
    wanted = {_display(g) for g in genotypes} | {"WT"}
    frames = {}
    for variant, rel in sources.items():
        d = pd.read_csv(ROOT / rel)
        # Match on the DISPLAY name, not safe_label. Callers hand us whichever
        # form their summary sheet uses, and filtering "K103N+M230L" against a
        # safe_label of "K103N_M230L" silently drops every compound genotype.
        d = d.assign(genotype=d.safe_label.map(_display))
        frames[variant] = (d[d.genotype.isin(wanted)]
                           .set_index(["genotype", "replicate"]))

    # Only rows present under BOTH variants can be compared; one scored under
    # a single variant would appear as a half-empty row that reads as a missing
    # value rather than as work not yet done.
    common = frames[names[0]].index.intersection(frames[names[1]].index)
    # Sort ONCE and index everything with the sorted order. Building the frame
    # from common.sort_values() while pulling values with the unsorted common
    # silently pairs each row with another row's numbers.
    common = common.sort_values()
    out = pd.DataFrame(index=common).reset_index()
    out["dor_fold_reduction"] = frames[names[0]].loc[common, "fold_reduction"].values

    # Dropping genotypes quietly is how the compound ones went missing before.
    # Say which, and why, rather than shipping a shorter sheet than the summary.
    got = set(out.genotype)
    if wanted - got:
        missing = {g: [v for v in names
                       if g not in frames[v].index.get_level_values("genotype")]
                   for g in sorted(wanted - got)}
        raise AssertionError(
            "per-replicate sheet would omit genotypes the summary reports: "
            + "; ".join(f"{g} (absent from {', '.join(v)})" for g, v in missing.items()))

    cols = ["genotype", "replicate", "dor_fold_reduction"]
    for col, stem in _ENERGY_COLUMNS.items():
        for variant in names:
            name = f"{variant}_{stem}_kcal"
            out[name] = (frames[variant].loc[common, col].values / _KJ_PER_KCAL).round(3)
            cols.append(name)

    # Guard against the row-pairing bug above ever returning: every row must
    # carry the fold reduction its own genotype has in the manifest.
    for variant, f in frames.items():
        got = pd.Series(f.loc[common, "fold_reduction"].values)
        if not got.equals(pd.Series(out.dor_fold_reduction.values)):
            raise AssertionError(
                f"{variant} rows are not aligned with the output rows; the "
                "per-replicate sheet would pair genotypes with the wrong numbers")
    return out[cols].sort_values(["genotype", "replicate"])


def _check_sheets_agree(summary: pd.DataFrame, per_rep: pd.DataFrame,
                        variants: tuple[str, str], tol: float = 0.01) -> None:
    """The summary and per-replicate sheets must tell the same story.

    One is WT-referenced and the other absolute, built from the checkpoints by
    different routes, so nothing otherwise forces them to agree. Re-derive the
    first from the second rather than publishing two sheets that could
    silently diverge.
    """
    wt = per_rep[per_rep.genotype == "WT"]
    if wt.empty:
        raise AssertionError("per-replicate sheet has no WT rows to reference against")
    for variant in variants:
        ref = wt[f"{variant}_total_kcal"].mean()
        for g, rows in per_rep[per_rep.genotype != "WT"].groupby("genotype"):
            derived = rows[f"{variant}_total_kcal"].mean() - ref
            stated = summary[(summary.genotype == g)
                             & (summary.energy_component == "Total")][f"{variant}_ddE_kcal"]
            if stated.empty:
                continue
            if abs(derived - float(stated.iloc[0])) > tol:
                raise AssertionError(
                    f"{g} {variant} ddE Total: per-replicate sheet gives {derived:+.3f} "
                    f"but the summary states {float(stated.iloc[0]):+.3f}")


def table5() -> None:
    """Snapshot relaxation protocol: hydrogens-only against full minimisation."""
    panel = pd.read_csv(REL / "relaxation_panel_long.csv")
    agree = pd.read_csv(REL / "relaxation_protocol_agreement.csv")
    pheno = pd.read_csv(REL / "relaxation_phenotype_correlation.csv")

    summary = pd.DataFrame({
        "genotype": panel.genotype.map(_display),
        "energy_component": _stem(panel.component),
        "h_relax_ddE_kcal": panel.h_relax_mean,
        "h_relax_sem_kcal": panel.h_relax_sem,
        "full_min_ddE_kcal": panel.full_min_mean,
        "full_min_sem_kcal": panel.full_min_sem,
    })
    summary["full_min_minus_h_relax_kcal"] = (
        summary.full_min_ddE_kcal - summary.h_relax_ddE_kcal).round(3)

    agreement = pd.DataFrame({
        "energy_component": _stem(agree.component),
        "n_genotypes": agree.n,
        "spearman_h_relax_vs_full_min": agree.protocol_spearman,
        "spearman_p": agree.protocol_p,
        "sign_agreement": agree.sign_agreement,
    })

    phenotype = pd.DataFrame({
        "energy_component": _stem(pheno.component),
        "relaxation_protocol": pheno.protocol,
        "n_genotypes": pheno.n,
        "R2_vs_log10_fold": pheno.R2_vs_log10_fold,
        "ols_p": pheno.ols_p,
        "spearman_vs_log10_fold": pheno.spearman,
        "spearman_p": pheno.spearman_p,
    })

    per_rep = _per_replicate_pair(
        {"h_relax": _CHECKPOINTS["h_relax"], "full_min": _CHECKPOINTS["full_min"]},
        sorted(panel.genotype.unique()))
    _check_sheets_agree(summary, per_rep, ("h_relax", "full_min"))

    _write(OUT / "Supplementary-Table-5.xlsx",
           {"Table 5": summary,
            "Table 5b agreement": agreement,
            "Table 5b phenotype": phenotype,
            "Per replicate": per_rep})


def table6() -> None:
    """Ligand charge model: Gasteiger (used throughout) against AM1-BCC."""
    comp = pd.read_csv(REL / "charge_model_comparison.csv")
    charges = pd.read_csv(REL / "dor_charge_models.csv")

    order = ["ddE_Total", "ddE_vdW", "ddE_elec", "ddE_GB"]
    comp["component"] = pd.Categorical(comp.component, order, ordered=True)
    comp = comp.sort_values(["genotype", "component"])
    # n_gast and n_am1bcc were 3 on every row by construction; a column that
    # never varies tells the reader nothing.
    summary = pd.DataFrame({
        "genotype": comp.genotype.map(_display),
        "energy_component": _stem(comp.component),
        "gasteiger_ddE_kcal": comp.mean_gast,
        "gasteiger_sem_kcal": comp.sem_gast,
        "am1bcc_ddE_kcal": comp.mean_am1bcc,
        "am1bcc_sem_kcal": comp.sem_am1bcc,
        "am1bcc_minus_gasteiger_kcal": comp.difference,
    })

    per_rep = _per_replicate_pair(
        {"gasteiger": _CHECKPOINTS["gasteiger"], "am1bcc": _CHECKPOINTS["am1bcc"]},
        sorted(comp.genotype.unique()))
    _check_sheets_agree(summary, per_rep, ("gasteiger", "am1bcc"))

    partial_charges = pd.DataFrame({
        "atom_index": charges["index"],
        "element": charges.element,
        "gasteiger_charge_e": charges.gasteiger.round(4),
        "am1bcc_charge_e": charges.am1bcc.round(4),
        "am1bcc_minus_gasteiger_e": charges.difference.round(4),
    })

    _write(OUT / "Supplementary-Table-6.xlsx",
           {"Table 6": summary,
            "Per replicate": per_rep,
            "DOR partial charges": partial_charges})


def table7() -> None:
    """Unliganded pocket volume beside the bound-state volume."""
    apo = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_matched_window_summary.csv")
    holo = pd.read_csv(ROOT / "results/analysis/pocket_volume/tables/pocket_volume_genotype.csv")
    crystal = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_crystal_reference.csv")
    frames = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/apo_frame_metrics.csv")

    apo = apo.assign(mutation=apo.safe_label.map(_display))
    t = apo.merge(holo, on="mutation", how="left")
    summary = pd.DataFrame({
        "genotype": t.mutation,
        "n_replicates": t.n_rep,
        "apo_volume_A3": t.pocket_volume.round(1),
        "apo_volume_sd_A3": t.pocket_volume_sd.round(1),
        "holo_volume_A3": t.volume_mean.round(1),
        "holo_volume_sem_A3": t.volume_sem.round(1),
        "induced_fit_opening_A3": (t.volume_mean - t.pocket_volume).round(1),
        "apo_chi1_Y181_deg": t.chi1_Y181.round(1),
        "apo_chi1_Y188_deg": t.chi1_Y188.round(1),
        "apo_volume_vs_WT_A3": t.dVolume_vs_WT.round(1),
    }).sort_values("genotype").reset_index(drop=True)

    per_rep = (frames[frames.window == "matched"]
               .groupby(["safe_label", "replicate"], as_index=False)
               .agg(apo_volume_mean_A3=("pocket_volume", "mean"),
                    apo_volume_sd_A3=("pocket_volume", "std"),
                    n_frames=("pocket_volume", "size"),
                    apo_chi1_Y181_deg=("chi1_Y181", "mean"),
                    apo_chi1_Y188_deg=("chi1_Y188", "mean")))
    per_rep.insert(0, "genotype", per_rep.pop("safe_label").map(_display))
    for c in per_rep.columns:
        if c.endswith(("_A3", "_deg")):
            per_rep[c] = per_rep[c].round(1)

    crystal_ref = pd.DataFrame({
        "source": crystal.source,
        "state": crystal.state,
        "pocket_volume_A3": crystal.pocket_volume.round(1),
        "n_ca_atoms": crystal.n_ca,
    })

    _write(OUT / "Supplementary-Table-7.xlsx",
           {"Table 7": summary,
            "Per replicate apo": per_rep,
            "Crystal reference": crystal_ref})


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    table5()
    table6()
    table7()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
