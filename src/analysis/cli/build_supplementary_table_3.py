#!/usr/bin/env python3
"""Build Supplementary Table 3 from WT-referenced MM/GBSA outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter


SUMMARY_COMPONENTS = [
    ("ddg", "Total shift"),
    ("ddg_vdw", "van der Waals shift"),
    ("ddg_electrostatic", "Electrostatic shift"),
    ("ddg_gb", "GB polar solvation shift"),
    ("ddg_sa", "Nonpolar SA shift"),
]

DETAIL_COLUMNS = [
    ("mutation", "Mutation"),
    ("replicate", "Replicate"),
    ("dor_fold_reduction", "DOR fold reduction"),
    ("mmgbsa_snapshots", "MM/GBSA snapshots"),
    ("binding_dg", "Mutant total energy (kJ/mol)"),
    ("wt_binding_dg", "Matched WT total energy (kJ/mol)"),
    ("ddg", "Total shift (kJ/mol)"),
    ("binding_dg_vdw", "Mutant van der Waals (kJ/mol)"),
    ("wt_binding_dg_vdw", "Matched WT van der Waals (kJ/mol)"),
    ("ddg_vdw", "van der Waals shift (kJ/mol)"),
    ("binding_dg_electrostatic", "Mutant electrostatic (kJ/mol)"),
    ("wt_binding_dg_electrostatic", "Matched WT electrostatic (kJ/mol)"),
    ("ddg_electrostatic", "Electrostatic shift (kJ/mol)"),
    ("binding_dg_gb", "Mutant GB polar solvation (kJ/mol)"),
    ("wt_binding_dg_gb", "Matched WT GB polar solvation (kJ/mol)"),
    ("ddg_gb", "GB polar solvation shift (kJ/mol)"),
    ("binding_dg_sa", "Mutant nonpolar SA (kJ/mol)"),
    ("wt_binding_dg_sa", "Matched WT nonpolar SA (kJ/mol)"),
    ("ddg_sa", "Nonpolar SA shift (kJ/mol)"),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(
        description="Generate manuscript/Supplementary-Table-3.xlsx from WT-referenced energetic shifts."
    )
    parser.add_argument(
        "--ddg-csv",
        type=Path,
        default=root / "results/analysis/binding_energy/tables/ddg_full.csv",
        help="Replicate-level WT-referenced MM/GBSA table.",
    )
    parser.add_argument(
        "--panel-csv",
        type=Path,
        default=root / "results/analysis/dor_susceptibility_bar_chart/tables/dor_susceptibility_values.csv",
        help="Current manuscript mutation panel and DOR fold reductions.",
    )
    parser.add_argument(
        "--output-xlsx",
        type=Path,
        default=root / "manuscript/Supplementary-Table-3.xlsx",
        help="Output workbook path.",
    )
    parser.add_argument(
        "--expected-replicates",
        type=int,
        default=3,
        help="Expected number of production replicates per mutation.",
    )
    return parser.parse_args()


def sem(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if numeric.size <= 1:
        return 0.0
    return float(np.std(numeric, ddof=1) / np.sqrt(numeric.size))


def load_panel(panel_csv: Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_csv)
    required = {"mutation", "dor_fold_reduction"}
    missing = required.difference(panel.columns)
    if missing:
        raise ValueError(f"{panel_csv} is missing required columns: {sorted(missing)}")
    panel = panel[["mutation", "dor_fold_reduction"]].copy()
    panel["mutation"] = panel["mutation"].astype(str)
    return panel


def load_details(ddg_csv: Path, panel: pd.DataFrame, expected_replicates: int) -> pd.DataFrame:
    ddg = pd.read_csv(ddg_csv)
    required = {"mutation", "replicate", "mmgbsa_snapshots", *[column for column, _label in SUMMARY_COMPONENTS]}
    required.update(column for column, _label in DETAIL_COLUMNS if column in ddg.columns)
    missing = required.difference(ddg.columns)
    if missing:
        raise ValueError(f"{ddg_csv} is missing required columns: {sorted(missing)}")

    panel_order = panel["mutation"].tolist()
    filtered = ddg[ddg["mutation"].isin(panel_order)].copy()
    missing_mutations = sorted(set(panel_order).difference(filtered["mutation"].astype(str)))
    if missing_mutations:
        raise ValueError(f"No energetic rows found for panel mutations: {missing_mutations}")

    replicate_counts = filtered.groupby("mutation")["replicate"].nunique()
    bad_counts = replicate_counts[replicate_counts != expected_replicates]
    if not bad_counts.empty:
        bad = ", ".join(f"{mutation}={count}" for mutation, count in bad_counts.items())
        raise ValueError(f"Unexpected replicate counts; expected {expected_replicates}: {bad}")

    fold_lookup = panel.set_index("mutation")["dor_fold_reduction"]
    filtered["dor_fold_reduction"] = filtered["mutation"].map(fold_lookup)
    filtered["mutation"] = pd.Categorical(filtered["mutation"], categories=panel_order, ordered=True)
    filtered = filtered.sort_values(["mutation", "replicate"]).reset_index(drop=True)
    filtered["mutation"] = filtered["mutation"].astype(str)
    return filtered


def build_summary(details: pd.DataFrame) -> pd.DataFrame:
    grouped = details.groupby("mutation", sort=False)
    summary = grouped.agg(
        **{
            "DOR fold reduction": ("dor_fold_reduction", "first"),
            "n replicates": ("replicate", "nunique"),
        }
    ).reset_index()
    summary = summary.rename(columns={"mutation": "Mutation"})
    for column, label in SUMMARY_COMPONENTS:
        summary[f"{label} mean (kJ/mol)"] = grouped[column].mean().to_numpy(dtype=float)
        summary[f"{label} SEM (kJ/mol)"] = grouped[column].apply(sem).to_numpy(dtype=float)
    return summary


def build_detail_table(details: pd.DataFrame) -> pd.DataFrame:
    return details[[column for column, _label in DETAIL_COLUMNS]].rename(columns=dict(DETAIL_COLUMNS))


def style_workbook(output_xlsx: Path, summary_rows: int, detail_rows: int) -> None:
    from openpyxl import load_workbook

    wb = load_workbook(output_xlsx)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(bold=True, color="FFFFFF")
    note_font = Font(color="374151")
    title_font = Font(bold=True, size=13, color="1F2937")
    thin_border = Border(bottom=Side(style="thin", color="D1D5DB"))

    ws = wb["Summary"]
    ws["A1"] = "Notes"
    ws["A2"] = (
        "Values are WT-referenced MM/GBSA energetic shifts for doravirine-bound RT variants, reported in "
        "kJ/mol. For each mutant replicate, the energetic component from that replicate was compared with "
        "the corresponding WT replicate, so shifts are calculated as mutant minus matched WT."
    )
    ws["A3"] = (
        "The total energy is the MM/GBSA endpoint score from uniformly sampled bound-state trajectory "
        "snapshots. Component shifts are reported for van der Waals, electrostatic, generalized Born polar "
        "solvation, and nonpolar surface area terms."
    )
    ws["A4"] = (
        "Summary values are the mean across three independent production replicates. SEM is the sample "
        "standard deviation divided by the square root of the number of replicates. Positive shifts indicate "
        "a less favorable energetic score relative to WT, whereas negative shifts indicate a more favorable score."
    )
    for row in (2, 3, 4):
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=13)
        ws.cell(row=row, column=1).font = note_font
        ws.cell(row=row, column=1).alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[row].height = 42
    ws["A1"].font = title_font

    header_row = 7
    for cell in ws[header_row]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin_border
    for row in ws.iter_rows(min_row=header_row + 1, max_row=header_row + summary_rows):
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            cell.border = thin_border
            if cell.column > 1:
                cell.number_format = "0.000"
        row[2].number_format = "0"
    ws.freeze_panes = "A8"
    summary_widths = {
        "A": 18,
        "B": 13,
        "C": 12,
        "D": 16,
        "E": 16,
        "F": 18,
        "G": 18,
        "H": 18,
        "I": 18,
        "J": 20,
        "K": 20,
        "L": 18,
        "M": 18,
    }
    for column, width in summary_widths.items():
        ws.column_dimensions[column].width = width

    ws = wb["Details"]
    for cell in ws[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin_border
    for row in ws.iter_rows(min_row=2, max_row=detail_rows + 1):
        for cell in row:
            cell.border = thin_border
            if cell.column > 1:
                cell.number_format = "0.0000000000"
        row[1].number_format = "0"
        row[2].number_format = "0.000"
        row[3].number_format = "0"
    ws.freeze_panes = "A2"
    for col_idx in range(1, ws.max_column + 1):
        ws.column_dimensions[get_column_letter(col_idx)].width = 18 if col_idx == 1 else 15

    wb.save(output_xlsx)


def write_workbook(summary: pd.DataFrame, details: pd.DataFrame, output_xlsx: Path) -> None:
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        pd.DataFrame({"Notes": ["", "", ""]}).to_excel(writer, sheet_name="Summary", index=False, startrow=0)
        summary.to_excel(writer, sheet_name="Summary", index=False, startrow=6)
        details.to_excel(writer, sheet_name="Details", index=False)
    style_workbook(output_xlsx, summary_rows=len(summary), detail_rows=len(details))


def main() -> int:
    args = parse_args()
    panel = load_panel(args.panel_csv)
    details_source = load_details(args.ddg_csv, panel, args.expected_replicates)
    summary = build_summary(details_source)
    details = build_detail_table(details_source)
    write_workbook(summary, details, args.output_xlsx)
    print(f"Wrote {args.output_xlsx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
