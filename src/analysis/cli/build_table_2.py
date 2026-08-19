#!/usr/bin/env python3
"""Build the manuscript Table 2 reference CSV: ddE (MM/GBSA) beside ddG (pmx FEP).

ddE is the MM/GBSA endpoint interaction energy; ddG is the pmx non-equilibrium
FEP binding free energy. They are different observables and the table keeps them
in separate, explicitly labelled columns so neither is mistaken for the other.

Row order and the mutation-category labels are taken from an existing manuscript
.docx so the output can be pasted straight in; the .docx is read only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

#: Labels used in the draft that differ from the analysis keys.
DOCX_ALIASES = {"K103N+L100I": "L100I+K103N"}

#: (ddE column, header, decimal places)
DDE_COLUMNS = [
    ("ddg", "∆∆E_Total (kcal/mol)", 2),
    ("ddg_vdw", "∆∆E_vdW (kcal/mol)", 2),
    ("ddg_electrostatic", "∆∆E_elec (kcal/mol)", 2),
    ("ddg_gb", "∆∆E_GB (kcal/mol)", 2),
    ("ddg_sa", "∆∆E_SA (kcal/mol)", 3),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ddg-csv", type=Path,
                        default=root / "results/analysis/binding_energy/tables/ddg_full.csv")
    parser.add_argument("--fep-csv", type=Path,
                        default=root / "results/analysis/fep_pmx/panel_discussion_tiers.csv")
    parser.add_argument("--docx", type=Path,
                        default=root / "manuscript/DorDRM-MD-06-21-26.docx",
                        help="Read-only source of row order and mutation-category labels.")
    parser.add_argument("--table-index", type=int, default=1,
                        help="Which table in the .docx carries the energetics panel.")
    parser.add_argument("--output-csv", type=Path,
                        default=root / "manuscript/Table-2-energetics.csv")
    return parser.parse_args()


def row_order(docx: Path, table_index: int) -> list[tuple[str, str]]:
    from docx import Document

    rows = [[c.text.strip() for c in r.cells] for r in Document(str(docx)).tables[table_index].rows][1:]
    # A blank category cell in the draft is inherited from the row above.
    out, last = [], ""
    for r in rows:
        cat = r[0] or last
        last = cat
        out.append((cat, r[1]))
    return out


def main() -> int:
    args = parse_args()
    ddg = pd.read_csv(args.ddg_csv)
    grouped = ddg[ddg["mutation"] != "WT"].groupby("mutation")
    stats = {
        col: (grouped[col].mean(), grouped[col].apply(lambda s: s.std(ddof=1) / np.sqrt(s.size)))
        for col, _h, _d in DDE_COLUMNS
    }
    fep = pd.read_csv(args.fep_csv).set_index("genotype") if args.fep_csv.exists() else pd.DataFrame()

    records = []
    for category, label in row_order(args.docx, args.table_index):
        key = DOCX_ALIASES.get(label, label)
        rec = {"Mutation category": category, "Mutation": label}
        for col, header, dp in DDE_COLUMNS:
            mean, sem = stats[col]
            rec[header] = f"{mean[key]:.{dp}f} ± {sem[key]:.{dp}f}" if key in mean.index else ""
        if key in fep.index and pd.notna(fep.loc[key, "ddg_bind_kcal"]):
            rec["∆∆G_bind (kcal/mol)"] = f"{fep.loc[key,'ddg_bind_kcal']:.2f} ± {fep.loc[key,'sem_kcal']:.2f}"
            rec["∆∆G tier"] = fep.loc[key, "tier"]
        else:
            rec["∆∆G_bind (kcal/mol)"] = "not determined"
            rec["∆∆G tier"] = ""
        records.append(rec)

    table = pd.DataFrame(records)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_csv, index=False)
    print(table.to_string(index=False))
    print(f"\nWrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
