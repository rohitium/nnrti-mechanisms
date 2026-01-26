from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_metrics_xlsx(df: pd.DataFrame, xlsx_path: Path) -> None:
    cols_to_drop = [
        c for c in ["mutation_order", "chain", "subunit", "structure"] if c in df.columns
    ]
    with pd.ExcelWriter(xlsx_path) as writer:
        for structure in ["RPV", "DOR"]:
            sheet_df = df[df["structure"] == structure].copy()
            if sheet_df.empty:
                continue
            wt = sheet_df[sheet_df["state"] == "WT"].copy()
            mut = sheet_df[sheet_df["state"] == "MUT"].copy()
            wt_single = wt.sort_values(["metric"]).drop_duplicates(subset=["metric"])
            out_df = pd.concat([wt_single, mut], ignore_index=True)
            if cols_to_drop:
                out_df = out_df.drop(columns=cols_to_drop)
            out_df.to_excel(writer, sheet_name=structure, index=False)
