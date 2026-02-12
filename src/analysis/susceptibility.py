from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def _normalize_mutation_label(label: str) -> str:
    tokens = [t.strip().upper() for t in re.split(r"[+,]", str(label)) if t.strip()]
    if not tokens:
        raise ValueError(f"Empty mutation label: {label!r}")
    return "+".join(tokens)


def load_dor_susceptibilities(
    excel_path: Path,
    default_chain: str = "A",
) -> pd.DataFrame:
    df = pd.read_excel(excel_path, header=1)
    required = {"Mutations", "DOR (Median-fold reduction)"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in {excel_path}: {', '.join(sorted(missing))}"
        )

    dor = df[df["DOR (Median-fold reduction)"].notna()].copy()
    dor["mutation"] = dor["Mutations"].map(_normalize_mutation_label)
    dor["dor_fold_reduction"] = dor["DOR (Median-fold reduction)"].astype(float)
    dor["chain"] = default_chain
    dor["drug"] = "DOR"
    dor = dor.drop_duplicates(subset=["mutation"]).reset_index(drop=True)
    dor["order"] = range(len(dor))
    return dor[["drug", "mutation", "chain", "dor_fold_reduction", "order"]]
