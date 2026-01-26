from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_drms(drms_path: Path) -> pd.DataFrame:
    df = pd.read_csv(drms_path)
    df["drug"] = df["drug"].str.upper()
    df["mutation"] = df["mutation"].str.strip()
    if "chain" not in df.columns:
        raise ValueError("DRM table must include a 'chain' column.")
    df["chain"] = df["chain"].astype(str).str.strip().str.upper()
    df["order"] = range(len(df))
    return df
