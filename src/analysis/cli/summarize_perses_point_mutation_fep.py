from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize Perses point-mutation FEP cycle results.")
    parser.add_argument("--results-dir", type=Path, default=Path("results/analysis/perses_point_mutation_fep"))
    parser.add_argument("--mmgbsa", type=Path, default=Path("results/analysis/binding_energy/tables/ddg_full.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/analysis/perses_point_mutation_fep/perses_fep_vs_mmgbsa.csv"))
    args = parser.parse_args()

    rows = []
    for path in sorted(args.results_dir.glob("*/**/summary.json")):
        row = json.loads(path.read_text())
        row["summary_json"] = str(path)
        rows.append(row)
    if not rows:
        raise SystemExit(f"No summary.json files found under {args.results_dir}")

    out = pd.DataFrame(rows)
    if args.mmgbsa.exists():
        mm = pd.read_csv(args.mmgbsa)
        if "ddg" in mm.columns and "safe_label" in mm.columns:
            mm_mean = (
                mm.groupby("safe_label", as_index=False)
                .agg(mmgbsa_ddg_kj_mol=("ddg", "mean"), mmgbsa_ddg_sem_kj_mol=("ddg", "sem"))
            )
            out = out.merge(mm_mean, left_on="mutation", right_on="safe_label", how="left")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
