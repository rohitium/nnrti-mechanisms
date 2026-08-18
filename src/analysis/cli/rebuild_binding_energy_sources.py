#!/usr/bin/env python3
"""Rebuild canonical binding-energy source tables from current MM/GBSA results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ..result_collector import compute_binding_ddg, merge_with_structural_metrics
from ..units import KCAL_UNITS, KJ_UNITS, convert_energy_columns, frame_energy_units


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Promote current MM/GBSA replicate metrics into canonical binding-energy CSVs."
    )
    parser.add_argument(
        "--mmgbsa-csv",
        type=Path,
        default=Path("results/.checkpoints/.checkpoint_mmgbsa_replicate_metrics.csv"),
    )
    parser.add_argument("--structural-csv", type=Path, default=Path("results/structural_metrics.csv"))
    parser.add_argument(
        "--susceptibility-csv",
        type=Path,
        default=Path("results/analysis/dor_susceptibility_bar_chart/tables/dor_susceptibility_values.csv"),
    )
    parser.add_argument("--output-mmgbsa", type=Path, default=Path("results/analysis/binding_energy/tables/mmgbsa_replicate_metrics.csv"))
    parser.add_argument("--output-ddg", type=Path, default=Path("results/analysis/binding_energy/tables/ddg_full.csv"))
    parser.add_argument(
        "--energy-units",
        choices=(KCAL_UNITS, KJ_UNITS),
        default=KCAL_UNITS,
        help=(
            "Units for all energy columns in the canonical outputs. Default kcal/mol, "
            "matching the pmx FEP outputs. Conversion is idempotent."
        ),
    )
    parser.add_argument(
        "--wt-reference",
        choices=("unmatched", "matched"),
        default="unmatched",
        help=(
            "How mutant replicates are referenced to WT. 'unmatched' (default) subtracts the "
            "WT replicate mean; 'matched' subtracts the same-index WT replicate."
        ),
    )
    parser.add_argument(
        "--config-json",
        type=Path,
        default=Path("results/analysis/binding_energy/config/source_rebuild_config.json"),
    )
    return parser.parse_args()


def _add_fold_change_alias(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "fold_reduction" in out.columns:
        if "fold_change" in out.columns:
            out["fold_change"] = out["fold_reduction"]
        else:
            insert_at = list(out.columns).index("fold_reduction") + 1
            out.insert(insert_at, "fold_change", out["fold_reduction"])
    return out


def _apply_current_fold_values(df: pd.DataFrame, susceptibility_csv: Path) -> pd.DataFrame:
    out = df.copy()
    if not susceptibility_csv.exists() or "fold_reduction" not in out.columns:
        return _add_fold_change_alias(out)

    susceptibility = pd.read_csv(susceptibility_csv)
    fold_lookup = dict(
        zip(
            susceptibility["mutation"].astype(str),
            pd.to_numeric(susceptibility["dor_fold_reduction"], errors="coerce"),
        )
    )
    labels = out["mutation"].astype(str)
    replacement = labels.map(fold_lookup)
    out.loc[replacement.notna(), "fold_reduction"] = replacement[replacement.notna()].astype(float)
    return _add_fold_change_alias(out)


def main() -> int:
    args = _parse_args()
    if not args.mmgbsa_csv.exists():
        raise FileNotFoundError(args.mmgbsa_csv)
    if not args.structural_csv.exists():
        raise FileNotFoundError(args.structural_csv)

    mmgbsa = pd.read_csv(args.mmgbsa_csv)
    structural = pd.read_csv(args.structural_csv)

    key_cols = ["structure", "mutation", "safe_label", "replicate"]
    mmgbsa = mmgbsa.drop_duplicates(subset=key_cols, keep="last").copy()
    mmgbsa["replicate"] = pd.to_numeric(mmgbsa["replicate"], errors="raise").astype(int)
    mmgbsa = _apply_current_fold_values(mmgbsa, args.susceptibility_csv)
    source_units = frame_energy_units(mmgbsa)
    mmgbsa = convert_energy_columns(mmgbsa, args.energy_units)
    mmgbsa = mmgbsa.sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)

    structural = structural.drop_duplicates(subset=key_cols, keep="last").copy()
    structural["replicate"] = pd.to_numeric(structural["replicate"], errors="raise").astype(int)
    structural = _apply_current_fold_values(structural, args.susceptibility_csv)

    ddg = compute_binding_ddg(mmgbsa, wt_reference=args.wt_reference)
    ddg = merge_with_structural_metrics(ddg, structural)
    ddg = _add_fold_change_alias(ddg)

    args.output_mmgbsa.parent.mkdir(parents=True, exist_ok=True)
    args.output_ddg.parent.mkdir(parents=True, exist_ok=True)
    mmgbsa.to_csv(args.output_mmgbsa, index=False)
    ddg.to_csv(args.output_ddg, index=False)

    if args.config_json:
        args.config_json.parent.mkdir(parents=True, exist_ok=True)
        args.config_json.write_text(
            json.dumps(
                {
                    "mmgbsa_csv": str(args.mmgbsa_csv),
                    "structural_csv": str(args.structural_csv),
                    "susceptibility_csv": str(args.susceptibility_csv),
                    "wt_reference": str(args.wt_reference),
                    "energy_units": str(args.energy_units),
                    "source_energy_units": str(source_units),
                    "output_mmgbsa": str(args.output_mmgbsa),
                    "output_ddg": str(args.output_ddg),
                    "n_mmgbsa_rows": int(len(mmgbsa)),
                    "n_ddg_rows": int(len(ddg)),
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
