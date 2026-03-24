#!/usr/bin/env python3
"""Create ablated mutation-level feature matrices for logistic classifier comparisons."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


META_COLUMNS = ["drug", "mutation", "chain", "target_fold_reduction"]


def _load_contribution_keep_set(contribution_csv: Path, threshold: float) -> set[str]:
    contrib = pd.read_csv(contribution_csv)
    if "feature" not in contrib.columns or "contribution" not in contrib.columns:
        raise ValueError(f"Expected feature/contribution columns in {contribution_csv}")
    agg = (
        contrib.assign(abs_contribution=contrib["contribution"].abs())
        .groupby("feature", as_index=False)["abs_contribution"]
        .max()
        .rename(columns={"abs_contribution": "max_abs_contribution"})
    )
    keep = agg.loc[agg["max_abs_contribution"] >= float(threshold), "feature"].astype(str).tolist()
    return set(keep)


def _subset_feature_columns(
    feature_columns: list[str],
    *,
    mode: str,
    contribution_keep: set[str] | None,
) -> list[str]:
    if mode == "drop_residue":
        return [c for c in feature_columns if not str(c).startswith("residue_min_distance_")]
    if mode == "only_residue":
        return [c for c in feature_columns if str(c).startswith("residue_min_distance_")]
    if mode == "contribution_threshold":
        if contribution_keep is None:
            raise ValueError("contribution_keep is required for contribution_threshold mode")
        return [c for c in feature_columns if str(c) in contribution_keep]
    raise ValueError(f"Unsupported mode: {mode}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create ablated logistic feature matrices.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--mode", type=str, required=True, choices=["drop_residue", "only_residue", "contribution_threshold"])
    parser.add_argument(
        "--contribution-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression/tables/full_model_feature_contributions.csv"),
    )
    parser.add_argument("--contribution-threshold", type=float, default=0.5)
    args = parser.parse_args()

    feat = pd.read_csv(args.input_csv)
    feature_columns = [c for c in feat.columns if c not in META_COLUMNS]
    contribution_keep = None
    if str(args.mode) == "contribution_threshold":
        if not args.contribution_csv.exists():
            raise FileNotFoundError(args.contribution_csv)
        contribution_keep = _load_contribution_keep_set(args.contribution_csv, float(args.contribution_threshold))
    kept_features = _subset_feature_columns(
        feature_columns,
        mode=str(args.mode),
        contribution_keep=contribution_keep,
    )
    out = feat[META_COLUMNS + kept_features].copy()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)

    config = {
        "input_csv": str(args.input_csv),
        "output_csv": str(args.output_csv),
        "mode": str(args.mode),
        "n_features_input": int(len(feature_columns)),
        "n_features_output": int(len(kept_features)),
        "feature_columns_kept": kept_features,
        "contribution_csv": str(args.contribution_csv) if str(args.mode) == "contribution_threshold" else None,
        "contribution_threshold": float(args.contribution_threshold) if str(args.mode) == "contribution_threshold" else None,
    }
    args.output_csv.with_suffix(".json").write_text(json.dumps(config, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
