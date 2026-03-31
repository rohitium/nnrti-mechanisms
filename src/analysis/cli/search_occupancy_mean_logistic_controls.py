#!/usr/bin/env python3
"""Train occupancy_mean-only logistic combo models on curated controls and score uncertain genotypes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .plot_best_combo_logistic import main as plot_best_combo_logistic_main
from .search_feature_combo_logistic import (
    _binary_labels,
    _feature_target_correlations,
    _logistic_pipeline,
    _parse_float_list,
    _parse_int_list,
)
from .search_feature_combo_logistic_controls import (
    NEGATIVE_CONTROLS,
    POSITIVE_CONTROLS,
    UNCERTAIN_LIMITED,
    _parallel_results,
    _tasks,
    _write_combo_outputs,
)
from .search_feature_combo_logistic import _evaluate_combo as _evaluate_combo_standard


def _control_category(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation in NEGATIVE_CONTROLS:
        return "negative_control"
    if mutation in POSITIVE_CONTROLS:
        return "positive_control"
    if mutation in UNCERTAIN_LIMITED:
        return "uncertain_limited"
    if mutation == "WT":
        return "wt_reference"
    raise ValueError(f"Mutation is not assigned to a supported occupancy-model category: {label}")


def _prepare_feature_matrix(path: Path) -> pd.DataFrame:
    feat = pd.read_csv(path).copy()
    if "mutation" not in feat.columns or "dor_fold_reduction" not in feat.columns:
        raise ValueError(f"Unexpected occupancy feature matrix format: {path}")
    feat["mutation"] = feat["mutation"].astype(str)
    feat["target_fold_reduction"] = pd.to_numeric(feat["dor_fold_reduction"], errors="coerce").astype(float)
    feat["control_category"] = feat["mutation"].map(_control_category)
    feat["drug"] = "DOR"
    feat["chain"] = "A"
    feat["target_binary_class"] = _binary_labels(feat["target_fold_reduction"], low_max=10.0)
    return feat


def _score_holdout(
    *,
    feat_train: pd.DataFrame,
    feat_score: pd.DataFrame,
    summary_path: Path,
    coef_path: Path,
    low_max_fold: float,
    output_csv: Path,
) -> pd.DataFrame:
    summary = pd.read_csv(summary_path).iloc[0]
    coefs = pd.read_csv(coef_path)
    coef_sub = coefs[
        (coefs["penalty"].astype(str) == str(summary["penalty"]))
        & (coefs["combo_size"].astype(int) == int(summary["combo_size"]))
        & (coefs["feature_combo"].astype(str) == str(summary["feature_combo"]))
    ].copy()

    features = str(summary["feature_combo"]).split("|")
    penalty = str(summary["penalty"])
    c_value = float(coef_sub["fullfit_best_c"].iloc[0])
    decision_threshold = 0.5

    x_train = feat_train[features].copy()
    y_train = _binary_labels(
        pd.to_numeric(feat_train["target_fold_reduction"], errors="coerce").astype(float),
        low_max=float(low_max_fold),
    )
    fitted = _logistic_pipeline(random_state=0, penalty=penalty)
    fitted.set_params(model__C=float(c_value))
    fitted.fit(x_train, y_train)
    classes = list(fitted.named_steps["model"].classes_)
    pos_idx = int(classes.index("high"))

    rows: list[dict[str, object]] = []
    for _, row in feat_score.iterrows():
        x_row = row[features].to_frame().T
        prob_high = float(fitted.predict_proba(x_row)[:, pos_idx][0])
        rows.append(
            {
                "mutation": str(row["mutation"]),
                "control_category": str(row["control_category"]),
                "target_fold_reduction": float(row["target_fold_reduction"]),
                "observed_class_by_fold_cutoff": "low" if float(row["target_fold_reduction"]) < float(low_max_fold) else "high",
                "prob_high": float(prob_high),
                "prob_low": float(1.0 - prob_high),
                "predicted_class": "high" if prob_high >= decision_threshold else "low",
                "decision_threshold": float(decision_threshold),
                "penalty": penalty,
                "c_value": float(c_value),
                "feature_combo": "|".join(features),
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["control_category", "prob_high", "mutation"],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_csv, index=False)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Occupancy_mean-only logistic search on negative vs positive controls.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/contact_occupancy_feature_screen/tables/occupancy_mean_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/occupancy_mean_logistic_controls"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--min-r-squared", type=float, default=0.0)
    parser.add_argument("--combo-sizes", type=str, default="1,2,3,4")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--inner-scoring", type=str, default="balanced_accuracy")
    parser.add_argument("--l2-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--l1-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--parallel-backend", type=str, default="loky", choices=["threading", "loky"])
    parser.add_argument("--top-n-diagnostics", type=int, default=4)
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    combo_sizes = _parse_int_list(args.combo_sizes)
    l2_c_values = _parse_float_list(args.l2_c_values)
    l1_c_values = _parse_float_list(args.l1_c_values)
    penalty_c_map = {"l2": l2_c_values, "l1": l1_c_values}

    output_root = args.output_dir
    out_tables = output_root / "tables"
    out_config = output_root / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feat_all = _prepare_feature_matrix(args.feature_matrix_csv)
    feat_train = feat_all[feat_all["control_category"].isin({"negative_control", "positive_control"})].copy().reset_index(drop=True)
    feat_holdout = feat_all[feat_all["control_category"].isin({"uncertain_limited", "wt_reference"})].copy().reset_index(drop=True)

    feature_cols = [c for c in feat_train.columns if c.startswith("occupancy_mean_")]
    corr_df = _feature_target_correlations(feat_train, feature_cols=feature_cols, target_col=str(args.target_col))
    filtered_df = corr_df[corr_df["r_squared"] >= float(args.min_r_squared)].copy().reset_index(drop=True)
    filtered_features = filtered_df["feature"].astype(str).tolist()
    tasks = _tasks(filtered_features, combo_sizes)

    category_table = feat_all[["mutation", "target_fold_reduction", "control_category"]].sort_values(
        ["control_category", "target_fold_reduction", "mutation"],
        ascending=[True, True, True],
        kind="stable",
    )
    category_table.to_csv(out_tables / "mutation_control_categories.csv", index=False)
    feat_train.to_csv(out_tables / "control_training_feature_matrix.csv", index=False)
    feat_holdout.to_csv(out_tables / "heldout_feature_matrix.csv", index=False)

    standard_dir = output_root / "standard"
    standard_results = _parallel_results(
        _evaluate_combo_standard,
        feat_train,
        tasks=tasks,
        filtered_features=filtered_features,
        target_col=str(args.target_col),
        low_max_fold=float(args.low_max_fold),
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        penalty_c_map=penalty_c_map,
        inner_scoring=str(args.inner_scoring),
        n_jobs=int(args.n_jobs),
        parallel_backend=str(args.parallel_backend),
    )
    _write_combo_outputs(
        standard_dir,
        feat_train=feat_train,
        corr_df=corr_df,
        filtered_df=filtered_df,
        results=standard_results,
        sort_columns=["balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[False, False, False, False, False],
        config_extra={
            "training_design": "negative_vs_positive_controls_only",
            "feature_matrix_csv": str(args.feature_matrix_csv),
            "output_dir": str(standard_dir),
            "target_col": str(args.target_col),
            "min_r_squared": float(args.min_r_squared),
            "combo_sizes": combo_sizes,
            "cv_folds": int(args.cv_folds),
            "random_state": int(args.random_state),
            "low_max_fold": float(args.low_max_fold),
            "inner_scoring": str(args.inner_scoring),
            "l2_c_values": l2_c_values,
            "l1_c_values": l1_c_values,
            "n_jobs": int(args.n_jobs),
            "parallel_backend": str(args.parallel_backend),
            "n_training_mutations": int(len(feat_train)),
            "n_holdout_mutations": int(len(feat_holdout)),
            "filtered_features": filtered_features,
            "feature_family": "occupancy_mean_only",
        },
    )

    _score_holdout(
        feat_train=feat_train,
        feat_score=feat_holdout,
        summary_path=standard_dir / "tables" / "combo_model_summary.csv",
        coef_path=standard_dir / "tables" / "combo_fullfit_coefficients.csv",
        low_max_fold=float(args.low_max_fold),
        output_csv=standard_dir / "tables" / "heldout_uncertain_and_wt_predictions.csv",
    )

    import sys

    argv_prev = sys.argv[:]
    try:
        sys.argv = [
            "plot_best_combo_logistic",
            "--input-dir",
            str(standard_dir),
            "--output-dir",
            str(standard_dir / "diagnostics"),
            "--top-n",
            str(int(args.top_n_diagnostics)),
        ]
        plot_best_combo_logistic_main()
    finally:
        sys.argv = argv_prev

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "training_design": "negative_vs_positive_controls_only",
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "target_col": str(args.target_col),
                "min_r_squared": float(args.min_r_squared),
                "combo_sizes": combo_sizes,
                "cv_folds": int(args.cv_folds),
                "random_state": int(args.random_state),
                "low_max_fold": float(args.low_max_fold),
                "inner_scoring": str(args.inner_scoring),
                "l2_c_values": l2_c_values,
                "l1_c_values": l1_c_values,
                "n_jobs": int(args.n_jobs),
                "parallel_backend": str(args.parallel_backend),
                "n_training_mutations": int(len(feat_train)),
                "n_holdout_mutations": int(len(feat_holdout)),
                "top_n_diagnostics": int(args.top_n_diagnostics),
                "negative_controls": sorted(NEGATIVE_CONTROLS),
                "positive_controls": sorted(POSITIVE_CONTROLS),
                "uncertain_limited": sorted(UNCERTAIN_LIMITED),
                "feature_family": "occupancy_mean_only",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
