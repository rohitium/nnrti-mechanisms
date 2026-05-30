#!/usr/bin/env python3
"""Train combo logistic models on curated control mutations and score uncertain genotypes."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .plot_best_combo_logistic import main as plot_best_combo_logistic_main
from .run_logistic_regression_parsimonious import _wt_feature_row
from .search_feature_combo_logistic import (
    META_COLUMNS,
    _binary_labels,
    _evaluate_combo as _evaluate_combo_standard,
    _feature_target_correlations,
    _logistic_pipeline,
    _parse_float_list,
    _parse_int_list,
)
from .search_feature_combo_logistic_sensitivity import (
    _evaluate_combo as _evaluate_combo_sensitivity,
)


NEGATIVE_CONTROLS = {
    "K103N",
    "Y181C",
    "G190A",
    "V106I",
    "F227C",
}

POSITIVE_CONTROLS = {
    "V106A",
    "Y188L",
    "Y318F",
    "A98G+F227C",
    "V106A+F227L",
    "V106A+L234I",
    "V106A+P225H",
    "V106I+F227C",
    "K103N+M230L",
}

UNCERTAIN_LIMITED = {
    "L100I+K103N",
    "K103N+P225H",
    "V106M",
    "G190E",
    "G190S",
}


def _category_for_mutation(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation in NEGATIVE_CONTROLS:
        return "negative_control"
    if mutation in POSITIVE_CONTROLS:
        return "positive_control"
    if mutation in UNCERTAIN_LIMITED:
        return "uncertain_limited"
    raise ValueError(f"Mutation is not assigned to a control category: {label}")


def _ensure_control_categories(feat: pd.DataFrame) -> pd.DataFrame:
    out = feat.copy()
    out["control_category"] = out["mutation"].astype(str).map(_category_for_mutation)
    return out


def _tasks(filtered_features: list[str], combo_sizes: list[int]) -> list[tuple[tuple[str, ...], int]]:
    tasks: list[tuple[tuple[str, ...], int]] = []
    for combo_size in sorted(set(int(size) for size in combo_sizes)):
        if combo_size <= 0 or combo_size > len(filtered_features):
            continue
        for combo in itertools.combinations(filtered_features, combo_size):
            tasks.append((combo, combo_size))
    return tasks


def _parallel_results(
    evaluator,
    feat: pd.DataFrame,
    *,
    tasks: list[tuple[tuple[str, ...], int]],
    filtered_features: list[str],
    target_col: str,
    low_max_fold: float,
    cv_folds: int,
    random_state: int,
    penalty_c_map: dict[str, list[float]],
    inner_scoring: str | None,
    n_jobs: int,
    parallel_backend: str,
) -> list[tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]]:
    parallel_kwargs: dict[str, object] = {"n_jobs": int(n_jobs), "verbose": 10}
    if str(parallel_backend) == "threading":
        parallel_kwargs["prefer"] = "threads"

    def _job(combo: tuple[str, ...], combo_size: int):
        kwargs = dict(
            feat=feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(len(filtered_features)),
            target_col=str(target_col),
            low_max_fold=float(low_max_fold),
            cv_folds=int(cv_folds),
            random_state=int(random_state),
            penalty_c_map={key: list(values) for key, values in penalty_c_map.items()},
        )
        if inner_scoring is not None:
            kwargs["inner_scoring"] = str(inner_scoring)
        return evaluator(**kwargs)

    return Parallel(**parallel_kwargs)(delayed(_job)(combo, combo_size) for combo, combo_size in tasks)


def _write_combo_outputs(
    output_dir: Path,
    *,
    feat_train: pd.DataFrame,
    corr_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    results,
    sort_columns: list[str],
    ascending: list[bool],
    config_extra: dict[str, object],
) -> None:
    out_tables = output_dir / "tables"
    out_config = output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    summary_rows = [row for combo_summary, _pred, _coef, _cm in results for row in combo_summary]
    pred_rows = [row for _summary, combo_pred, _coef, _cm in results for row in combo_pred]
    coef_rows = [row for _summary, _pred, combo_coef, _cm in results for row in combo_coef]
    cm_rows = [row for _summary, _pred, _coef, combo_cm in results for row in combo_cm]

    summary_df = pd.DataFrame(summary_rows).sort_values(sort_columns, ascending=ascending, kind="stable").reset_index(drop=True)
    pred_df = pd.DataFrame(pred_rows)
    coef_df = pd.DataFrame(coef_rows).sort_values(
        ["penalty", "combo_size", "feature_combo", "abs_coefficient"],
        ascending=[True, True, True, False],
        kind="stable",
    ).reset_index(drop=True)
    cm_df = pd.DataFrame(cm_rows)

    corr_df.to_csv(out_tables / "feature_target_correlations.csv", index=False)
    filtered_df.to_csv(out_tables / "filtered_features.csv", index=False)
    summary_df.to_csv(out_tables / "combo_model_summary.csv", index=False)
    pred_df.to_csv(out_tables / "combo_cv_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "combo_fullfit_coefficients.csv", index=False)
    cm_df.to_csv(out_tables / "combo_confusion_matrices.csv", index=False)
    summary_df.head(25).to_csv(out_tables / "top_combo_models.csv", index=False)

    (out_config / "run_config.json").write_text(json.dumps(config_extra, indent=2))


def _score_heldout_and_wt(
    *,
    feat_train: pd.DataFrame,
    feat_holdout: pd.DataFrame,
    summary_path: Path,
    coef_path: Path,
    frame_feature_csv: Path,
    mmgbsa_replicate_csv: Path,
    susceptibility_xlsx: Path,
    low_max_fold: float,
    output_csv: Path,
    threshold_mode: str,
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
    if threshold_mode == "standard":
        decision_threshold = 0.5
    elif threshold_mode == "sensitivity":
        decision_threshold = float(coef_sub["fullfit_best_threshold"].iloc[0])
    else:
        raise ValueError(f"Unexpected threshold_mode: {threshold_mode}")

    x_train = feat_train[features].copy()
    y_train = _binary_labels(pd.to_numeric(feat_train["target_fold_reduction"], errors="coerce").astype(float), low_max=float(low_max_fold))
    fitted = _logistic_pipeline(random_state=0, penalty=penalty)
    fitted.set_params(model__C=float(c_value))
    fitted.fit(x_train, y_train)
    classes = list(fitted.named_steps["model"].classes_)
    pos_idx = int(classes.index("high"))

    rows: list[dict[str, object]] = []
    for _, row in feat_holdout.iterrows():
        x_row = row[features].to_frame().T
        prob_high = float(fitted.predict_proba(x_row)[:, pos_idx][0])
        rows.append(
            {
                "mutation": str(row["mutation"]),
                "control_category": str(row["control_category"]),
                "target_fold_reduction": float(row["target_fold_reduction"]),
                "prob_high": float(prob_high),
                "prob_low": float(1.0 - prob_high),
                "predicted_class": "high" if prob_high >= decision_threshold else "low",
                "decision_threshold": float(decision_threshold),
                "penalty": penalty,
                "c_value": float(c_value),
                "feature_combo": "|".join(features),
            }
        )

    wt_df = _wt_feature_row(
        frame_feature_csv=frame_feature_csv,
        mmgbsa_replicate_csv=mmgbsa_replicate_csv,
        susceptibility_xlsx=susceptibility_xlsx,
    )
    if not wt_df.empty:
        x_wt = wt_df[features].copy()
        wt_prob = float(fitted.predict_proba(x_wt)[:, pos_idx][0])
        rows.append(
            {
                "mutation": "WT",
                "control_category": "wt_reference",
                "target_fold_reduction": 1.0,
                "prob_high": float(wt_prob),
                "prob_low": float(1.0 - wt_prob),
                "predicted_class": "high" if wt_prob >= decision_threshold else "low",
                "decision_threshold": float(decision_threshold),
                "penalty": penalty,
                "c_value": float(c_value),
                "feature_combo": "|".join(features),
            }
        )

    out = pd.DataFrame(rows).sort_values(["control_category", "prob_high", "mutation"], ascending=[True, False, True], kind="stable").reset_index(drop=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_csv, index=False)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Train exhaustive logistic models on curated control mutations only.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_including_energy/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression_controls"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--min-r-squared", type=float, default=0.1)
    parser.add_argument("--combo-sizes", type=str, default="1,2,3,4,5")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--inner-scoring", type=str, default="balanced_accuracy")
    parser.add_argument("--l2-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--l1-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--parallel-backend", type=str, default="loky", choices=["threading", "loky"])
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--mmgbsa-replicate-csv",
        type=Path,
        default=Path("results/mmgbsa_replicate_metrics.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument("--top-n-diagnostics", type=int, default=4)
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    combo_sizes = _parse_int_list(args.combo_sizes)
    l2_c_values = _parse_float_list(args.l2_c_values)
    l1_c_values = _parse_float_list(args.l1_c_values)

    output_root = args.output_dir
    out_tables = output_root / "tables"
    out_config = output_root / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feat_all = _ensure_control_categories(pd.read_csv(args.feature_matrix_csv))
    feat_train = feat_all[feat_all["control_category"].isin({"negative_control", "positive_control"})].copy().reset_index(drop=True)
    feat_holdout = feat_all[feat_all["control_category"].astype(str) == "uncertain_limited"].copy().reset_index(drop=True)

    feature_cols = [column for column in feat_train.columns if column not in META_COLUMNS | {"control_category"}]
    corr_df = _feature_target_correlations(feat_train, feature_cols=feature_cols, target_col=str(args.target_col))
    filtered_df = corr_df[corr_df["r_squared"] >= float(args.min_r_squared)].copy().reset_index(drop=True)
    filtered_features = filtered_df["feature"].astype(str).tolist()
    penalty_c_map = {"l2": l2_c_values, "l1": l1_c_values}
    tasks = _tasks(filtered_features, combo_sizes)

    category_table = feat_all[["mutation", "target_fold_reduction", "control_category"]].sort_values(
        ["control_category", "target_fold_reduction", "mutation"], ascending=[True, True, True], kind="stable"
    )
    category_table.to_csv(out_tables / "mutation_control_categories.csv", index=False)

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
        },
    )

    sensitivity_dir = output_root / "sensitivity"
    sensitivity_results = _parallel_results(
        _evaluate_combo_sensitivity,
        feat_train,
        tasks=tasks,
        filtered_features=filtered_features,
        target_col=str(args.target_col),
        low_max_fold=float(args.low_max_fold),
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        penalty_c_map=penalty_c_map,
        inner_scoring=None,
        n_jobs=int(args.n_jobs),
        parallel_backend=str(args.parallel_backend),
    )
    _write_combo_outputs(
        sensitivity_dir,
        feat_train=feat_train,
        corr_df=corr_df,
        filtered_df=filtered_df,
        results=sensitivity_results,
        sort_columns=["fn", "balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[True, False, False, False, False, False],
        config_extra={
            "training_design": "negative_vs_positive_controls_only",
            "feature_matrix_csv": str(args.feature_matrix_csv),
            "output_dir": str(sensitivity_dir),
            "target_col": str(args.target_col),
            "min_r_squared": float(args.min_r_squared),
            "combo_sizes": combo_sizes,
            "cv_folds": int(args.cv_folds),
            "random_state": int(args.random_state),
            "low_max_fold": float(args.low_max_fold),
            "threshold_objective": "min_fn_then_max_balanced_accuracy",
            "l2_c_values": l2_c_values,
            "l1_c_values": l1_c_values,
            "n_jobs": int(args.n_jobs),
            "parallel_backend": str(args.parallel_backend),
            "n_training_mutations": int(len(feat_train)),
            "n_holdout_mutations": int(len(feat_holdout)),
            "filtered_features": filtered_features,
        },
    )

    _score_heldout_and_wt(
        feat_train=feat_train,
        feat_holdout=feat_holdout,
        summary_path=standard_dir / "tables" / "combo_model_summary.csv",
        coef_path=standard_dir / "tables" / "combo_fullfit_coefficients.csv",
        frame_feature_csv=args.frame_feature_csv,
        mmgbsa_replicate_csv=args.mmgbsa_replicate_csv,
        susceptibility_xlsx=args.susceptibility_xlsx,
        low_max_fold=float(args.low_max_fold),
        output_csv=standard_dir / "tables" / "heldout_uncertain_and_wt_predictions.csv",
        threshold_mode="standard",
    )
    _score_heldout_and_wt(
        feat_train=feat_train,
        feat_holdout=feat_holdout,
        summary_path=sensitivity_dir / "tables" / "combo_model_summary.csv",
        coef_path=sensitivity_dir / "tables" / "combo_fullfit_coefficients.csv",
        frame_feature_csv=args.frame_feature_csv,
        mmgbsa_replicate_csv=args.mmgbsa_replicate_csv,
        susceptibility_xlsx=args.susceptibility_xlsx,
        low_max_fold=float(args.low_max_fold),
        output_csv=sensitivity_dir / "tables" / "heldout_uncertain_and_wt_predictions.csv",
        threshold_mode="sensitivity",
    )

    plot_best_combo_logistic_main_args = [
        ("standard", standard_dir, standard_dir / "diagnostics"),
        ("sensitivity", sensitivity_dir, sensitivity_dir / "diagnostics"),
    ]
    import sys

    argv_prev = sys.argv[:]
    try:
        for _name, input_dir, output_dir in plot_best_combo_logistic_main_args:
            sys.argv = [
                "plot_best_combo_logistic",
                "--input-dir",
                str(input_dir),
                "--output-dir",
                str(output_dir),
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
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
