#!/usr/bin/env python3
"""Exhaustive small-combo logistic search with inner-CV threshold tuning for sensitivity."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold

from .search_feature_combo_logistic import (
    META_COLUMNS,
    _binary_labels,
    _classification_metrics,
    _combo_count,
    _feature_target_correlations,
    _logistic_pipeline,
    _parse_float_list,
    _parse_int_list,
)


def _threshold_metrics(y_true: np.ndarray, prob_high: np.ndarray, threshold: float) -> dict[str, float]:
    preds = np.where(prob_high >= float(threshold), "high", "low")
    metrics = _classification_metrics(y_true.astype(str), preds.astype(str))
    cm = confusion_matrix(y_true.astype(str), preds.astype(str), labels=["low", "high"])
    return {
        "threshold": float(threshold),
        "fn": int(cm[1, 0]),
        "fp": int(cm[0, 1]),
        "tn": int(cm[0, 0]),
        "tp": int(cm[1, 1]),
        **metrics,
    }


def _best_threshold_from_oof(y_true: np.ndarray, prob_high: np.ndarray) -> dict[str, float]:
    candidates = sorted(set(float(x) for x in prob_high.tolist()))
    if 0.5 not in candidates:
        candidates.append(0.5)
    rows = [_threshold_metrics(y_true, prob_high, threshold) for threshold in sorted(candidates)]
    threshold_df = pd.DataFrame(rows).sort_values(
        ["fn", "balanced_accuracy", "macro_f1", "accuracy", "threshold"],
        ascending=[True, False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)
    return threshold_df.iloc[0].to_dict()


def _fit_probabilities(
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_test: pd.DataFrame,
    *,
    penalty: str,
    c_value: float,
    random_state: int,
) -> np.ndarray:
    pipeline = _logistic_pipeline(random_state=int(random_state), penalty=str(penalty))
    pipeline.set_params(model__C=float(c_value))
    pipeline.fit(x_train, y_train)
    classes = list(pipeline.named_steps["model"].classes_)
    pos_idx = int(classes.index("high"))
    return pipeline.predict_proba(x_test)[:, pos_idx].astype(float)


def _choose_c_and_threshold(
    x_train: pd.DataFrame,
    y_train: pd.Series,
    *,
    penalty: str,
    c_values: list[float],
    random_state: int,
    combo_label: str,
    fold_idx: int,
) -> dict[str, float]:
    y_train_np = y_train.to_numpy(dtype=str)
    inner_splits = max(2, min(3, int(y_train.value_counts().min())))
    inner = StratifiedKFold(
        n_splits=inner_splits,
        shuffle=True,
        random_state=int(random_state) + hash((combo_label, penalty, fold_idx, "threshold")) % 100000,
    )

    candidate_rows: list[dict[str, float]] = []
    for c_value in c_values:
        oof_prob = np.full(len(y_train), np.nan, dtype=float)
        for inner_train_idx, inner_val_idx in inner.split(x_train, y_train):
            x_inner_train = x_train.iloc[inner_train_idx]
            y_inner_train = y_train.iloc[inner_train_idx]
            x_inner_val = x_train.iloc[inner_val_idx]
            oof_prob[inner_val_idx] = _fit_probabilities(
                x_inner_train,
                y_inner_train,
                x_inner_val,
                penalty=str(penalty),
                c_value=float(c_value),
                random_state=int(random_state),
            )
        best_threshold = _best_threshold_from_oof(y_train_np, oof_prob)
        candidate_rows.append(
            {
                "c_value": float(c_value),
                "threshold": float(best_threshold["threshold"]),
                "inner_fn": int(best_threshold["fn"]),
                "inner_fp": int(best_threshold["fp"]),
                "inner_tn": int(best_threshold["tn"]),
                "inner_tp": int(best_threshold["tp"]),
                "inner_accuracy": float(best_threshold["accuracy"]),
                "inner_balanced_accuracy": float(best_threshold["balanced_accuracy"]),
                "inner_macro_f1": float(best_threshold["macro_f1"]),
                "inner_macro_precision": float(best_threshold["macro_precision"]),
                "inner_macro_recall": float(best_threshold["macro_recall"]),
            }
        )

    candidate_df = pd.DataFrame(candidate_rows).sort_values(
        ["inner_fn", "inner_balanced_accuracy", "inner_macro_f1", "inner_accuracy", "threshold", "c_value"],
        ascending=[True, False, False, False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    return candidate_df.iloc[0].to_dict()


def _fullfit_coefficients(
    x: pd.DataFrame,
    y: pd.Series,
    *,
    combo: tuple[str, ...],
    penalty: str,
    random_state: int,
    c_value: float,
) -> tuple[list[dict[str, object]], float]:
    fitted = _logistic_pipeline(random_state=int(random_state), penalty=str(penalty))
    fitted.set_params(model__C=float(c_value))
    fitted.fit(x, y)
    model = fitted.named_steps["model"]
    coef_raw = np.asarray(model.coef_, dtype=float).reshape(-1)
    intercept_raw = float(np.asarray(model.intercept_, dtype=float).reshape(-1)[0])
    if list(model.classes_)[1] == "high":
        coef = coef_raw.astype(float)
        intercept = float(intercept_raw)
    else:
        coef = (-coef_raw).astype(float)
        intercept = float(-intercept_raw)

    rows: list[dict[str, object]] = []
    combo_label = "|".join(combo)
    for feature_name, coef_value in zip(combo, coef.tolist()):
        rows.append(
            {
                "penalty": str(penalty),
                "combo_size": int(len(combo)),
                "feature_combo": combo_label,
                "feature": str(feature_name),
                "coefficient": float(coef_value),
                "abs_coefficient": float(abs(coef_value)),
                "direction": "toward_high" if coef_value >= 0.0 else "toward_low",
                "fullfit_best_c": float(c_value),
                "intercept": float(intercept),
            }
        )
    return rows, intercept


def _evaluate_combo_penalty(
    feat: pd.DataFrame,
    *,
    combo: tuple[str, ...],
    combo_size: int,
    feature_count_filtered: int,
    target_col: str,
    low_max_fold: float,
    cv_folds: int,
    random_state: int,
    penalty: str,
    c_values: list[float],
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    combo_label = "|".join(combo)
    x = feat[list(combo)].copy()
    y_value = pd.to_numeric(feat[target_col], errors="coerce").astype(float)
    y = _binary_labels(y_value, low_max=float(low_max_fold))
    effective_folds = max(2, min(int(cv_folds), int(y.value_counts().min())))
    outer = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(random_state))

    preds = np.empty(len(y), dtype=object)
    prob_high = np.full(len(y), np.nan, dtype=float)
    chosen_c: list[float] = []
    chosen_thresholds: list[float] = []
    inner_fn: list[int] = []
    inner_fp: list[int] = []
    inner_bal_acc: list[float] = []
    pred_rows: list[dict[str, object]] = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
        x_train = x.iloc[train_idx]
        y_train = y.iloc[train_idx]
        x_test = x.iloc[test_idx]
        selection = _choose_c_and_threshold(
            x_train,
            y_train,
            penalty=str(penalty),
            c_values=list(c_values),
            random_state=int(random_state),
            combo_label=combo_label,
            fold_idx=int(fold_idx),
        )
        fitted_prob = _fit_probabilities(
            x_train,
            y_train,
            x_test,
            penalty=str(penalty),
            c_value=float(selection["c_value"]),
            random_state=int(random_state),
        )
        fold_pred = np.where(fitted_prob >= float(selection["threshold"]), "high", "low")
        preds[test_idx] = fold_pred.astype(str)
        prob_high[test_idx] = fitted_prob.astype(float)
        chosen_c.append(float(selection["c_value"]))
        chosen_thresholds.append(float(selection["threshold"]))
        inner_fn.append(int(selection["inner_fn"]))
        inner_fp.append(int(selection["inner_fp"]))
        inner_bal_acc.append(float(selection["inner_balanced_accuracy"]))

        for sample_idx, sample_prob, sample_pred in zip(test_idx, fitted_prob.tolist(), fold_pred.tolist()):
            pred_rows.append(
                {
                    "penalty": str(penalty),
                    "combo_size": int(combo_size),
                    "feature_combo": combo_label,
                    "fold": int(fold_idx),
                    "decision_threshold": float(selection["threshold"]),
                    "selected_c": float(selection["c_value"]),
                    "mutation": str(feat.iloc[sample_idx]["mutation"]),
                    "target_value": float(y_value.iloc[sample_idx]),
                    "observed_class": str(y.iloc[sample_idx]),
                    "predicted_class": str(sample_pred),
                    "prob_high": float(sample_prob),
                    "prob_low": float(1.0 - sample_prob),
                }
            )

    y_np = y.to_numpy(dtype=str)
    y_bin = (y == "high").astype(int).to_numpy(dtype=int)
    metrics = _classification_metrics(y_np, preds.astype(str))
    from sklearn.metrics import average_precision_score, roc_auc_score

    metrics["roc_auc"] = float(roc_auc_score(y_bin, prob_high))
    metrics["average_precision"] = float(average_precision_score(y_bin, prob_high))
    cm = confusion_matrix(y_np, preds.astype(str), labels=["low", "high"])
    summary_row = {
        "penalty": str(penalty),
        "combo_size": int(combo_size),
        "feature_combo": combo_label,
        "n_features_filtered": int(feature_count_filtered),
        "cv_folds": int(effective_folds),
        "low_max_fold": float(low_max_fold),
        "threshold_objective": "min_fn_then_max_balanced_accuracy",
        "c_by_fold": json.dumps(chosen_c),
        "threshold_by_fold": json.dumps(chosen_thresholds),
        "c_mean": float(np.mean(chosen_c)),
        "c_median": float(np.median(chosen_c)),
        "threshold_mean": float(np.mean(chosen_thresholds)),
        "threshold_median": float(np.median(chosen_thresholds)),
        "threshold_min": float(np.min(chosen_thresholds)),
        "threshold_max": float(np.max(chosen_thresholds)),
        "inner_fn_mean": float(np.mean(inner_fn)),
        "inner_fp_mean": float(np.mean(inner_fp)),
        "inner_balanced_accuracy_mean": float(np.mean(inner_bal_acc)),
        "tn": int(cm[0, 0]),
        "fp": int(cm[0, 1]),
        "fn": int(cm[1, 0]),
        "tp": int(cm[1, 1]),
        **metrics,
    }

    full_selection = _choose_c_and_threshold(
        x,
        y,
        penalty=str(penalty),
        c_values=list(c_values),
        random_state=int(random_state),
        combo_label=combo_label,
        fold_idx=0,
    )
    coef_rows, _intercept = _fullfit_coefficients(
        x,
        y,
        combo=combo,
        penalty=str(penalty),
        random_state=int(random_state),
        c_value=float(full_selection["c_value"]),
    )
    for row in coef_rows:
        row["fullfit_best_threshold"] = float(full_selection["threshold"])
        row["fullfit_inner_fn"] = int(full_selection["inner_fn"])
        row["fullfit_inner_fp"] = int(full_selection["inner_fp"])
        row["fullfit_inner_balanced_accuracy"] = float(full_selection["inner_balanced_accuracy"])

    cm_rows = [
        {
            "penalty": str(penalty),
            "combo_size": int(combo_size),
            "feature_combo": combo_label,
            "observed": "obs_low",
            "pred_low": int(cm[0, 0]),
            "pred_high": int(cm[0, 1]),
        },
        {
            "penalty": str(penalty),
            "combo_size": int(combo_size),
            "feature_combo": combo_label,
            "observed": "obs_high",
            "pred_low": int(cm[1, 0]),
            "pred_high": int(cm[1, 1]),
        },
    ]
    return summary_row, pred_rows, coef_rows, cm_rows


def _evaluate_combo(
    feat: pd.DataFrame,
    *,
    combo: tuple[str, ...],
    combo_size: int,
    feature_count_filtered: int,
    target_col: str,
    low_max_fold: float,
    cv_folds: int,
    random_state: int,
    penalty_c_map: dict[str, list[float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    summary_rows: list[dict[str, object]] = []
    pred_rows: list[dict[str, object]] = []
    coef_rows: list[dict[str, object]] = []
    cm_rows: list[dict[str, object]] = []
    for penalty, c_values in penalty_c_map.items():
        summary_row, combo_pred_rows, combo_coef_rows, combo_cm_rows = _evaluate_combo_penalty(
            feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(feature_count_filtered),
            target_col=str(target_col),
            low_max_fold=float(low_max_fold),
            cv_folds=int(cv_folds),
            random_state=int(random_state),
            penalty=str(penalty),
            c_values=list(c_values),
        )
        summary_rows.append(summary_row)
        pred_rows.extend(combo_pred_rows)
        coef_rows.extend(combo_coef_rows)
        cm_rows.extend(combo_cm_rows)
    return summary_rows, pred_rows, coef_rows, cm_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Exhaustive small-combo logistic search with sensitivity-prioritized threshold tuning.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_including_energy/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_logistic_sensitivity"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--min-r-squared", type=float, default=0.1)
    parser.add_argument("--combo-sizes", type=str, default="2,3,4,5")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--l2-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--l1-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--parallel-backend", type=str, default="loky", choices=["threading", "loky"])
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    combo_sizes = _parse_int_list(args.combo_sizes)
    l2_c_values = _parse_float_list(args.l2_c_values)
    l1_c_values = _parse_float_list(args.l1_c_values)

    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feat = pd.read_csv(args.feature_matrix_csv)
    feature_cols = [column for column in feat.columns if column not in META_COLUMNS]
    corr_df = _feature_target_correlations(
        feat,
        feature_cols=feature_cols,
        target_col=str(args.target_col),
    )
    filtered_df = corr_df[corr_df["r_squared"] >= float(args.min_r_squared)].copy().reset_index(drop=True)
    filtered_features = filtered_df["feature"].astype(str).tolist()

    penalty_c_map = {
        "l2": l2_c_values,
        "l1": l1_c_values,
    }

    tasks: list[tuple[tuple[str, ...], int]] = []
    for combo_size in sorted(set(int(size) for size in combo_sizes)):
        if combo_size <= 0 or combo_size > len(filtered_features):
            continue
        for combo in itertools.combinations(filtered_features, combo_size):
            tasks.append((combo, combo_size))

    parallel_kwargs = {"n_jobs": int(args.n_jobs), "verbose": 10}
    if str(args.parallel_backend) == "threading":
        parallel_kwargs["prefer"] = "threads"
    results = Parallel(**parallel_kwargs)(
        delayed(_evaluate_combo)(
            feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(len(filtered_features)),
            target_col=str(args.target_col),
            low_max_fold=float(args.low_max_fold),
            cv_folds=int(args.cv_folds),
            random_state=int(args.random_state),
            penalty_c_map={key: list(values) for key, values in penalty_c_map.items()},
        )
        for combo, combo_size in tasks
    )

    summary_rows = [row for combo_summary, _pred, _coef, _cm in results for row in combo_summary]
    pred_rows = [row for _summary, combo_pred, _coef, _cm in results for row in combo_pred]
    coef_rows = [row for _summary, _pred, combo_coef, _cm in results for row in combo_coef]
    cm_rows = [row for _summary, _pred, _coef, combo_cm in results for row in combo_cm]

    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["fn", "balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[True, False, False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)
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

    top_rows: list[pd.DataFrame] = []
    if not summary_df.empty:
        for penalty in sorted(summary_df["penalty"].astype(str).unique()):
            for combo_size in sorted(summary_df["combo_size"].astype(int).unique()):
                sub = summary_df[(summary_df["penalty"].astype(str) == penalty) & (summary_df["combo_size"].astype(int) == combo_size)]
                if not sub.empty:
                    top_rows.append(sub.head(10).copy())
    top_df = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()
    top_df.to_csv(out_tables / "top_combo_models.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "target_col": str(args.target_col),
                "min_r_squared": float(args.min_r_squared),
                "combo_sizes": combo_sizes,
                "cv_folds": int(args.cv_folds),
                "random_state": int(args.random_state),
                "low_max_fold": float(args.low_max_fold),
                "threshold_objective": "min_fn_then_max_balanced_accuracy",
                "l2_c_values": l2_c_values,
                "l1_c_values": l1_c_values,
                "n_mutations": int(len(feat)),
                "n_features_input": int(len(feature_cols)),
                "n_features_filtered": int(len(filtered_df)),
                "filtered_features": filtered_features,
                "n_feature_combinations": _combo_count(len(filtered_features), combo_sizes),
                "penalties": ["l2", "l1"],
                "n_jobs": int(args.n_jobs),
                "parallel_backend": str(args.parallel_backend),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
