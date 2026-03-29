#!/usr/bin/env python3
"""Screen correlated features and exhaustively benchmark small logistic-regression combos."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


META_COLUMNS = {"drug", "mutation", "chain", "target_fold_reduction", "target_binary_class"}


def _parse_int_list(text: str) -> list[int]:
    values = [int(token.strip()) for token in str(text).split(",") if token.strip()]
    if not values:
        raise ValueError(f"Expected at least one integer in {text!r}")
    return values


def _parse_float_list(text: str) -> list[float]:
    values = [float(token.strip()) for token in str(text).split(",") if token.strip()]
    if not values:
        raise ValueError(f"Expected at least one float in {text!r}")
    return values


def _classification_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
    }


def _binary_labels(y: pd.Series, *, low_max: float) -> pd.Series:
    return y.map(lambda value: "low" if float(value) < float(low_max) else "high")


def _feature_target_correlations(
    feat: pd.DataFrame,
    *,
    feature_cols: list[str],
    target_col: str,
) -> pd.DataFrame:
    y = pd.to_numeric(feat[target_col], errors="coerce").astype(float)
    rows: list[dict[str, object]] = []
    for feature in feature_cols:
        x = pd.to_numeric(feat[feature], errors="coerce").astype(float)
        sub = pd.DataFrame({"feature_value": x, "target_value": y}).dropna().copy()
        if len(sub) >= 3:
            pearson_r, pearson_p = stats.pearsonr(sub["target_value"], sub["feature_value"])
            spearman_rho, spearman_p = stats.spearmanr(sub["target_value"], sub["feature_value"])
            slope, intercept, r_value, p_value, _stderr = stats.linregress(sub["target_value"], sub["feature_value"])
        else:
            pearson_r = pearson_p = spearman_rho = spearman_p = slope = intercept = r_value = p_value = np.nan
        rows.append(
            {
                "feature": str(feature),
                "n_mutations": int(len(sub)),
                "pearson_r": float(pearson_r) if np.isfinite(pearson_r) else np.nan,
                "pearson_pvalue": float(pearson_p) if np.isfinite(pearson_p) else np.nan,
                "spearman_rho": float(spearman_rho) if np.isfinite(spearman_rho) else np.nan,
                "spearman_pvalue": float(spearman_p) if np.isfinite(spearman_p) else np.nan,
                "r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan,
                "slope": float(slope) if np.isfinite(slope) else np.nan,
                "intercept": float(intercept) if np.isfinite(intercept) else np.nan,
                "linregress_pvalue": float(p_value) if np.isfinite(p_value) else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["r_squared", "pearson_r"], ascending=[False, False], kind="stable").reset_index(drop=True)


def _logistic_pipeline(*, random_state: int, penalty: str) -> Pipeline:
    penalty_text = str(penalty).strip().lower()
    solver = "liblinear" if penalty_text == "l1" else "lbfgs"
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=int(random_state),
                    penalty=penalty_text,
                    solver=solver,
                ),
            ),
        ]
    )


def _grid_for_penalty(penalty: str, c_values: list[float]) -> dict[str, list[float]]:
    del penalty
    return {"model__C": [float(value) for value in c_values]}


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
    inner_scoring: str,
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
    best_scores: list[float] = []
    pred_rows: list[dict[str, object]] = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
        x_train = x.iloc[train_idx]
        y_train = y.iloc[train_idx]
        x_test = x.iloc[test_idx]
        inner_splits = max(2, min(3, int(y_train.value_counts().min())))
        inner = StratifiedKFold(
            n_splits=inner_splits,
            shuffle=True,
            random_state=int(random_state) + hash((combo_label, penalty, fold_idx)) % 100000,
        )
        search = GridSearchCV(
            _logistic_pipeline(random_state=int(random_state), penalty=str(penalty)),
            _grid_for_penalty(str(penalty), c_values),
            cv=inner.split(x_train, y_train),
            scoring=str(inner_scoring),
            n_jobs=1,
            refit=True,
        )
        search.fit(x_train, y_train)
        fitted = search.best_estimator_
        fold_pred = fitted.predict(x_test).astype(str)
        classes = list(fitted.named_steps["model"].classes_)
        pos_idx = int(classes.index("high"))
        fold_prob = fitted.predict_proba(x_test)[:, pos_idx].astype(float)
        preds[test_idx] = fold_pred
        prob_high[test_idx] = fold_prob
        chosen_c.append(float(search.best_params_["model__C"]))
        best_scores.append(float(search.best_score_))
        for sample_idx in test_idx:
            pred_rows.append(
                {
                    "penalty": str(penalty),
                    "combo_size": int(combo_size),
                    "feature_combo": combo_label,
                    "fold": int(fold_idx),
                    "mutation": str(feat.iloc[sample_idx]["mutation"]),
                    "target_value": float(y_value.iloc[sample_idx]),
                    "observed_class": str(y.iloc[sample_idx]),
                    "predicted_class": str(preds[sample_idx]),
                    "prob_high": float(prob_high[sample_idx]),
                    "prob_low": float(1.0 - prob_high[sample_idx]),
                }
            )

    y_np = y.to_numpy(dtype=str)
    y_bin = (y == "high").astype(int).to_numpy(dtype=int)
    metrics = _classification_metrics(y_np, preds.astype(str))
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
        "inner_scoring": str(inner_scoring),
        "c_by_fold": json.dumps(chosen_c),
        "c_mean": float(np.mean(chosen_c)),
        "c_median": float(np.median(chosen_c)),
        "inner_best_score_mean": float(np.mean(best_scores)),
        "inner_best_score_std": float(np.std(best_scores, ddof=1)) if len(best_scores) > 1 else 0.0,
        "tn": int(cm[0, 0]),
        "fp": int(cm[0, 1]),
        "fn": int(cm[1, 0]),
        "tp": int(cm[1, 1]),
        **metrics,
    }

    fullfit = GridSearchCV(
        _logistic_pipeline(random_state=int(random_state), penalty=str(penalty)),
        _grid_for_penalty(str(penalty), c_values),
        cv=outer.split(x, y),
        scoring=str(inner_scoring),
        n_jobs=1,
        refit=True,
    )
    fullfit.fit(x, y)
    fitted = fullfit.best_estimator_
    model = fitted.named_steps["model"]
    coef_raw = np.asarray(model.coef_, dtype=float).reshape(-1)
    intercept_raw = float(np.asarray(model.intercept_, dtype=float).reshape(-1)[0])
    if list(model.classes_)[1] == "high":
        coef = coef_raw.astype(float)
        intercept = float(intercept_raw)
    else:
        coef = (-coef_raw).astype(float)
        intercept = float(-intercept_raw)

    coef_rows = []
    for feature_name, coef_value in zip(combo, coef.tolist()):
        coef_rows.append(
            {
                "penalty": str(penalty),
                "combo_size": int(combo_size),
                "feature_combo": combo_label,
                "feature": str(feature_name),
                "coefficient": float(coef_value),
                "abs_coefficient": float(abs(coef_value)),
                "direction": "toward_high" if coef_value >= 0.0 else "toward_low",
                "fullfit_best_c": float(fullfit.best_params_["model__C"]),
                "fullfit_best_score": float(fullfit.best_score_),
                "intercept": float(intercept),
            }
        )

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
    inner_scoring: str,
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
            inner_scoring=str(inner_scoring),
        )
        summary_rows.append(summary_row)
        pred_rows.extend(combo_pred_rows)
        coef_rows.extend(combo_coef_rows)
        cm_rows.extend(combo_cm_rows)
    return summary_rows, pred_rows, coef_rows, cm_rows


def _combo_count(n_features: int, combo_sizes: list[int]) -> int:
    total = 0
    for size in sorted(set(int(size) for size in combo_sizes)):
        if size <= n_features:
            total += int(len(list(itertools.combinations(range(n_features), size))))
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description="Exhaustive small-combo logistic-regression search for DOR susceptibility.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_including_energy/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_logistic"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--min-r-squared", type=float, default=0.1)
    parser.add_argument("--combo-sizes", type=str, default="2,3,4,5")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--inner-scoring", type=str, default="balanced_accuracy")
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
            inner_scoring=str(args.inner_scoring),
            penalty_c_map={key: list(values) for key, values in penalty_c_map.items()},
        )
        for combo, combo_size in tasks
    )

    summary_rows = [row for combo_summary, _pred, _coef, _cm in results for row in combo_summary]
    pred_rows = [row for _summary, combo_pred, _coef, _cm in results for row in combo_pred]
    coef_rows = [row for _summary, _pred, combo_coef, _cm in results for row in combo_coef]
    cm_rows = [row for _summary, _pred, _coef, combo_cm in results for row in combo_cm]

    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[False, False, False, False, False],
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
                "inner_scoring": str(args.inner_scoring),
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
