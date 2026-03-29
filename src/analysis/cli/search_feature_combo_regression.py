#!/usr/bin/env python3
"""Screen correlated features and exhaustively benchmark small regression combos."""
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
from sklearn.linear_model import Lasso, Ridge
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .benchmark_resistance_models import _regression_metrics, _regression_strata


META_COLUMNS = {"drug", "mutation", "chain", "target_fold_reduction", "target_binary_class"}


def _parse_int_list(text: str) -> list[int]:
    values: list[int] = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        values.append(int(token))
    if not values:
        raise ValueError(f"Expected at least one integer in {text!r}")
    return values


def _parse_float_list(text: str) -> list[float]:
    values: list[float] = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    if not values:
        raise ValueError(f"Expected at least one float in {text!r}")
    return values


def _transform_target(y: pd.Series, *, mode: str) -> pd.Series:
    values = pd.to_numeric(y, errors="coerce").astype(float)
    if str(mode) == "raw":
        return values
    if np.any(values <= 0.0):
        raise ValueError("log10 target transform requires strictly positive fold values")
    return np.log10(values)


def _feature_target_correlations(
    feat: pd.DataFrame,
    *,
    feature_cols: list[str],
    target_col: str,
    target_transform: str,
) -> pd.DataFrame:
    y_raw = pd.to_numeric(feat[target_col], errors="coerce").astype(float)
    y = _transform_target(y_raw, mode=target_transform)
    rows: list[dict[str, object]] = []
    for feature in feature_cols:
        x = pd.to_numeric(feat[feature], errors="coerce").astype(float)
        sub = pd.DataFrame({"feature_value": x, "target_value": y, "target_value_raw": y_raw}).dropna().copy()
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
                "target_transform": str(target_transform),
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


def _regression_pipeline(model_name: str) -> Pipeline:
    if str(model_name) == "ridge":
        model = Ridge(random_state=0)
    elif str(model_name) == "lasso":
        model = Lasso(max_iter=20000, random_state=0)
    else:
        raise ValueError(f"Unsupported model: {model_name}")
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def _grid_for_model(model_name: str, alphas: list[float]) -> dict[str, list[float]]:
    return {"model__alpha": [float(alpha) for alpha in alphas]}


def _evaluate_combo_model(
    feat: pd.DataFrame,
    *,
    combo: tuple[str, ...],
    combo_size: int,
    feature_count_filtered: int,
    target_col: str,
    target_transform: str,
    cv_folds: int,
    random_state: int,
    model_name: str,
    alphas: list[float],
    inner_scoring: str,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    combo_label = "|".join(combo)
    y_raw = pd.to_numeric(feat[target_col], errors="coerce").astype(float)
    y = _transform_target(y_raw, mode=target_transform).astype(float)
    x = feat[list(combo)].copy()
    effective_folds = max(2, min(int(cv_folds), int(len(y))))
    strata = _regression_strata(pd.Series(y), n_bins=3)
    outer = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(random_state))

    preds = np.full(len(y), np.nan, dtype=float)
    chosen_alphas: list[float] = []
    best_scores: list[float] = []
    prediction_rows: list[dict[str, object]] = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, strata), start=1):
        x_train = x.iloc[train_idx]
        y_train = y.iloc[train_idx]
        x_test = x.iloc[test_idx]
        train_strata = strata[train_idx]
        inner_splits = max(2, min(3, int(pd.Series(train_strata).value_counts().min())))
        inner = StratifiedKFold(
            n_splits=inner_splits,
            shuffle=True,
            random_state=int(random_state) + hash((combo_label, model_name, fold_idx)) % 100000,
        )
        search = GridSearchCV(
            _regression_pipeline(model_name),
            _grid_for_model(model_name, alphas),
            cv=inner.split(x_train, train_strata),
            scoring=str(inner_scoring),
            n_jobs=1,
            refit=True,
        )
        search.fit(x_train, y_train)
        preds[test_idx] = search.best_estimator_.predict(x_test).astype(float)
        chosen_alphas.append(float(search.best_params_["model__alpha"]))
        best_scores.append(float(search.best_score_))
        for sample_idx in test_idx:
            prediction_rows.append(
                {
                    "model": str(model_name),
                    "combo_size": int(combo_size),
                    "feature_combo": combo_label,
                    "fold": int(fold_idx),
                    "mutation": str(feat.iloc[sample_idx]["mutation"]),
                    "target_value_raw": float(y_raw.iloc[sample_idx]),
                    "target_value_model": float(y.iloc[sample_idx]),
                    "predicted_value_model": float(preds[sample_idx]),
                }
            )

    metrics = _regression_metrics(y.to_numpy(dtype=float), preds)
    summary_row = {
        "model": str(model_name),
        "combo_size": int(combo_size),
        "feature_combo": combo_label,
        "n_features_filtered": int(feature_count_filtered),
        "cv_folds": int(effective_folds),
        "target_transform": str(target_transform),
        "inner_scoring": str(inner_scoring),
        "alpha_by_fold": json.dumps(chosen_alphas),
        "alpha_mean": float(np.mean(chosen_alphas)),
        "alpha_median": float(np.median(chosen_alphas)),
        "inner_best_score_mean": float(np.mean(best_scores)),
        "inner_best_score_std": float(np.std(best_scores, ddof=1)) if len(best_scores) > 1 else 0.0,
        **metrics,
    }

    fullfit = GridSearchCV(
        _regression_pipeline(model_name),
        _grid_for_model(model_name, alphas),
        cv=outer.split(x, strata),
        scoring=str(inner_scoring),
        n_jobs=1,
        refit=True,
    )
    fullfit.fit(x, y)
    estimator = fullfit.best_estimator_
    coef = np.asarray(estimator.named_steps["model"].coef_, dtype=float).reshape(-1)
    intercept = float(np.asarray(estimator.named_steps["model"].intercept_, dtype=float).reshape(-1)[0])
    coef_rows = [
        {
            "model": str(model_name),
            "combo_size": int(combo_size),
            "feature_combo": combo_label,
            "feature": str(feature_name),
            "coefficient": float(coef_value),
            "abs_coefficient": float(abs(coef_value)),
            "fullfit_best_alpha": float(fullfit.best_params_["model__alpha"]),
            "fullfit_best_score": float(fullfit.best_score_),
            "intercept": float(intercept),
            "target_transform": str(target_transform),
        }
        for feature_name, coef_value in zip(combo, coef.tolist())
    ]
    return summary_row, prediction_rows, coef_rows


def _evaluate_combo(
    feat: pd.DataFrame,
    *,
    combo: tuple[str, ...],
    combo_size: int,
    feature_count_filtered: int,
    target_col: str,
    target_transform: str,
    cv_folds: int,
    random_state: int,
    inner_scoring: str,
    model_alpha_map: dict[str, list[float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    summary_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    coef_rows: list[dict[str, object]] = []
    for model_name, alphas in model_alpha_map.items():
        summary_row, pred_rows, combo_coef_rows = _evaluate_combo_model(
            feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(feature_count_filtered),
            target_col=str(target_col),
            target_transform=str(target_transform),
            cv_folds=int(cv_folds),
            random_state=int(random_state),
            model_name=str(model_name),
            alphas=list(alphas),
            inner_scoring=str(inner_scoring),
        )
        summary_rows.append(summary_row)
        prediction_rows.extend(pred_rows)
        coef_rows.extend(combo_coef_rows)
    return summary_rows, prediction_rows, coef_rows


def _nested_cv_combo_search(
    feat: pd.DataFrame,
    *,
    feature_cols: list[str],
    target_col: str,
    target_transform: str,
    combo_sizes: list[int],
    min_r_squared: float,
    cv_folds: int,
    random_state: int,
    ridge_alphas: list[float],
    lasso_alphas: list[float],
    inner_scoring: str,
    n_jobs: int,
    parallel_backend: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    correlations = _feature_target_correlations(
        feat,
        feature_cols=feature_cols,
        target_col=target_col,
        target_transform=target_transform,
    )
    filtered_features = (
        correlations.loc[correlations["r_squared"] >= float(min_r_squared), "feature"]
        .astype(str)
        .tolist()
    )

    model_alpha_map = {
        "ridge": list(ridge_alphas),
        "lasso": list(lasso_alphas),
    }

    tasks: list[tuple[tuple[str, ...], int]] = []
    for combo_size in sorted(set(int(size) for size in combo_sizes)):
        if combo_size <= 0 or combo_size > len(filtered_features):
            continue
        for combo in itertools.combinations(filtered_features, combo_size):
            tasks.append((combo, combo_size))

    parallel_kwargs = {"n_jobs": int(n_jobs), "verbose": 10}
    if str(parallel_backend) == "threading":
        parallel_kwargs["prefer"] = "threads"
    results = Parallel(**parallel_kwargs)(
        delayed(_evaluate_combo)(
            feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(len(filtered_features)),
            target_col=str(target_col),
            target_transform=str(target_transform),
            cv_folds=int(cv_folds),
            random_state=int(random_state),
            inner_scoring=str(inner_scoring),
            model_alpha_map={k: list(v) for k, v in model_alpha_map.items()},
        )
        for combo, combo_size in tasks
    )

    combo_rows = [row for summary_rows, _pred_rows, _coef_rows in results for row in summary_rows]
    prediction_rows = [row for _summary_rows, pred_rows, _coef_rows in results for row in pred_rows]
    coef_rows = [row for _summary_rows, _pred_rows, coef_rows in results for row in coef_rows]

    filtered_df = correlations[correlations["feature"].isin(filtered_features)].copy().reset_index(drop=True)
    combo_df = pd.DataFrame(combo_rows).sort_values(
        ["r2", "pearson_r", "spearman_rho", "rmse"],
        ascending=[False, False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    pred_df = pd.DataFrame(prediction_rows)
    coef_df = pd.DataFrame(coef_rows).sort_values(
        ["model", "combo_size", "feature_combo", "abs_coefficient"],
        ascending=[True, True, True, False],
        kind="stable",
    ).reset_index(drop=True)
    return correlations, filtered_df, combo_df, pred_df, coef_df


def main() -> int:
    parser = argparse.ArgumentParser(description="Exhaustive small-combo Ridge/Lasso search for susceptibility regression.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_including_energy/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_regression"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--target-transform", type=str, default="raw", choices=["raw", "log10"])
    parser.add_argument("--min-r-squared", type=float, default=0.1)
    parser.add_argument("--combo-sizes", type=str, default="2,3,4,5")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--inner-scoring", type=str, default="r2")
    parser.add_argument(
        "--ridge-alphas",
        type=str,
        default="0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100",
    )
    parser.add_argument(
        "--lasso-alphas",
        type=str,
        default="0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3,1",
    )
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--parallel-backend", type=str, default="threading", choices=["threading", "loky"])
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    combo_sizes = _parse_int_list(args.combo_sizes)
    ridge_alphas = _parse_float_list(args.ridge_alphas)
    lasso_alphas = _parse_float_list(args.lasso_alphas)

    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feat = pd.read_csv(args.feature_matrix_csv)
    feature_cols = [c for c in feat.columns if c not in META_COLUMNS]

    corr_df, filtered_df, combo_df, pred_df, coef_df = _nested_cv_combo_search(
        feat,
        feature_cols=feature_cols,
        target_col=str(args.target_col),
        target_transform=str(args.target_transform),
        combo_sizes=combo_sizes,
        min_r_squared=float(args.min_r_squared),
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        ridge_alphas=ridge_alphas,
        lasso_alphas=lasso_alphas,
        inner_scoring=str(args.inner_scoring),
        n_jobs=int(args.n_jobs),
        parallel_backend=str(args.parallel_backend),
    )

    corr_df.to_csv(out_tables / "feature_target_correlations.csv", index=False)
    filtered_df.to_csv(out_tables / "filtered_features.csv", index=False)
    combo_df.to_csv(out_tables / "combo_model_summary.csv", index=False)
    pred_df.to_csv(out_tables / "combo_cv_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "combo_fullfit_coefficients.csv", index=False)

    top_rows: list[pd.DataFrame] = []
    if not combo_df.empty:
        for model_name in sorted(combo_df["model"].astype(str).unique()):
            for combo_size in sorted(combo_df["combo_size"].astype(int).unique()):
                sub = combo_df[(combo_df["model"].astype(str) == model_name) & (combo_df["combo_size"].astype(int) == combo_size)]
                if sub.empty:
                    continue
                top_rows.append(sub.head(10).copy())
    top_df = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()
    top_df.to_csv(out_tables / "top_combo_models.csv", index=False)

    total_combos = int(
        sum(
            len(list(itertools.combinations(filtered_df["feature"].astype(str).tolist(), size)))
            for size in sorted(set(combo_sizes))
            if size <= len(filtered_df)
        )
    )
    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "target_col": str(args.target_col),
                "target_transform": str(args.target_transform),
                "min_r_squared": float(args.min_r_squared),
                "combo_sizes": combo_sizes,
                "cv_folds": int(args.cv_folds),
                "random_state": int(args.random_state),
                "inner_scoring": str(args.inner_scoring),
                "ridge_alphas": ridge_alphas,
                "lasso_alphas": lasso_alphas,
                "n_mutations": int(len(feat)),
                "n_features_input": int(len(feature_cols)),
                "n_features_filtered": int(len(filtered_df)),
                "filtered_features": filtered_df["feature"].astype(str).tolist(),
                "n_feature_combinations": total_combos,
                "models": ["ridge", "lasso"],
                "n_jobs": int(args.n_jobs),
                "parallel_backend": str(args.parallel_backend),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
