#!/usr/bin/env python3
"""Benchmark tuned regression and classification models for DOR susceptibility."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import (
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
)
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.svm import SVR, LinearSVC, SVC


def _safe_corr(fn, a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2:
        return float("nan")
    try:
        value, _ = fn(a, b)
        return float(value)
    except Exception:
        return float("nan")


def _regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
        "pearson_r": _safe_corr(stats.pearsonr, y_true, y_pred),
        "spearman_rho": _safe_corr(stats.spearmanr, y_true, y_pred),
    }


def _classification_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
    }


def _regression_strata(y: pd.Series, n_bins: int = 3) -> np.ndarray:
    q = max(2, min(int(n_bins), int(len(y))))
    return pd.qcut(pd.Series(y).rank(method="first"), q=q, labels=False, duplicates="drop").to_numpy(dtype=int)


def _resistance_class_label(value: float, *, low_max: float, high_min: float) -> str:
    v = float(value)
    if v < float(low_max):
        return "low"
    if v < float(high_min):
        return "medium"
    return "high"


def _classification_labels(y: pd.Series, *, low_max: float, high_min: float) -> pd.Series:
    return y.map(lambda v: _resistance_class_label(v, low_max=low_max, high_min=high_min))


def _binary_labels(y: pd.Series, *, low_max: float) -> pd.Series:
    return y.map(lambda v: "low" if float(v) < float(low_max) else "high")


def _selection_k_grid(n_features: int) -> list[object]:
    values: list[object] = []
    for k in [3, 5, 8, 12]:
        if int(k) < int(n_features):
            values.append(int(k))
    values.append("all")
    return values


class OrdinalThresholdClassifier(BaseEstimator, ClassifierMixin):
    """Ordinal low/medium/high classifier via cumulative binary logits."""

    def __init__(
        self,
        C: float = 1.0,
        class_weight: str | dict[str, float] | None = "balanced",
        max_iter: int = 5000,
        random_state: int = 0,
    ) -> None:
        self.C = float(C)
        self.class_weight = class_weight
        self.max_iter = int(max_iter)
        self.random_state = int(random_state)

    def fit(self, X, y):
        labels = pd.Series(y, dtype="object")
        self.classes_ = np.asarray(["low", "medium", "high"], dtype=object)
        cat = pd.Categorical(labels, categories=self.classes_.tolist(), ordered=True)
        y_codes = np.asarray(cat.codes, dtype=int)
        if np.any(y_codes < 0):
            raise ValueError("OrdinalThresholdClassifier expects only low/medium/high labels.")
        self.models_: list[LogisticRegression] = []
        for threshold in range(len(self.classes_) - 1):
            target = (y_codes > threshold).astype(int)
            model = LogisticRegression(
                C=float(self.C),
                class_weight=self.class_weight,
                max_iter=int(self.max_iter),
                random_state=int(self.random_state),
            )
            model.fit(X, target)
            self.models_.append(model)
        return self

    @property
    def coef_(self) -> np.ndarray:
        return np.vstack([np.asarray(model.coef_, dtype=float).reshape(1, -1) for model in self.models_])

    def predict_proba(self, X) -> np.ndarray:
        p_gt = np.column_stack([model.predict_proba(X)[:, 1] for model in self.models_]).astype(float)
        if p_gt.shape[1] > 1:
            p_gt = np.minimum.accumulate(p_gt, axis=1)
        probs = np.zeros((p_gt.shape[0], len(self.classes_)), dtype=float)
        probs[:, 0] = 1.0 - p_gt[:, 0]
        for idx in range(1, len(self.classes_) - 1):
            probs[:, idx] = p_gt[:, idx - 1] - p_gt[:, idx]
        probs[:, -1] = p_gt[:, -1]
        probs = np.clip(probs, 0.0, 1.0)
        row_sums = probs.sum(axis=1, keepdims=True)
        probs = np.divide(probs, np.where(row_sums == 0.0, 1.0, row_sums))
        return probs

    def predict(self, X) -> np.ndarray:
        probs = self.predict_proba(X)
        return self.classes_[np.argmax(probs, axis=1)]


def _selected_feature_names(fitted_pipeline: Pipeline, feature_names: list[str]) -> list[str]:
    selector = fitted_pipeline.named_steps.get("selector")
    if selector is None or selector == "passthrough":
        return list(feature_names)
    if hasattr(selector, "get_support"):
        mask = np.asarray(selector.get_support(), dtype=bool)
        return [feature for feature, keep in zip(feature_names, mask.tolist()) if keep]
    return list(feature_names)


def _feature_importance_from_estimator(
    model_name: str,
    estimator,
    feature_names: list[str],
    *,
    selected_feature_names: list[str] | None = None,
) -> pd.DataFrame | None:
    raw: np.ndarray | None = None
    if hasattr(estimator, "feature_importances_"):
        raw = np.asarray(estimator.feature_importances_, dtype=float)
    elif hasattr(estimator, "coef_"):
        coef = np.asarray(estimator.coef_, dtype=float)
        raw = np.mean(np.abs(coef), axis=0) if coef.ndim > 1 else np.abs(coef)
    if raw is None:
        return None
    selected = list(selected_feature_names) if selected_feature_names is not None else list(feature_names)
    if len(raw) != len(selected):
        raise ValueError(f"Importance length mismatch for {model_name}: {len(raw)} vs {len(selected)}")
    raw_by_feature = {feature: 0.0 for feature in feature_names}
    for feature, value in zip(selected, raw.tolist()):
        raw_by_feature[str(feature)] = float(value)
    raw = np.asarray([raw_by_feature[str(feature)] for feature in feature_names], dtype=float)
    total = float(np.sum(raw))
    norm = raw / total if total > 0 else np.zeros_like(raw)
    df = pd.DataFrame(
        {
            "model": str(model_name),
            "feature": feature_names,
            "importance_raw": raw.astype(float),
            "importance_norm": norm.astype(float),
        }
    ).sort_values("importance_norm", ascending=False, kind="stable")
    df["rank"] = np.arange(1, len(df) + 1, dtype=int)
    return df.reset_index(drop=True)


def _permutation_importance_rows(
    model_name: str,
    fitted_pipeline: Pipeline,
    x_eval: pd.DataFrame,
    y_eval: np.ndarray,
    *,
    scoring: str,
    random_state: int,
) -> pd.DataFrame:
    result = permutation_importance(
        fitted_pipeline,
        x_eval,
        y_eval,
        n_repeats=50,
        random_state=int(random_state),
        scoring=scoring,
    )
    raw = np.asarray(result.importances_mean, dtype=float)
    total = float(np.sum(np.abs(raw)))
    norm = np.abs(raw) / total if total > 0 else np.zeros_like(raw)
    out = pd.DataFrame(
        {
            "model": str(model_name),
            "feature": list(x_eval.columns),
            "importance_raw": raw.astype(float),
            "importance_norm": norm.astype(float),
            "importance_std": np.asarray(result.importances_std, dtype=float),
        }
    ).sort_values("importance_norm", ascending=False, kind="stable")
    out["rank"] = np.arange(1, len(out) + 1, dtype=int)
    return out.reset_index(drop=True)


def _aggregate_importance(df: pd.DataFrame, model_name: str) -> pd.DataFrame:
    out = (
        df.groupby("feature", as_index=False)
        .agg(
            importance_norm=("importance_norm", "mean"),
            importance_norm_std=("importance_norm", "std"),
            rank_mean=("rank", "mean"),
            rank_std=("rank", "std"),
            n_folds=("fold", "nunique"),
        )
        .sort_values(["importance_norm", "rank_mean"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    out.insert(0, "model", str(model_name))
    return out


def _plot_regression_predictions(df: pd.DataFrame, output_png: Path, target_label: str) -> None:
    models = [c.replace("pred_", "") for c in df.columns if c.startswith("pred_")]
    n = len(models)
    fig, axes = plt.subplots(1, n, figsize=(4.7 * n, 4.2), sharex=True, sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])
    x = df["target_value"].to_numpy(dtype=float)
    lo = float(np.nanmin(x))
    hi = float(np.nanmax(x))
    pad = max(0.25, 0.08 * (hi - lo))
    palette = ["#1f77b4", "#2a9d8f", "#d62828", "#6a4c93"]
    for ax, model_name, color in zip(axes, models, palette):
        col = f"pred_{model_name}"
        ax.scatter(df["target_value"], df[col], s=46, color=color, alpha=0.9, linewidths=0)
        for _, row in df.iterrows():
            ax.text(float(row["target_value"]), float(row[col]), str(row["mutation"]), fontsize=7, alpha=0.8)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", color="#666666", linewidth=1.0)
        ax.set_title(model_name.replace("_", " ").title())
        ax.set_xlabel(f"Observed {target_label}")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(f"CV Predicted {target_label}")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_confusion_matrix(cm: np.ndarray, labels: list[str], title: str, output_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(cm, cmap="Blues", aspect="auto")
    ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(labels)), labels=labels)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Observed")
    ax.set_title(title)
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, str(int(cm[i, j])), ha="center", va="center", color="#111111", fontsize=10)
    ax.invert_yaxis()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _regression_model_specs() -> list[dict[str, object]]:
    return [
        {
            "name": "ridge",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("model", Ridge(random_state=0)),
                ]
            ),
            "param_grid": {"model__alpha": [0.1, 1.0, 10.0, 100.0]},
            "scoring": "neg_mean_absolute_error",
        },
        {
            "name": "svr_rbf",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("model", SVR(kernel="rbf")),
                ]
            ),
            "param_grid": {
                "model__C": [0.1, 1.0, 10.0],
                "model__gamma": ["scale", 0.1, 1.0],
                "model__epsilon": [0.1, 1.0, 5.0],
            },
            "scoring": "neg_mean_absolute_error",
        },
        {
            "name": "random_forest",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    (
                        "model",
                        RandomForestRegressor(random_state=0),
                    ),
                ]
            ),
            "param_grid": {
                "model__n_estimators": [100, 200],
                "model__max_depth": [None, 3],
                "model__min_samples_leaf": [1, 2],
            },
            "scoring": "neg_mean_absolute_error",
        },
        {
            "name": "gradient_boosting",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("model", GradientBoostingRegressor(random_state=0)),
                ]
            ),
            "param_grid": {
                "model__n_estimators": [100, 200],
                "model__learning_rate": [0.03, 0.1],
                "model__max_depth": [2, 3],
                "model__subsample": [0.7, 1.0],
            },
            "scoring": "neg_mean_absolute_error",
        },
    ]


def _classification_model_specs(n_features: int) -> list[dict[str, object]]:
    selector_grid = _selection_k_grid(n_features)
    return [
        {
            "name": "logistic_regression",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("selector", SelectKBest(score_func=f_classif)),
                    (
                        "model",
                        LogisticRegression(
                            max_iter=5000,
                            class_weight="balanced",
                            random_state=0,
                        ),
                    ),
                ]
            ),
            "param_grid": {
                "selector__k": selector_grid,
                "model__C": [0.1, 1.0, 10.0],
            },
            "scoring": "balanced_accuracy",
        },
        {
            "name": "linear_svm",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("selector", SelectKBest(score_func=f_classif)),
                    ("model", LinearSVC(class_weight="balanced", dual="auto", random_state=0)),
                ]
            ),
            "param_grid": {
                "selector__k": selector_grid,
                "model__C": [0.1, 1.0, 10.0],
            },
            "scoring": "balanced_accuracy",
        },
        {
            "name": "rbf_svm",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("selector", SelectKBest(score_func=f_classif)),
                    ("model", SVC(kernel="rbf", class_weight="balanced", random_state=0)),
                ]
            ),
            "param_grid": {
                "selector__k": selector_grid,
                "model__C": [0.1, 1.0, 10.0],
                "model__gamma": ["scale", 0.1, 1.0],
            },
            "scoring": "balanced_accuracy",
        },
    ]


def _ordinal_model_specs(n_features: int) -> list[dict[str, object]]:
    selector_grid = _selection_k_grid(n_features)
    return [
        {
            "name": "ordinal_logistic",
            "pipeline": Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("selector", SelectKBest(score_func=f_classif)),
                    ("model", OrdinalThresholdClassifier(class_weight="balanced", max_iter=5000, random_state=0)),
                ]
            ),
            "param_grid": {
                "selector__k": selector_grid,
                "model__C": [0.1, 1.0, 10.0],
            },
            "scoring": "balanced_accuracy",
        }
    ]


def _run_regression(
    feat: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    output_dir: Path,
    cv_folds: int,
    random_state: int,
) -> None:
    x = feat[feature_cols].copy()
    y = feat[target_col].astype(float).copy()
    strata = _regression_strata(y, n_bins=3)
    outer = StratifiedKFold(n_splits=int(cv_folds), shuffle=True, random_state=int(random_state))

    pred_df = feat[["mutation", target_col]].copy().rename(columns={target_col: "target_value"})
    summary_rows: list[dict[str, object]] = []
    param_rows: list[dict[str, object]] = []
    importance_rows: list[pd.DataFrame] = []

    for spec in _regression_model_specs():
        preds = np.full(len(y), np.nan, dtype=float)
        for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, strata), start=1):
            x_train = x.iloc[train_idx]
            y_train = y.iloc[train_idx]
            x_test = x.iloc[test_idx]
            train_strata = strata[train_idx]
            inner_splits = max(2, min(3, int(pd.Series(train_strata).value_counts().min())))
            inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=int(random_state) + fold_idx)
            search = GridSearchCV(
                spec["pipeline"],
                spec["param_grid"],
                cv=inner.split(x_train, train_strata),
                scoring=str(spec["scoring"]),
                n_jobs=1,
                refit=True,
            )
            search.fit(x_train, y_train)
            preds[test_idx] = search.best_estimator_.predict(x_test).astype(float)
            param_rows.append(
                {
                    "task": "regression",
                    "model": str(spec["name"]),
                    "fold": int(fold_idx),
                    "best_params": json.dumps(search.best_params_, sort_keys=True),
                    "best_score": float(search.best_score_),
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                }
            )
            imp = _feature_importance_from_estimator(
                str(spec["name"]),
                search.best_estimator_.named_steps["model"],
                feature_cols,
            )
            if imp is not None:
                imp["fold"] = int(fold_idx)
                importance_rows.append(imp)
        pred_df[f"pred_{spec['name']}"] = preds
        summary_rows.append(
            {
                "task": "regression",
                "model": str(spec["name"]),
                "target": str(target_col),
                "n_mutations": int(len(y)),
                "cv_folds": int(cv_folds),
                **_regression_metrics(y.to_numpy(dtype=float), preds),
            }
        )

    pd.DataFrame(summary_rows).sort_values("mae", ascending=True).to_csv(output_dir / "tables" / "regression_cv_summary.csv", index=False)
    pred_df.to_csv(output_dir / "tables" / "regression_cv_predictions.csv", index=False)
    pd.DataFrame(param_rows).to_csv(output_dir / "tables" / "regression_best_params_by_fold.csv", index=False)
    _plot_regression_predictions(pred_df, output_dir / "plots" / "regression_cv_predictions.png", "Fold Reduction")

    if importance_rows:
        imp_df = pd.concat(importance_rows, ignore_index=True)
        for model_name, grp in imp_df.groupby("model"):
            grp.to_csv(output_dir / "tables" / f"regression_feature_importance_{model_name}_by_fold.csv", index=False)
            _aggregate_importance(grp, str(model_name)).to_csv(
                output_dir / "tables" / f"regression_feature_importance_{model_name}.csv",
                index=False,
            )


def _run_classification_task(
    feat: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    task_name: str,
    y: pd.Series,
    classes: list[str],
    model_specs: list[dict[str, object]],
    output_dir: Path,
    cv_folds: int,
    random_state: int,
    low_max: float,
    high_min: float,
) -> None:
    x = feat[feature_cols].copy()
    y = y.astype(str).copy()
    min_class_count = int(y.value_counts().min())
    effective_folds = max(2, min(int(cv_folds), int(min_class_count)))
    outer = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(random_state))

    pred_df = feat[["mutation", target_col]].copy().rename(columns={target_col: "target_value"})
    pred_df["observed_class"] = y.to_numpy(dtype=str)
    summary_rows: list[dict[str, object]] = []
    param_rows: list[dict[str, object]] = []
    importance_rows: list[pd.DataFrame] = []

    for spec in model_specs:
        preds = np.full(len(y), "", dtype=object)
        for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
            x_train = x.iloc[train_idx]
            y_train = y.iloc[train_idx]
            x_test = x.iloc[test_idx]
            inner_splits = max(2, min(3, int(y_train.value_counts().min())))
            inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=int(random_state) + fold_idx)
            search = GridSearchCV(
                spec["pipeline"],
                spec["param_grid"],
                cv=inner.split(x_train, y_train),
                scoring=str(spec["scoring"]),
                n_jobs=1,
                refit=True,
            )
            search.fit(x_train, y_train)
            preds[test_idx] = search.best_estimator_.predict(x_test).astype(str)
            param_rows.append(
                {
                    "task": str(task_name),
                    "model": str(spec["name"]),
                    "fold": int(fold_idx),
                    "best_params": json.dumps(search.best_params_, sort_keys=True),
                    "best_score": float(search.best_score_),
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                }
            )
            selected_features = _selected_feature_names(search.best_estimator_, feature_cols)
            imp = _feature_importance_from_estimator(
                str(spec["name"]),
                search.best_estimator_.named_steps["model"],
                feature_cols,
                selected_feature_names=selected_features,
            )
            if imp is None:
                imp = _permutation_importance_rows(
                    str(spec["name"]),
                    search.best_estimator_,
                    x_test,
                    y.iloc[test_idx].to_numpy(dtype=str),
                    scoring="balanced_accuracy",
                    random_state=int(random_state) + fold_idx,
                )
            imp["fold"] = int(fold_idx)
            imp["task"] = str(task_name)
            importance_rows.append(imp)

        pred_df[f"pred_{spec['name']}"] = preds
        metrics = _classification_metrics(y.to_numpy(dtype=str), preds.astype(str))
        summary_rows.append(
            {
                "task": str(task_name),
                "model": str(spec["name"]),
                "target": str(task_name),
                "n_mutations": int(len(y)),
                "cv_folds_requested": int(cv_folds),
                "cv_folds_effective": int(effective_folds),
                "low_max_fold": float(low_max),
                "high_min_fold": float(high_min),
                **metrics,
            }
        )
        cm = confusion_matrix(y, preds, labels=classes)
        cm_df = pd.DataFrame(cm, index=[f"obs_{c}" for c in classes], columns=[f"pred_{c}" for c in classes]).reset_index().rename(columns={"index": "observed"})
        cm_df.to_csv(output_dir / "tables" / f"{task_name}_confusion_matrix_{spec['name']}.csv", index=False)
        _plot_confusion_matrix(
            cm,
            classes,
            f"{task_name.replace('_', ' ').title()}: {spec['name'].replace('_', ' ').title()}",
            output_dir / "plots" / f"{task_name}_confusion_matrix_{spec['name']}.png",
        )

    pd.DataFrame(summary_rows).sort_values("balanced_accuracy", ascending=False).to_csv(
        output_dir / "tables" / f"{task_name}_cv_summary.csv",
        index=False,
    )
    pred_df.to_csv(output_dir / "tables" / f"{task_name}_cv_predictions.csv", index=False)
    pd.DataFrame(param_rows).to_csv(output_dir / "tables" / f"{task_name}_best_params_by_fold.csv", index=False)
    if importance_rows:
        imp_df = pd.concat(importance_rows, ignore_index=True)
        for model_name, grp in imp_df.groupby("model"):
            grp.to_csv(output_dir / "tables" / f"{task_name}_feature_importance_{model_name}_by_fold.csv", index=False)
            _aggregate_importance(grp, str(model_name)).to_csv(
                output_dir / "tables" / f"{task_name}_feature_importance_{model_name}.csv",
                index=False,
            )


def _run_classification_benchmarks(
    feat: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    output_dir: Path,
    cv_folds: int,
    random_state: int,
    low_max: float,
    high_min: float,
) -> None:
    y_reg = feat[target_col].astype(float).copy()
    multiclass = _classification_labels(y_reg, low_max=low_max, high_min=high_min)
    binary = _binary_labels(y_reg, low_max=low_max)
    _run_classification_task(
        feat,
        feature_cols,
        target_col,
        task_name="classification_multiclass",
        y=multiclass,
        classes=["low", "medium", "high"],
        model_specs=_classification_model_specs(len(feature_cols)),
        output_dir=output_dir,
        cv_folds=int(cv_folds),
        random_state=int(random_state),
        low_max=float(low_max),
        high_min=float(high_min),
    )
    _run_classification_task(
        feat,
        feature_cols,
        target_col,
        task_name="classification_binary",
        y=binary,
        classes=["low", "high"],
        model_specs=_classification_model_specs(len(feature_cols)),
        output_dir=output_dir,
        cv_folds=int(cv_folds),
        random_state=int(random_state),
        low_max=float(low_max),
        high_min=float(high_min),
    )
    _run_classification_task(
        feat,
        feature_cols,
        target_col,
        task_name="classification_ordinal",
        y=multiclass,
        classes=["low", "medium", "high"],
        model_specs=_ordinal_model_specs(len(feature_cols)),
        output_dir=output_dir,
        cv_folds=int(cv_folds),
        random_state=int(random_state),
        low_max=float(low_max),
        high_min=float(high_min),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark tuned susceptibility models from mutation-level features.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/susceptibility_ml/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/resistance_model_benchmark"),
    )
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--high-min-fold", type=float, default=50.0)
    parser.add_argument("--skip-regression", action="store_true")
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    feat = pd.read_csv(args.feature_matrix_csv)
    target_col = "target_fold_reduction"
    non_feature_cols = {
        "drug",
        "mutation",
        "chain",
        "dor_fold_reduction",
        "order",
        "target_fold_reduction",
        "target_log10_fold_reduction",
        "target_ddg_exp_kj",
    }
    feature_cols = [c for c in feat.columns if c not in non_feature_cols]

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    if not args.skip_regression:
        _run_regression(
            feat,
            feature_cols,
            target_col,
            output_dir=args.output_dir,
            cv_folds=int(args.cv_folds),
            random_state=int(args.random_state),
        )
    _run_classification_benchmarks(
        feat,
        feature_cols,
        target_col,
        output_dir=args.output_dir,
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        low_max=float(args.low_max_fold),
        high_min=float(args.high_min_fold),
    )

    config = {
        "feature_matrix_csv": str(args.feature_matrix_csv),
        "target_col": str(target_col),
        "cv_folds_requested": int(args.cv_folds),
        "random_state": int(args.random_state),
        "low_max_fold": float(args.low_max_fold),
        "high_min_fold": float(args.high_min_fold),
        "n_mutations": int(len(feat)),
        "n_features": int(len(feature_cols)),
        "feature_columns": feature_cols,
        "feature_selection_inside_cv": {
            "selector": "SelectKBest",
            "score_func": "f_classif",
            "k_grid": _selection_k_grid(len(feature_cols)),
        },
        "classification_tasks": [
            "classification_multiclass",
            "classification_binary",
            "classification_ordinal",
        ],
        "skip_regression": bool(args.skip_regression),
        "xgboost_available": False,
    }
    (out_config / "run_config.json").write_text(json.dumps(config, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
