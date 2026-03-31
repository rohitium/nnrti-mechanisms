#!/usr/bin/env python3
"""Fit one occupancy_mean-only logistic model on curated controls and score held-out mutations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import pearsonr, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_dor_susceptibility_bars import NEGATIVE_CONTROLS, POSITIVE_CONTROLS, UNCERTAIN_LIMITED


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit a single occupancy_mean logistic regression model on curated controls.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/contact_occupancy_feature_screen/tables/occupancy_mean_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/occupancy_mean_logistic_single_model"),
    )
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--c-value", type=float, default=1.0)
    parser.add_argument("--penalty", type=str, default="l2", choices=["l2"])
    parser.add_argument("--solver", type=str, default="lbfgs")
    parser.add_argument("--max-iter", type=int, default=5000)
    parser.add_argument("--decision-threshold", type=float, default=0.5)
    return parser.parse_args()


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
    raise ValueError(f"Unsupported mutation category: {label}")


def _binary_class(fold: float, low_max_fold: float) -> str:
    return "low" if float(fold) < float(low_max_fold) else "high"


def _pipeline(c_value: float, solver: str, max_iter: int) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="l2",
                    C=float(c_value),
                    solver=str(solver),
                    max_iter=int(max_iter),
                    random_state=0,
                ),
            ),
        ]
    )


def _safe_logit(prob: np.ndarray) -> np.ndarray:
    return logit(np.clip(prob.astype(float), 1e-6, 1.0 - 1e-6))


def main() -> int:
    args = _parse_args()
    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)

    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.feature_matrix_csv).copy()
    df["mutation"] = df["mutation"].astype(str)
    df["control_category"] = df["mutation"].map(_control_category)
    df["target_binary_class"] = df["dor_fold_reduction"].map(lambda x: _binary_class(float(x), float(args.low_max_fold)))

    feature_cols = [c for c in df.columns if c.startswith("occupancy_mean_")]
    if not feature_cols:
        raise ValueError("No occupancy_mean features found in feature matrix.")

    train_df = df[df["control_category"].isin({"negative_control", "positive_control"})].copy().reset_index(drop=True)
    holdout_df = df[df["control_category"].isin({"uncertain_limited", "wt_reference"})].copy().reset_index(drop=True)
    train_df.to_csv(out_tables / "control_training_feature_matrix.csv", index=False)
    holdout_df.to_csv(out_tables / "heldout_feature_matrix.csv", index=False)

    x_train = train_df[feature_cols].copy()
    y_train = train_df["target_binary_class"].astype(str).copy()
    loo = LeaveOneOut()

    cv_rows: list[dict[str, object]] = []
    prob_high = np.full(len(train_df), np.nan, dtype=float)
    for train_idx, test_idx in loo.split(x_train, y_train):
        fitted = _pipeline(args.c_value, args.solver, args.max_iter)
        fitted.fit(x_train.iloc[train_idx], y_train.iloc[train_idx])
        classes = list(fitted.named_steps["model"].classes_)
        pos_idx = int(classes.index("high"))
        probs = fitted.predict_proba(x_train.iloc[test_idx])[:, pos_idx]
        preds = np.where(probs >= float(args.decision_threshold), "high", "low")
        for row_idx, p_hi, pred in zip(test_idx.tolist(), probs.tolist(), preds.tolist()):
            prob_high[row_idx] = float(p_hi)
            cv_rows.append(
                {
                    "mutation": str(train_df.loc[row_idx, "mutation"]),
                    "control_category": str(train_df.loc[row_idx, "control_category"]),
                    "target_fold_reduction": float(train_df.loc[row_idx, "dor_fold_reduction"]),
                    "observed_class": str(train_df.loc[row_idx, "target_binary_class"]),
                    "prob_high": float(p_hi),
                    "prob_low": float(1.0 - p_hi),
                    "predicted_class": str(pred),
                    "decision_threshold": float(args.decision_threshold),
                }
            )

    cv_pred = pd.DataFrame(cv_rows).sort_values("target_fold_reduction", kind="stable").reset_index(drop=True)
    cv_pred["logit_prob_high"] = _safe_logit(cv_pred["prob_high"].to_numpy(dtype=float))
    cv_pred.to_csv(out_tables / "cv_predictions.csv", index=False)

    y_true = train_df["target_binary_class"].to_numpy(dtype=str)
    y_pred = np.where(prob_high >= float(args.decision_threshold), "high", "low")
    cm = confusion_matrix(y_true, y_pred, labels=["low", "high"])
    tn, fp, fn, tp = [int(v) for v in cm.ravel()]
    metrics = {
        "n_training_mutations": int(len(train_df)),
        "n_features": int(len(feature_cols)),
        "cv_scheme": "leave_one_out",
        "penalty": str(args.penalty),
        "solver": str(args.solver),
        "c_value": float(args.c_value),
        "decision_threshold": float(args.decision_threshold),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro")),
        "roc_auc": float(roc_auc_score((y_true == "high").astype(int), prob_high)),
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }
    logits = _safe_logit(prob_high)
    fold = train_df["dor_fold_reduction"].to_numpy(dtype=float)
    log10_fold = np.log10(fold)
    metrics["pearson_r_logit_vs_fold"] = float(pearsonr(logits, fold).statistic)
    metrics["pearson_p_logit_vs_fold"] = float(pearsonr(logits, fold).pvalue)
    metrics["pearson_r_logit_vs_log10_fold"] = float(pearsonr(logits, log10_fold).statistic)
    metrics["pearson_p_logit_vs_log10_fold"] = float(pearsonr(logits, log10_fold).pvalue)
    metrics["spearman_rho_logit_vs_fold"] = float(spearmanr(logits, fold).statistic)
    metrics["spearman_p_logit_vs_fold"] = float(spearmanr(logits, fold).pvalue)
    pd.DataFrame([metrics]).to_csv(out_tables / "model_summary.csv", index=False)

    full_fit = _pipeline(args.c_value, args.solver, args.max_iter)
    full_fit.fit(x_train, y_train)
    coefs = full_fit.named_steps["model"].coef_[0]
    coef_df = pd.DataFrame(
        {
            "feature": feature_cols,
            "coefficient": coefs,
            "abs_coefficient": np.abs(coefs),
        }
    ).sort_values("abs_coefficient", ascending=False, kind="stable")
    coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)

    x_holdout = holdout_df[feature_cols].copy()
    classes = list(full_fit.named_steps["model"].classes_)
    pos_idx = int(classes.index("high"))
    holdout_prob = full_fit.predict_proba(x_holdout)[:, pos_idx]
    holdout_pred = np.where(holdout_prob >= float(args.decision_threshold), "high", "low")
    holdout_out = holdout_df[
        ["mutation", "control_category", "dor_fold_reduction", "target_binary_class"]
    ].copy()
    holdout_out = holdout_out.rename(
        columns={
            "dor_fold_reduction": "target_fold_reduction",
            "target_binary_class": "observed_class_by_fold_cutoff",
        }
    )
    holdout_out["prob_high"] = holdout_prob
    holdout_out["prob_low"] = 1.0 - holdout_prob
    holdout_out["logit_prob_high"] = _safe_logit(holdout_prob)
    holdout_out["predicted_class"] = holdout_pred
    holdout_out["decision_threshold"] = float(args.decision_threshold)
    holdout_out.to_csv(out_tables / "heldout_uncertain_and_wt_predictions.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "training_design": "single_logistic_fit_on_all_occupancy_mean_features_using_curated_controls_only",
                "low_max_fold": float(args.low_max_fold),
                "penalty": str(args.penalty),
                "solver": str(args.solver),
                "c_value": float(args.c_value),
                "decision_threshold": float(args.decision_threshold),
                "max_iter": int(args.max_iter),
                "n_features": int(len(feature_cols)),
                "features": feature_cols,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
