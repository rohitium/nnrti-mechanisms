#!/usr/bin/env python3
"""Fit and plot a fixed parsimonious logistic DOR classifier."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from .benchmark_resistance_models import _plot_confusion_matrix
from .focus_binary_logistic_classifier import (
    _plot_cv_probability_ranked,
    _plot_probability_vs_fold,
    _plot_selected_coefficients,
)
from .model_susceptibility_from_state_features import _mutation_feature_matrix
from .search_feature_combo_logistic import _binary_labels, _classification_metrics, _logistic_pipeline
from ..susceptibility import load_dor_susceptibilities


FEATURES = [
    "residue_min_distance_LYS101_angstrom_repstd",
    "ligand_pose_rmsd_angstrom_mean",
    "ligand_palm_distance_angstrom_repstd",
]


def _wt_feature_row(
    *,
    frame_feature_csv: Path,
    mmgbsa_replicate_csv: Path,
    susceptibility_xlsx: Path,
) -> pd.DataFrame:
    frame_df = pd.read_csv(frame_feature_csv)
    mmgbsa_df = pd.read_csv(mmgbsa_replicate_csv) if mmgbsa_replicate_csv.exists() else None
    target_df = load_dor_susceptibilities(susceptibility_xlsx)
    wt_row = pd.DataFrame(
        [{"drug": "DOR", "mutation": "WT", "chain": "A", "dor_fold_reduction": 1.0, "order": -1}]
    )
    target_aug = pd.concat([target_df, wt_row], ignore_index=True)
    feat = _mutation_feature_matrix(
        frame_df,
        target_df=target_aug,
        temperature_k=300.0,
        dispersion_mode="replicate_sd",
        mmgbsa_df=mmgbsa_df,
    )
    return feat[feat["mutation"].astype(str) == "WT"].copy()


def main() -> int:
    parser = argparse.ArgumentParser(description="Run fixed parsimonious logistic regression diagnostics.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_including_energy/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression_parsimonious"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--penalty", type=str, default="l1", choices=["l1", "l2"])
    parser.add_argument("--c-value", type=float, default=1.0)
    parser.add_argument("--decision-threshold", type=float, default=0.4)
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
    args = parser.parse_args()

    feat = pd.read_csv(args.feature_matrix_csv)
    x = feat[FEATURES].copy()
    y_value = pd.to_numeric(feat[str(args.target_col)], errors="coerce").astype(float)
    y = _binary_labels(y_value, low_max=float(args.low_max_fold))

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    effective_folds = max(2, min(int(args.cv_folds), int(y.value_counts().min())))
    outer = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(args.random_state))

    prob_high = np.full(len(y), np.nan, dtype=float)
    preds = np.empty(len(y), dtype=object)
    pred_rows: list[dict[str, object]] = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
        pipeline = _logistic_pipeline(random_state=int(args.random_state), penalty=str(args.penalty))
        pipeline.set_params(model__C=float(args.c_value))
        pipeline.fit(x.iloc[train_idx], y.iloc[train_idx])
        classes = list(pipeline.named_steps["model"].classes_)
        pos_idx = int(classes.index("high"))
        fold_prob = pipeline.predict_proba(x.iloc[test_idx])[:, pos_idx].astype(float)
        fold_pred = np.where(fold_prob >= float(args.decision_threshold), "high", "low")
        prob_high[test_idx] = fold_prob
        preds[test_idx] = fold_pred
        for sample_idx, sample_prob, sample_pred in zip(test_idx, fold_prob.tolist(), fold_pred.tolist()):
            pred_rows.append(
                {
                    "fold": int(fold_idx),
                    "mutation": str(feat.iloc[sample_idx]["mutation"]),
                    "target_value": float(y_value.iloc[sample_idx]),
                    "observed_class": str(y.iloc[sample_idx]),
                    "predicted_class": str(sample_pred),
                    "prob_high": float(sample_prob),
                    "prob_low": float(1.0 - sample_prob),
                    "decision_threshold": float(args.decision_threshold),
                }
            )

    pred_df = pd.DataFrame(pred_rows).sort_values(["prob_high", "target_value", "mutation"], ascending=[True, True, True], kind="stable").reset_index(drop=True)
    y_np = y.to_numpy(dtype=str)
    y_bin = (y == "high").astype(int).to_numpy(dtype=int)
    metrics = _classification_metrics(y_np, preds.astype(str))
    metrics["roc_auc"] = float(roc_auc_score(y_bin, prob_high))
    metrics["average_precision"] = float(average_precision_score(y_bin, prob_high))
    cm = confusion_matrix(y_np, preds.astype(str), labels=["low", "high"])

    fitted = _logistic_pipeline(random_state=int(args.random_state), penalty=str(args.penalty))
    fitted.set_params(model__C=float(args.c_value))
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
    coef_df = pd.DataFrame(
        [
            {
                "feature": str(feature_name),
                "coefficient": float(coef_value),
                "abs_coefficient": float(abs(coef_value)),
                "direction": "toward_high" if coef_value >= 0.0 else "toward_low",
                "intercept": float(intercept),
                "penalty": str(args.penalty),
                "c_value": float(args.c_value),
                "decision_threshold": float(args.decision_threshold),
            }
            for feature_name, coef_value in zip(FEATURES, coef.tolist())
        ]
    ).sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)

    cm_df = pd.DataFrame(cm, index=["obs_low", "obs_high"], columns=["pred_low", "pred_high"]).reset_index().rename(columns={"index": "observed"})
    summary_row = {
        "model_name": "parsimonious_fixed_threshold_logistic",
        "feature_combo": "|".join(FEATURES),
        "penalty": str(args.penalty),
        "c_value": float(args.c_value),
        "decision_threshold": float(args.decision_threshold),
        "cv_folds": int(effective_folds),
        "low_max_fold": float(args.low_max_fold),
        "tn": int(cm[0, 0]),
        "fp": int(cm[0, 1]),
        "fn": int(cm[1, 0]),
        "tp": int(cm[1, 1]),
        **metrics,
    }

    pred_df.to_csv(out_tables / "cv_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)
    cm_df.to_csv(out_tables / "confusion_matrix.csv", index=False)
    pd.DataFrame([summary_row]).to_csv(out_tables / "model_summary.csv", index=False)

    _plot_cv_probability_ranked(pred_df, out_plots / "cv_probability_ranked.png")
    _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_fold.png", use_log10_x=False)
    _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_log10_fold.png", use_log10_x=True)
    _plot_selected_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png")
    _plot_confusion_matrix(
        cm,
        ["low", "high"],
        f"Parsimonious Logistic Confusion Matrix\nThreshold = {float(args.decision_threshold):.3f}",
        out_plots / "confusion_matrix.png",
    )

    imputer = fitted.named_steps["imputer"]
    scaler = fitted.named_steps["scaler"]
    raw_intercept = float(intercept)
    for mean_value, scale_value, coef_value in zip(scaler.mean_, scaler.scale_, coef):
        raw_intercept -= float(coef_value) * float(mean_value) / float(scale_value)
    raw_coef_df = pd.DataFrame(
        [
            {
                "feature": str(feature_name),
                "median_impute_value": float(median_value),
                "scaler_mean": float(mean_value),
                "scaler_scale": float(scale_value),
                "standardized_coefficient": float(coef_value),
                "raw_space_coefficient": float(coef_value) / float(scale_value),
            }
            for feature_name, median_value, mean_value, scale_value, coef_value in zip(
                FEATURES,
                imputer.statistics_,
                scaler.mean_,
                scaler.scale_,
                coef.tolist(),
            )
        ]
    )
    raw_coef_df.to_csv(out_tables / "equation_terms.csv", index=False)

    wt_df = _wt_feature_row(
        frame_feature_csv=args.frame_feature_csv,
        mmgbsa_replicate_csv=args.mmgbsa_replicate_csv,
        susceptibility_xlsx=args.susceptibility_xlsx,
    )
    if not wt_df.empty:
        x_wt = wt_df[FEATURES].copy()
        classes = list(fitted.named_steps["model"].classes_)
        pos_idx = int(classes.index("high"))
        wt_prob = float(fitted.predict_proba(x_wt)[:, pos_idx][0])
        wt_pred = "high" if wt_prob >= float(args.decision_threshold) else "low"
        wt_row = {
            "fold": 0,
            "mutation": "WT",
            "target_value": 1.0,
            "observed_class": "low",
            "predicted_class": wt_pred,
            "prob_high": wt_prob,
            "prob_low": float(1.0 - wt_prob),
            "decision_threshold": float(args.decision_threshold),
            "prediction_source": "full_fit_reference",
        }
        pred_with_wt = pred_df.copy()
        pred_with_wt["prediction_source"] = "cv_out_of_fold"
        pred_with_wt = pd.concat([pred_with_wt, pd.DataFrame([wt_row])], ignore_index=True)
        pred_with_wt = pred_with_wt.sort_values(["prob_high", "target_value", "mutation"], ascending=[True, True, True], kind="stable").reset_index(drop=True)
        pred_with_wt.to_csv(out_tables / "predictions_with_wt.csv", index=False)
        _plot_cv_probability_ranked(pred_with_wt, out_plots / "probability_ranked_with_wt.png")
        _plot_probability_vs_fold(pred_with_wt, out_plots / "probability_vs_fold_with_wt.png", use_log10_x=False)
        _plot_probability_vs_fold(pred_with_wt, out_plots / "probability_vs_log10_fold_with_wt.png", use_log10_x=True)
        cm_with_wt = confusion_matrix(
            pred_with_wt["observed_class"],
            pred_with_wt["predicted_class"],
            labels=["low", "high"],
        )
        cm_with_wt_df = (
            pd.DataFrame(cm_with_wt, index=["obs_low", "obs_high"], columns=["pred_low", "pred_high"])
            .reset_index()
            .rename(columns={"index": "observed"})
        )
        cm_with_wt_df.to_csv(out_tables / "confusion_matrix_with_wt.csv", index=False)
        _plot_confusion_matrix(
            cm_with_wt,
            ["low", "high"],
            f"Parsimonious Logistic Confusion Matrix With WT\n19 CV predictions + 1 full-fit WT reference, Threshold = {float(args.decision_threshold):.3f}",
            out_plots / "confusion_matrix_with_wt.png",
        )

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "features": FEATURES,
                "target_col": str(args.target_col),
                "low_max_fold": float(args.low_max_fold),
                "cv_folds": int(effective_folds),
                "random_state": int(args.random_state),
                "penalty": str(args.penalty),
                "c_value": float(args.c_value),
                "decision_threshold": float(args.decision_threshold),
                "frame_feature_csv": str(args.frame_feature_csv),
                "mmgbsa_replicate_csv": str(args.mmgbsa_replicate_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "fullfit_intercept_standardized_space": float(intercept),
                "fullfit_intercept_raw_space": float(raw_intercept),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
