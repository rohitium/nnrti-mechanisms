#!/usr/bin/env python3
"""Rank control-trained logistic models by how well OOF logits track susceptibility."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .model_susceptibility_from_state_features import _mutation_feature_matrix
from ..susceptibility import load_dor_susceptibilities


def _safe_corr(x: pd.Series, y: pd.Series) -> dict[str, float]:
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return {
            "pearson_r": float("nan"),
            "pearson_pvalue": float("nan"),
            "spearman_rho": float("nan"),
            "spearman_pvalue": float("nan"),
        }
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)
    return {
        "pearson_r": float(pearson_r),
        "pearson_pvalue": float(pearson_p),
        "spearman_rho": float(spearman_rho),
        "spearman_pvalue": float(spearman_p),
    }


def _model_correlations(pred_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (penalty, combo_size, feature_combo), group in pred_df.groupby(
        ["penalty", "combo_size", "feature_combo"], sort=False
    ):
        y = pd.to_numeric(group["target_value"], errors="coerce").astype(float)
        p = pd.to_numeric(group["prob_high"], errors="coerce").astype(float)
        p_clip = p.clip(1e-6, 1 - 1e-6)
        logit = (p_clip / (1.0 - p_clip)).map(math.log)
        row = {
            "penalty": str(penalty),
            "combo_size": int(combo_size),
            "feature_combo": str(feature_combo),
            "n_controls": int(len(group)),
            "logit_raw_pearson_r": _safe_corr(logit, y)["pearson_r"],
            "logit_raw_pearson_pvalue": _safe_corr(logit, y)["pearson_pvalue"],
            "logit_raw_spearman_rho": _safe_corr(logit, y)["spearman_rho"],
            "logit_raw_spearman_pvalue": _safe_corr(logit, y)["spearman_pvalue"],
            "logit_log10_pearson_r": _safe_corr(logit, y.map(math.log10))["pearson_r"],
            "logit_log10_pearson_pvalue": _safe_corr(logit, y.map(math.log10))["pearson_pvalue"],
            "logit_log10_spearman_rho": _safe_corr(logit, y.map(math.log10))["spearman_rho"],
            "logit_log10_spearman_pvalue": _safe_corr(logit, y.map(math.log10))["spearman_pvalue"],
        }
        for class_name in ("low", "high"):
            sub = group[group["observed_class"].astype(str) == class_name].copy()
            y_sub = pd.to_numeric(sub["target_value"], errors="coerce").astype(float)
            p_sub = pd.to_numeric(sub["prob_high"], errors="coerce").astype(float).clip(1e-6, 1 - 1e-6)
            logit_sub = (p_sub / (1.0 - p_sub)).map(math.log)
            corr = _safe_corr(logit_sub, y_sub)
            row[f"{class_name}_n"] = int(len(sub))
            row[f"{class_name}_logit_raw_pearson_r"] = corr["pearson_r"]
            row[f"{class_name}_logit_raw_pearson_pvalue"] = corr["pearson_pvalue"]
            row[f"{class_name}_logit_raw_spearman_rho"] = corr["spearman_rho"]
            row[f"{class_name}_logit_raw_spearman_pvalue"] = corr["spearman_pvalue"]
        rows.append(row)
    return pd.DataFrame(rows)


def _logistic_pipeline(*, penalty: str, c_value: float) -> Pipeline:
    penalty_text = str(penalty).strip().lower()
    solver = "liblinear" if penalty_text == "l1" else "lbfgs"
    pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=0,
                    penalty=penalty_text,
                    solver=solver,
                    C=float(c_value),
                ),
            ),
        ]
    )
    return pipe


def _wt_feature_row(*, frame_feature_csv: Path, mmgbsa_replicate_csv: Path, susceptibility_xlsx: Path) -> pd.DataFrame:
    frame_df = pd.read_csv(frame_feature_csv)
    mmgbsa_df = pd.read_csv(mmgbsa_replicate_csv) if mmgbsa_replicate_csv.exists() else None
    target_df = load_dor_susceptibilities(susceptibility_xlsx)
    wt_row = pd.DataFrame([{"drug": "DOR", "mutation": "WT", "chain": "A", "dor_fold_reduction": 1.0, "order": -1}])
    target_aug = pd.concat([target_df, wt_row], ignore_index=True)
    feat = _mutation_feature_matrix(
        frame_df,
        target_df=target_aug,
        temperature_k=300.0,
        dispersion_mode="replicate_sd",
        mmgbsa_df=mmgbsa_df,
    )
    return feat[feat["mutation"].astype(str) == "WT"].copy()


def _score_candidates(
    *,
    feature_matrix_csv: Path,
    summary_df: pd.DataFrame,
    coef_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    frame_feature_csv: Path,
    mmgbsa_replicate_csv: Path,
    susceptibility_xlsx: Path,
    low_max_fold: float,
) -> pd.DataFrame:
    feat = pd.read_csv(feature_matrix_csv)
    feat_train = feat[feat["control_category"].astype(str).isin({"negative_control", "positive_control"})].copy()
    feat_hold = feat[feat["control_category"].astype(str).isin({"uncertain_limited"})].copy()
    rows: list[dict[str, object]] = []
    for record in candidate_df.itertuples(index=False):
        sub_summary = summary_df[
            (summary_df["penalty"].astype(str) == str(record.penalty))
            & (summary_df["combo_size"].astype(int) == int(record.combo_size))
            & (summary_df["feature_combo"].astype(str) == str(record.feature_combo))
        ].copy()
        penalty = str(sub_summary["penalty"].iloc[0])
        feature_combo = str(sub_summary["feature_combo"].iloc[0])
        features = feature_combo.split("|")
        sub_coef = coef_df[
            (coef_df["penalty"].astype(str) == penalty)
            & (coef_df["combo_size"].astype(int) == int(record.combo_size))
            & (coef_df["feature_combo"].astype(str) == feature_combo)
        ].copy()
        c_value = float(sub_coef["fullfit_best_c"].iloc[0])
        fitted = _logistic_pipeline(penalty=penalty, c_value=c_value)
        y_train = feat_train["control_category"].astype(str).map(
            lambda value: "low" if value == "negative_control" else "high"
        )
        fitted.fit(feat_train[features].copy(), y_train)
        classes = list(fitted.named_steps["model"].classes_)
        pos_idx = int(classes.index("high"))
        hold_rows: list[dict[str, object]] = []
        for _, sample in feat_hold.iterrows():
            sample_x = sample[features].to_frame().T
            prob_high = float(fitted.predict_proba(sample_x)[:, pos_idx][0])
            hold_rows.append(
                {
                    "mutation": str(sample["mutation"]),
                    "control_category": str(sample["control_category"]),
                    "target_fold_reduction": float(sample["target_fold_reduction"]),
                    "prob_high": float(prob_high),
                    "predicted_class": "high" if prob_high >= 0.5 else "low",
                }
            )
        wt_df = _wt_feature_row(
            frame_feature_csv=frame_feature_csv,
            mmgbsa_replicate_csv=mmgbsa_replicate_csv,
            susceptibility_xlsx=susceptibility_xlsx,
        )
        if not wt_df.empty:
            wt_x = wt_df[features].copy()
            wt_prob = float(fitted.predict_proba(wt_x)[:, pos_idx][0])
            hold_rows.append(
                {
                    "mutation": "WT",
                    "control_category": "wt_reference",
                    "target_fold_reduction": 1.0,
                    "prob_high": float(wt_prob),
                    "predicted_class": "high" if wt_prob >= 0.5 else "low",
                }
            )
        holdout = pd.DataFrame(hold_rows)
        wt_row = holdout[holdout["mutation"].astype(str) == "WT"].iloc[0]
        uncertain = holdout[holdout["control_category"].astype(str) == "uncertain_limited"].copy()
        rows.append(
            {
                "penalty": str(record.penalty),
                "combo_size": int(record.combo_size),
                "feature_combo": str(record.feature_combo),
                "wt_prob_high": float(wt_row["prob_high"]),
                "wt_predicted_class": str(wt_row["predicted_class"]),
                "uncertain_high_calls": int((uncertain["predicted_class"].astype(str) == "high").sum()),
                "uncertain_prediction_summary": "; ".join(
                    f"{row.mutation}={row.predicted_class}:{row.prob_high:.3f}" for row in uncertain.itertuples(index=False)
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Rank perfect control-trained logistic models by logit-fold correlation.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression_controls/standard"),
    )
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_controls/tables/feature_matrix_with_categories.csv"),
    )
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--mmgbsa-replicate-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/tables/mmgbsa_replicate_metrics.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--top-n", type=int, default=12)
    args = parser.parse_args()

    out_tables = args.input_dir / "tables"
    summary_df = pd.read_csv(out_tables / "combo_model_summary.csv")
    pred_df = pd.read_csv(out_tables / "combo_cv_predictions.csv")
    coef_df = pd.read_csv(out_tables / "combo_fullfit_coefficients.csv")

    corr_df = _model_correlations(pred_df)
    merged = summary_df.merge(
        corr_df,
        on=["penalty", "combo_size", "feature_combo"],
        how="left",
        validate="one_to_one",
    )
    merged.to_csv(out_tables / "combo_model_summary_with_correlations.csv", index=False)

    perfect = merged[(merged["fp"].astype(int) == 0) & (merged["fn"].astype(int) == 0)].copy()
    perfect = perfect.sort_values(
        ["logit_log10_pearson_r", "logit_raw_pearson_r", "combo_size"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    perfect.to_csv(out_tables / "perfect_models_correlation_scan.csv", index=False)

    top_candidates = pd.concat(
        [
            perfect.sort_values(["combo_size", "logit_log10_pearson_r"], ascending=[True, False]).head(1),
            perfect.sort_values(["logit_raw_pearson_r", "combo_size"], ascending=[False, True]).head(1),
            perfect.sort_values(["logit_log10_pearson_r", "combo_size"], ascending=[False, True]).head(1),
            perfect.head(max(1, int(args.top_n))),
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["penalty", "combo_size", "feature_combo"], keep="first")

    holdout_df = _score_candidates(
        feature_matrix_csv=args.feature_matrix_csv,
        summary_df=summary_df,
        coef_df=coef_df,
        candidate_df=top_candidates[["penalty", "combo_size", "feature_combo"]],
        frame_feature_csv=args.frame_feature_csv,
        mmgbsa_replicate_csv=args.mmgbsa_replicate_csv,
        susceptibility_xlsx=args.susceptibility_xlsx,
        low_max_fold=float(args.low_max_fold),
    )
    candidate_summary = top_candidates.merge(
        holdout_df,
        on=["penalty", "combo_size", "feature_combo"],
        how="left",
        validate="one_to_one",
    )
    candidate_summary.to_csv(out_tables / "top_correlation_candidates.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
