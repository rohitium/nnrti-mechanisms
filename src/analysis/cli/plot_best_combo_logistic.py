#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import confusion_matrix

from .benchmark_resistance_models import _plot_confusion_matrix
from .focus_binary_logistic_classifier import (
    _plot_cv_probability_ranked,
    _plot_probability_vs_fold,
    _plot_selected_coefficients,
)


def _model_slug(rank: int, penalty: str, combo_size: int) -> str:
    return f"rank_{rank:02d}_{str(penalty).lower()}_{int(combo_size)}f"


def _corr_stats(pred_df: pd.DataFrame, *, use_log10_x: bool) -> dict[str, float]:
    df = pred_df.copy()
    df["target_value"] = pd.to_numeric(df["target_value"], errors="coerce").astype(float)
    df["prob_high"] = pd.to_numeric(df["prob_high"], errors="coerce").astype(float)
    if use_log10_x:
        df = df[df["target_value"] > 0.0].copy()
        df["plot_x"] = np.log10(df["target_value"])
        prefix = "log10_fold"
    else:
        df["plot_x"] = df["target_value"]
        prefix = "raw_fold"
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["plot_x", "prob_high"]).reset_index(drop=True)
    if len(df) < 3:
        return {
            f"{prefix}_pearson_r": np.nan,
            f"{prefix}_pearson_pvalue": np.nan,
            f"{prefix}_spearman_rho": np.nan,
            f"{prefix}_spearman_pvalue": np.nan,
            f"{prefix}_r_squared": np.nan,
            f"{prefix}_slope": np.nan,
            f"{prefix}_intercept": np.nan,
            f"{prefix}_linregress_pvalue": np.nan,
        }
    pearson_r, pearson_p = stats.pearsonr(df["plot_x"], df["prob_high"])
    spearman_rho, spearman_p = stats.spearmanr(df["plot_x"], df["prob_high"])
    slope, intercept, r_value, p_value, _stderr = stats.linregress(df["plot_x"], df["prob_high"])
    return {
        f"{prefix}_pearson_r": float(pearson_r),
        f"{prefix}_pearson_pvalue": float(pearson_p),
        f"{prefix}_spearman_rho": float(spearman_rho),
        f"{prefix}_spearman_pvalue": float(spearman_p),
        f"{prefix}_r_squared": float(r_value**2),
        f"{prefix}_slope": float(slope),
        f"{prefix}_intercept": float(intercept),
        f"{prefix}_linregress_pvalue": float(p_value),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot diagnostics for top-performing combo logistic models.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_logistic"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_logistic/diagnostics"),
    )
    parser.add_argument("--top-n", type=int, default=4)
    args = parser.parse_args()

    in_tables = args.input_dir / "tables"
    summary = pd.read_csv(in_tables / "combo_model_summary.csv")
    preds = pd.read_csv(in_tables / "combo_cv_predictions.csv")
    coefs = pd.read_csv(in_tables / "combo_fullfit_coefficients.csv")

    ranked = summary.sort_values(
        ["balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[False, False, False, False, False],
        kind="stable",
    ).head(int(args.top_n)).reset_index(drop=True)

    out_root = args.output_dir
    out_root.mkdir(parents=True, exist_ok=True)
    overall_rows: list[dict[str, object]] = []

    for rank, row in enumerate(ranked.itertuples(index=False), start=1):
        penalty = str(row.penalty)
        combo_size = int(row.combo_size)
        feature_combo = str(row.feature_combo)
        model_dir = out_root / _model_slug(rank, penalty, combo_size)
        out_tables = model_dir / "tables"
        out_plots = model_dir / "plots"
        out_config = model_dir / "config"
        out_tables.mkdir(parents=True, exist_ok=True)
        out_plots.mkdir(parents=True, exist_ok=True)
        out_config.mkdir(parents=True, exist_ok=True)

        pred_df = preds[
            (preds["penalty"].astype(str) == penalty)
            & (preds["combo_size"].astype(int) == combo_size)
            & (preds["feature_combo"].astype(str) == feature_combo)
        ].copy().sort_values(["prob_high", "target_value", "mutation"], ascending=[True, True, True], kind="stable").reset_index(drop=True)
        coef_df = coefs[
            (coefs["penalty"].astype(str) == penalty)
            & (coefs["combo_size"].astype(int) == combo_size)
            & (coefs["feature_combo"].astype(str) == feature_combo)
        ].copy().sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)

        cm = confusion_matrix(pred_df["observed_class"], pred_df["predicted_class"], labels=["low", "high"])
        cm_df = pd.DataFrame(cm, index=["obs_low", "obs_high"], columns=["pred_low", "pred_high"]).reset_index().rename(columns={"index": "observed"})
        corr_raw = _corr_stats(pred_df, use_log10_x=False)
        corr_log = _corr_stats(pred_df, use_log10_x=True)

        title_stub = f"Top Combo Logistic #{rank}: {penalty.upper()} ({combo_size} Features)"
        _plot_cv_probability_ranked(pred_df, out_plots / "cv_probability_ranked.png")
        _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_fold.png", use_log10_x=False)
        _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_log10_fold.png", use_log10_x=True)
        _plot_selected_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png")
        _plot_confusion_matrix(
            cm,
            ["low", "high"],
            f"{title_stub}\nConfusion Matrix",
            out_plots / "confusion_matrix.png",
        )

        model_summary = {
            "rank": int(rank),
            "model_slug": model_dir.name,
            **row._asdict(),
            **corr_raw,
            **corr_log,
        }
        pd.DataFrame([model_summary]).to_csv(out_tables / "model_summary.csv", index=False)
        pred_df.to_csv(out_tables / "cv_predictions.csv", index=False)
        coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)
        cm_df.to_csv(out_tables / "confusion_matrix.csv", index=False)
        (out_config / "run_config.json").write_text(
            json.dumps(
                {
                    "input_dir": str(args.input_dir),
                    "output_dir": str(model_dir),
                    "rank": int(rank),
                    "penalty": penalty,
                    "combo_size": combo_size,
                    "feature_combo": feature_combo.split("|"),
                },
                indent=2,
            )
        )
        overall_rows.append(model_summary)

    pd.DataFrame(overall_rows).to_csv(out_root / "selected_models_summary.csv", index=False)
    (out_root / "config.json").write_text(
        json.dumps(
            {
                "input_dir": str(args.input_dir),
                "output_dir": str(args.output_dir),
                "top_n": int(args.top_n),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
