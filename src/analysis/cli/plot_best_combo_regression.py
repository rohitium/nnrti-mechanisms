#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score, precision_score, recall_score

from .benchmark_resistance_models import _plot_confusion_matrix


def _place_greedy_annotations(
    ax,
    x_values: np.ndarray,
    y_values: np.ndarray,
    labels: list[str],
    *,
    text_color: str = "#333333",
) -> None:
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    candidate_offsets = [
        (8, 8),
        (8, 18),
        (8, -18),
        (-8, 8),
        (-8, 18),
        (-8, -18),
        (14, 0),
        (-14, 0),
        (14, 22),
        (-14, 22),
        (14, -22),
        (-14, -22),
        (24, 8),
        (-24, 8),
        (24, -8),
        (-24, -8),
    ]
    placed_boxes: list[tuple[float, float, float, float]] = []

    order = np.argsort(x_values)
    for idx in order:
        x = float(x_values[idx])
        y = float(y_values[idx])
        label = str(labels[idx])
        anchor_px = ax.transData.transform((x, y))
        anchor_disp_x = float(anchor_px[0])
        anchor_disp_y = float(anchor_px[1])

        best = None
        best_score = None
        for dx_pt, dy_pt in candidate_offsets:
            ann = ax.annotate(
                label,
                xy=(x, y),
                xytext=(dx_pt, dy_pt),
                textcoords="offset points",
                fontsize=8,
                color=text_color,
                alpha=0.9,
                ha="left" if dx_pt >= 0 else "right",
                va="bottom" if dy_pt >= 0 else "top",
                arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
                bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
            )
            bbox = ann.get_window_extent(renderer=renderer)
            ann.remove()

            overlaps = 0
            overlap_area = 0.0
            for other in placed_boxes:
                x0 = max(float(bbox.x0), other[0])
                y0 = max(float(bbox.y0), other[1])
                x1 = min(float(bbox.x1), other[2])
                y1 = min(float(bbox.y1), other[3])
                if x1 > x0 and y1 > y0:
                    overlaps += 1
                    overlap_area += (x1 - x0) * (y1 - y0)

            dx_px = float(dx_pt) * fig.dpi / 72.0
            dy_px = float(dy_pt) * fig.dpi / 72.0
            dist_penalty = abs(dx_px) + 0.7 * abs(dy_px)
            center_x = 0.5 * (float(bbox.x0) + float(bbox.x1))
            center_y = 0.5 * (float(bbox.y0) + float(bbox.y1))
            anchor_penalty = abs(center_x - anchor_disp_x) + abs(center_y - anchor_disp_y)
            score = (overlaps * 1_000_000.0) + (overlap_area * 100.0) + dist_penalty + 0.15 * anchor_penalty
            if best_score is None or score < best_score:
                best_score = score
                best = (dx_pt, dy_pt, bbox)

        assert best is not None
        dx_pt, dy_pt, bbox = best
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(dx_pt, dy_pt),
            textcoords="offset points",
            fontsize=8,
            color=text_color,
            alpha=0.9,
            ha="left" if dx_pt >= 0 else "right",
            va="bottom" if dy_pt >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )
        placed_boxes.append((float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)))


def _plot_predictions(df: pd.DataFrame, output_png: Path, title: str) -> dict[str, float]:
    x = df["target_value_raw"].to_numpy(dtype=float)
    y = df["predicted_value_model"].to_numpy(dtype=float)
    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_rho, spearman_p = stats.spearmanr(x, y)
    slope, intercept, r_value, p_value, _stderr = stats.linregress(x, y)
    lo = float(min(np.min(x), np.min(y)))
    hi = float(max(np.max(x), np.max(y)))
    pad = max(1.0, 0.06 * (hi - lo))
    grid = np.linspace(lo - pad, hi + pad, 200)

    fig, ax = plt.subplots(figsize=(9.6, 7.4))
    ax.scatter(x, y, s=64, color="#1d3557", alpha=0.92, edgecolors="white", linewidths=0.8, zorder=3)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", color="#888888", linewidth=1.1, zorder=1)
    ax.plot(grid, slope * grid + intercept, linestyle="-", color="#d62828", linewidth=1.5, zorder=2)
    _place_greedy_annotations(ax, x, y, df["mutation"].astype(str).tolist())
    ax.set_xlabel("Observed Fold Reduction")
    ax.set_ylabel("CV Predicted Fold Reduction")
    ax.set_title(title)
    ax.grid(alpha=0.24)
    ax.text(
        0.02,
        0.98,
        (
            f"Pearson r = {pearson_r:.3f}\n"
            f"Spearman rho = {spearman_rho:.3f}\n"
            f"R^2 = {r_value**2:.3f}\n"
            f"p = {p_value:.3f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {
        "pearson_r": float(pearson_r),
        "pearson_pvalue": float(pearson_p),
        "spearman_rho": float(spearman_rho),
        "spearman_pvalue": float(spearman_p),
        "r_squared": float(r_value**2),
        "slope": float(slope),
        "intercept": float(intercept),
    }


def _plot_coefficients(df: pd.DataFrame, output_png: Path, title: str) -> None:
    coef_df = df.copy().sort_values("coefficient")
    colors = ["#1d3557" if value < 0 else "#d62828" for value in coef_df["coefficient"].tolist()]
    fig_h = max(4.6, 0.6 * len(coef_df) + 1.8)
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    ax.barh(coef_df["feature"], coef_df["coefficient"], color=colors)
    ax.axvline(0.0, color="#444444", linewidth=1.0)
    ax.set_xlabel("Coefficient")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _threshold_summary(df: pd.DataFrame, *, low_max_fold: float) -> tuple[pd.DataFrame, pd.Series]:
    obs = np.where(df["target_value_raw"].to_numpy(dtype=float) < float(low_max_fold), "low", "high")
    pred_values = df["predicted_value_model"].to_numpy(dtype=float)
    candidates = sorted(set(pred_values.tolist() + [float(low_max_fold)]))
    rows: list[dict[str, float | str]] = []
    for threshold in candidates:
        pred = np.where(pred_values < float(threshold), "low", "high")
        rows.append(
            {
                "prediction_threshold": float(threshold),
                "accuracy": float(accuracy_score(obs, pred)),
                "balanced_accuracy": float(balanced_accuracy_score(obs, pred)),
                "macro_f1": float(f1_score(obs, pred, average="macro", zero_division=0)),
                "macro_precision": float(precision_score(obs, pred, average="macro", zero_division=0)),
                "macro_recall": float(recall_score(obs, pred, average="macro", zero_division=0)),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["balanced_accuracy", "macro_f1", "accuracy", "prediction_threshold"],
        ascending=[False, False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    best = summary.iloc[0].copy()
    best["observed_low_max_fold"] = float(low_max_fold)
    return summary, best


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot diagnostics for the best-performing combo regression model.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_regression"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/feature_combo_regression/best_model"),
    )
    parser.add_argument("--rank-by", type=str, default="r2", choices=["r2", "pearson_r", "spearman_rho"])
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    args = parser.parse_args()

    in_tables = args.input_dir / "tables"
    summary = pd.read_csv(in_tables / "combo_model_summary.csv")
    preds = pd.read_csv(in_tables / "combo_cv_predictions.csv")
    coefs = pd.read_csv(in_tables / "combo_fullfit_coefficients.csv")

    best_row = summary.sort_values([str(args.rank_by), "rmse"], ascending=[False, True], kind="stable").iloc[0].copy()
    model_name = str(best_row["model"])
    combo_size = int(best_row["combo_size"])
    feature_combo = str(best_row["feature_combo"])

    pred_df = preds[
        (preds["model"].astype(str) == model_name)
        & (preds["combo_size"].astype(int) == combo_size)
        & (preds["feature_combo"].astype(str) == feature_combo)
    ].copy().sort_values("target_value_raw").reset_index(drop=True)
    coef_df = coefs[
        (coefs["model"].astype(str) == model_name)
        & (coefs["combo_size"].astype(int) == combo_size)
        & (coefs["feature_combo"].astype(str) == feature_combo)
    ].copy().sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    corr_stats = _plot_predictions(
        pred_df,
        out_plots / "cv_predicted_vs_observed_fold_change.png",
        title=f"Best Combo Regression: {model_name.title()} ({combo_size} Features)",
    )
    _plot_coefficients(
        coef_df,
        out_plots / "full_model_feature_coefficients.png",
        title=f"Best Combo Full-Fit Coefficients: {model_name.title()}",
    )

    threshold_df, best_threshold = _threshold_summary(pred_df, low_max_fold=float(args.low_max_fold))
    threshold = float(best_threshold["prediction_threshold"])
    observed = np.where(pred_df["target_value_raw"].to_numpy(dtype=float) < float(args.low_max_fold), "low", "high")
    predicted = np.where(pred_df["predicted_value_model"].to_numpy(dtype=float) < threshold, "low", "high")
    cm = confusion_matrix(observed, predicted, labels=["low", "high"])
    _plot_confusion_matrix(
        cm,
        ["low", "high"],
        f"Best Regression Confusion Matrix\nPrediction Threshold = {threshold:.2f}",
        out_plots / "confusion_matrix_best_threshold.png",
    )

    cm_df = pd.DataFrame(cm, index=["obs_low", "obs_high"], columns=["pred_low", "pred_high"]).reset_index().rename(columns={"index": "observed"})
    cm_df.to_csv(out_tables / "confusion_matrix_best_threshold.csv", index=False)
    threshold_df.to_csv(out_tables / "threshold_sweep.csv", index=False)
    pd.DataFrame([best_row.to_dict() | corr_stats | best_threshold.to_dict()]).to_csv(out_tables / "best_model_summary.csv", index=False)
    pred_df.assign(observed_class=observed, predicted_class=predicted).to_csv(out_tables / "best_model_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "best_model_coefficients.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "input_dir": str(args.input_dir),
                "output_dir": str(args.output_dir),
                "rank_by": str(args.rank_by),
                "low_max_fold": float(args.low_max_fold),
                "selected_model": model_name,
                "selected_combo_size": combo_size,
                "selected_feature_combo": feature_combo.split("|"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
