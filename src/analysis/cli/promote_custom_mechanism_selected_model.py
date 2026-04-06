#!/usr/bin/env python3
"""Promote the selected custom mechanism logistic model into a clean result folder with plots."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score, roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_dor_susceptibility_bars import CATEGORY_COLORS, NEGATIVE_CONTROLS, POSITIVE_CONTROLS, UNCERTAIN_LIMITED


SELECTED_FEATURES = [
    "ser105_dor_distance_angstrom_mean",
    "residue_min_distance_TYR188_angstrom_mean",
    "ligand_pose_rmsd_angstrom_mean",
]


FEATURE_LABELS = {
    "ser105_dor_distance_angstrom_mean": "SER105-DOR Distance",
    "residue_min_distance_TYR188_angstrom_mean": "TYR188-DOR Distance",
    "ligand_pose_rmsd_angstrom_mean": "Ligand Pose RMSD",
}

DISPLAY_TEST_SET_LABEL = "Test set"
DISPLAY_CATEGORY_COLORS = {
    "Negative control": CATEGORY_COLORS["Negative control"],
    "Positive control": CATEGORY_COLORS["Positive control"],
    DISPLAY_TEST_SET_LABEL: CATEGORY_COLORS["Uncertain/limited data"],
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Promote the selected custom mechanism logistic model.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("results/analysis/custom_mechanism_selected_model/tables/mechanism_panel_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/custom_mechanism_selected_model"),
    )
    parser.add_argument("--penalty", type=str, default="l1")
    parser.add_argument("--c-value", type=float, default=3.0)
    parser.add_argument("--decision-threshold", type=float, default=0.5)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    return parser.parse_args()


def _category_color(label: str) -> str:
    category = _display_category(label)
    if category == "WT":
        return "#333333"
    if category in DISPLAY_CATEGORY_COLORS:
        return DISPLAY_CATEGORY_COLORS[category]
    return "#777777"


def _display_category(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation == "WT":
        return "WT"
    if mutation in {"F227C", "V106I"}:
        return "Negative control"
    if mutation in UNCERTAIN_LIMITED:
        return DISPLAY_TEST_SET_LABEL
    if mutation in NEGATIVE_CONTROLS:
        return "Negative control"
    if mutation in POSITIVE_CONTROLS:
        return "Positive control"
    return "Other"


def _filter_selected_model_plot_df(df: pd.DataFrame) -> pd.DataFrame:
    plot_df = df.copy()
    if "display_category" not in plot_df.columns:
        plot_df["display_category"] = plot_df["mutation"].astype(str).map(_display_category)
    return plot_df.reset_index(drop=True)


def _display_control_cv_df(df: pd.DataFrame) -> pd.DataFrame:
    return _filter_selected_model_plot_df(df)


def _display_holdout_rank_df(df: pd.DataFrame) -> pd.DataFrame:
    plot_df = _filter_selected_model_plot_df(df)
    return plot_df[plot_df["display_category"].isin([DISPLAY_TEST_SET_LABEL, "WT"])].reset_index(drop=True)


def _pipeline(penalty: str, c_value: float) -> Pipeline:
    solver = "liblinear" if penalty == "l1" else "lbfgs"
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(penalty=penalty, C=float(c_value), solver=solver, max_iter=5000, random_state=0)),
        ]
    )


def _safe_logit(x: np.ndarray) -> np.ndarray:
    return logit(np.clip(np.asarray(x, dtype=float), 1e-6, 1.0 - 1e-6))


def _annotate_points(ax: plt.Axes, df: pd.DataFrame, overrides: dict[str, tuple[int, int]] | None = None) -> None:
    overrides = overrides or {}
    right_side_candidates = [(-10, 5), (-10, -8), (-14, 10), (-14, -12), (-18, 4), (-18, -10)]
    left_side_candidates = [(5, 5), (5, -8), (9, 10), (9, -12), (14, 4), (14, -10)]
    middle_candidates = [(5, 5), (5, -8), (-10, 5), (-10, -8), (9, 10), (9, -12)]
    renderer = ax.figure.canvas.get_renderer()
    placed: list[object] = []
    xvals = df["target_fold_reduction"].astype(float).to_numpy()
    xmin = float(np.min(xvals))
    xmax = float(np.max(xvals))
    lxmin = np.log10(xmin)
    lxmax = np.log10(xmax)
    for _, row in df.sort_values(["target_fold_reduction", "prob_high"], kind="stable").iterrows():
        label = str(row["mutation"])
        xy = (float(row["target_fold_reduction"]), float(row["prob_high"]))
        color = _category_color(label)
        if label in overrides:
            candidates = [overrides[label]]
        else:
            xfrac = (np.log10(xy[0]) - lxmin) / max(lxmax - lxmin, 1e-8)
            if xfrac > 0.72:
                candidates = right_side_candidates
            elif xfrac < 0.22:
                candidates = left_side_candidates
            else:
                candidates = middle_candidates
        best_ann = None
        best_score = None
        for dx, dy in candidates:
            ann = ax.annotate(
                label,
                xy,
                textcoords="offset points",
                xytext=(dx, dy),
                fontsize=8,
                ha="right" if dx < 0 else "left",
                color=color,
            )
            ann.set_path_effects([pe.withStroke(linewidth=3, foreground="white", alpha=0.92)])
            bbox = ann.get_window_extent(renderer=renderer).expanded(1.05, 1.15)
            overlap = sum(bbox.overlaps(prev) for prev in placed)
            score = (overlap, abs(dy), abs(dx))
            if best_score is None or score < best_score:
                if best_ann is not None:
                    best_ann.remove()
                best_ann = ann
                best_score = score
            else:
                ann.remove()
        if best_ann is not None:
            placed.append(best_ann.get_window_extent(renderer=renderer).expanded(1.05, 1.15))


def _plot_cv_vs_fold(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for _, row in df.iterrows():
        m = str(row["mutation"])
        x = float(row["target_fold_reduction"])
        y = float(row["prob_high"])
        c = _category_color(m)
        ax.scatter(x, y, s=48, color=c, edgecolor="white", linewidth=0.7, zorder=3)
    x = df["target_fold_reduction"].astype(float).to_numpy()
    y = df["prob_high"].astype(float).to_numpy()
    logx = np.log10(x)
    model = LinearRegression().fit(logx.reshape(-1, 1), y)
    x_line = np.geomspace(float(np.min(x)), float(np.max(x)), 200)
    y_line = model.predict(np.log10(x_line).reshape(-1, 1))
    y_line = np.clip(y_line, 0.0, 1.0)
    r, p_value = pearsonr(logx, y)
    ax.plot(x_line, y_line, color="#444444", linewidth=1.6, linestyle="--", zorder=2)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("DOR Fold-change")
    ax.set_ylabel("LOO CV P(high)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.22, linestyle=":")
    fig.canvas.draw()
    _annotate_points(ax, df)
    ax.set_title(f"$R^2$ = {r**2:.3f}    p = {p_value:.3g}", fontsize=10, pad=10)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_all_probability_vs_fold(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    marker_map = {"cv_control": "o", "holdout": "s", "wt_reference": "D"}
    for _, row in df.iterrows():
        m = str(row["mutation"])
        x = float(row["target_fold_reduction"])
        y = float(row["prob_high"])
        c = _category_color(m)
        marker = marker_map.get(str(row["prediction_source"]), "o")
        ax.scatter(x, y, s=52, marker=marker, color=c, edgecolor="white", linewidth=0.7, zorder=3)
    x = df["target_fold_reduction"].astype(float).to_numpy()
    y = df["prob_high"].astype(float).to_numpy()
    logx = np.log10(x)
    model = LinearRegression().fit(logx.reshape(-1, 1), y)
    x_line = np.geomspace(float(np.min(x)), float(np.max(x)), 200)
    y_line = model.predict(np.log10(x_line).reshape(-1, 1))
    y_line = np.clip(y_line, 0.0, 1.0)
    r, p_value = pearsonr(logx, y)
    ax.plot(x_line, y_line, color="#444444", linewidth=1.6, linestyle="--", zorder=2)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("DOR Fold-change")
    ax.set_ylabel("P(high)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.22, linestyle=":")
    fig.canvas.draw()
    _annotate_points(
        ax,
        df,
        overrides={
            "Y188L": (-8, 8),
            "V106A+P225H": (-12, -8),
            "V106A+F227L": (8, 6),
            "V106A+L234I": (8, 6),
            "V106I+F227C": (8, -8),
            "A98G+F227C": (8, 6),
        },
    )
    category_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Negative control"], markeredgecolor="white", markersize=8, label="Negative control"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Positive control"], markeredgecolor="white", markersize=8, label="Positive control"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Uncertain/limited data"], markeredgecolor="white", markersize=8, label=DISPLAY_TEST_SET_LABEL),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#333333", markeredgecolor="white", markersize=8, label="WT"),
    ]
    ax.legend(handles=category_handles, loc="center left", frameon=True, facecolor="white", framealpha=0.92)
    ax.set_title(f"$R^2$ = {r**2:.3f}    p = {p_value:.3g}", fontsize=10, pad=10)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_ranked_probabilities(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.sort_values(["prob_high", "target_fold_reduction"], ascending=[False, True], kind="stable").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    xs = np.arange(len(plot_df))
    colors = [_category_color(m) for m in plot_df["mutation"]]
    ax.bar(xs, plot_df["prob_high"], color=colors, edgecolor="white", linewidth=0.7)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.set_ylabel("P(high)")
    ax.set_ylim(0.0, 1.02)
    ax.set_xticks(xs)
    ax.set_xticklabels(plot_df["mutation"], rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_coefficients(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.sort_values("abs_coefficient", ascending=True, kind="stable")
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    colors = ["#c23b22" if v < 0 else "#2a6fbb" for v in plot_df["coefficient"]]
    labels = [FEATURE_LABELS.get(str(f), str(f)) for f in plot_df["feature"]]
    ax.barh(labels, plot_df["coefficient"], color=colors)
    ax.axvline(0.0, color="#444444", linewidth=1.0)
    ax.set_xlabel("Coefficient")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_confusion_matrix(cm: np.ndarray, out: Path) -> None:
    display_cm = cm[[1, 0], :]
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    im = ax.imshow(display_cm, cmap="Blues")
    for i in range(display_cm.shape[0]):
        for j in range(display_cm.shape[1]):
            ax.text(j, i, str(int(display_cm[i, j])), ha="center", va="center", color="#111111", fontsize=12)
    ax.set_xticks([0, 1], labels=["Pred low", "Pred high"])
    ax.set_yticks([0, 1], labels=["True high", "True low"])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_holdout(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.sort_values(["prob_high", "target_fold_reduction"], ascending=[False, True], kind="stable").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    xs = np.arange(len(plot_df))
    colors = [_category_color(m) for m in plot_df["mutation"]]
    ax.bar(xs, plot_df["prob_high"], color=colors, edgecolor="white", linewidth=0.7)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.set_ylabel("Full-fit P(high)")
    ax.set_ylim(0.0, 1.02)
    ax.set_xticks(xs)
    ax.set_xticklabels(plot_df["mutation"], rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def _plot_roc_curve(y_true: pd.Series, probs: np.ndarray, out: Path) -> None:
    y_bin = (y_true.astype(str) == "high").astype(int).to_numpy(dtype=int)
    fpr, tpr, _ = roc_curve(y_bin, probs)
    auc = roc_auc_score(y_bin, probs)

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.plot(fpr, tpr, color="#2a6fbb", linewidth=2.0, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], color="#888888", linestyle="--", linewidth=1.0)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Control LOO ROC Curve", fontsize=11, pad=10)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, facecolor="white", edgecolor="#cccccc")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    if not args.input_csv.exists():
        raise FileNotFoundError(args.input_csv)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv).copy()
    df.to_csv(out_tables / "mechanism_panel_feature_matrix.csv", index=False)
    train = df[df["control_category"].isin(["negative_control", "positive_control"])].copy().reset_index(drop=True)
    hold = df[df["control_category"].isin(["uncertain_limited", "wt_reference"])].copy().reset_index(drop=True)
    x_train = train[SELECTED_FEATURES].copy()
    y_train = train["target_binary_class"].astype(str).copy()

    loo = LeaveOneOut()
    cv_rows: list[dict[str, object]] = []
    probs = np.zeros(len(train), dtype=float)
    for tr_idx, te_idx in loo.split(x_train, y_train):
        fitted = _pipeline(str(args.penalty), float(args.c_value))
        fitted.fit(x_train.iloc[tr_idx], y_train.iloc[tr_idx])
        pos_idx = list(fitted.named_steps["model"].classes_).index("high")
        prob = float(fitted.predict_proba(x_train.iloc[te_idx])[:, pos_idx][0])
        probs[te_idx[0]] = prob
        cv_rows.append(
            {
                "mutation": str(train.loc[te_idx[0], "mutation"]),
                "control_category": str(train.loc[te_idx[0], "control_category"]),
                "target_fold_reduction": float(train.loc[te_idx[0], "target_fold_reduction"]),
                "observed_class": str(train.loc[te_idx[0], "target_binary_class"]),
                "prob_high": prob,
                "prob_low": float(1.0 - prob),
                "predicted_class": "high" if prob >= float(args.decision_threshold) else "low",
            }
        )
    cv_df = pd.DataFrame(cv_rows).sort_values("target_fold_reduction", kind="stable").reset_index(drop=True)
    cv_df["logit_prob_high"] = _safe_logit(cv_df["prob_high"].to_numpy(dtype=float))
    cv_df["prediction_source"] = "cv_control"
    cv_df.to_csv(out_tables / "cv_predictions.csv", index=False)

    y_pred = np.where(probs >= float(args.decision_threshold), "high", "low")
    cm = confusion_matrix(y_train, y_pred, labels=["low", "high"])
    tn, fp, fn, tp = [int(v) for v in cm.ravel()]
    summary = pd.DataFrame(
        [
            {
                "penalty": str(args.penalty),
                "c_value": float(args.c_value),
                "decision_threshold": float(args.decision_threshold),
                "n_features": int(len(SELECTED_FEATURES)),
                "features": "|".join(SELECTED_FEATURES),
                "cv_scheme": "leave_one_out",
                "accuracy": float(accuracy_score(y_train, y_pred)),
                "balanced_accuracy": float(balanced_accuracy_score(y_train, y_pred)),
                "macro_f1": float(f1_score(y_train, y_pred, average="macro")),
                "roc_auc": float(roc_auc_score((y_train == "high").astype(int), probs)),
                "pearson_r_logit_vs_fold": float(pearsonr(_safe_logit(probs), train["target_fold_reduction"].to_numpy(dtype=float)).statistic),
                "pearson_r_logit_vs_log10_fold": float(
                    pearsonr(_safe_logit(probs), np.log10(train["target_fold_reduction"].to_numpy(dtype=float))).statistic
                ),
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "tp": tp,
            }
        ]
    )
    summary.to_csv(out_tables / "model_summary.csv", index=False)

    full = _pipeline(str(args.penalty), float(args.c_value))
    full.fit(x_train, y_train)
    coef = full.named_steps["model"].coef_[0]
    coef_df = pd.DataFrame({"feature": SELECTED_FEATURES, "coefficient": coef, "abs_coefficient": np.abs(coef)}).sort_values(
        "abs_coefficient", ascending=False, kind="stable"
    )
    coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)

    pos_idx = list(full.named_steps["model"].classes_).index("high")
    hold_out = hold[["mutation", "control_category", "target_fold_reduction", "target_binary_class"]].copy()
    hold_out["prob_high"] = full.predict_proba(hold[SELECTED_FEATURES])[:, pos_idx]
    hold_out["prob_low"] = 1.0 - hold_out["prob_high"]
    hold_out["logit_prob_high"] = _safe_logit(hold_out["prob_high"].to_numpy(dtype=float))
    hold_out["predicted_class"] = np.where(hold_out["prob_high"] >= float(args.decision_threshold), "high", "low")
    hold_out["prediction_source"] = np.where(hold_out["control_category"] == "wt_reference", "wt_reference", "holdout")
    hold_out.to_csv(out_tables / "holdout_predictions.csv", index=False)

    all_prob = pd.concat(
        [
            cv_df[
                [
                    "mutation",
                    "control_category",
                    "target_fold_reduction",
                    "observed_class",
                    "prob_high",
                    "prob_low",
                    "logit_prob_high",
                    "predicted_class",
                    "prediction_source",
                ]
            ].rename(columns={"observed_class": "target_binary_class"}),
            hold_out[
                [
                    "mutation",
                    "control_category",
                    "target_fold_reduction",
                    "target_binary_class",
                    "prob_high",
                    "prob_low",
                    "logit_prob_high",
                    "predicted_class",
                    "prediction_source",
                ]
            ],
        ],
        ignore_index=True,
    ).sort_values("target_fold_reduction", kind="stable").reset_index(drop=True)
    all_prob.to_csv(out_tables / "all_mutation_probabilities.csv", index=False)

    cv_plot_df = _display_control_cv_df(cv_df)
    all_prob_plot_df = _filter_selected_model_plot_df(all_prob)
    hold_out_plot_df = _display_holdout_rank_df(all_prob)

    display_probs = cv_plot_df["prob_high"].to_numpy(dtype=float)
    display_y_true = cv_plot_df["observed_class"].astype(str)
    display_y_pred = np.where(display_probs >= float(args.decision_threshold), "high", "low")
    display_cm = confusion_matrix(display_y_true, display_y_pred, labels=["low", "high"])

    _plot_cv_vs_fold(cv_plot_df, out_plots / "cv_probability_vs_fold_change.png")
    _plot_all_probability_vs_fold(all_prob_plot_df, out_plots / "all_mutation_probability_vs_fold_change.png")
    _plot_ranked_probabilities(cv_plot_df, out_plots / "cv_probability_ranked.png")
    _plot_roc_curve(display_y_true, display_probs, out_plots / "cv_roc_curve.png")
    _plot_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png")
    _plot_confusion_matrix(display_cm, out_plots / "confusion_matrix.png")
    _plot_holdout(hold_out_plot_df, out_plots / "holdout_probability_ranked.png")

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "input_csv": str(args.input_csv),
                "output_dir": str(args.output_dir),
                "penalty": str(args.penalty),
                "c_value": float(args.c_value),
                "decision_threshold": float(args.decision_threshold),
                "low_max_fold": float(args.low_max_fold),
                "selected_features": SELECTED_FEATURES,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
