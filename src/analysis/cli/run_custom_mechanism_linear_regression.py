#!/usr/bin/env python3
"""Fit a small linear regression model on the custom mechanism control panel."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet, Lasso, LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_best_combo_regression import _place_greedy_annotations
from .plot_dor_susceptibility_bars import CATEGORY_COLORS


FEATURES = [
    "ser105_dor_distance_angstrom_mean",
    "residue_min_distance_TYR188_angstrom_mean",
    "ligand_pose_rmsd_angstrom_mean",
    "ddg_electrostatic_mean",
]


FEATURE_LABELS = {
    "ser105_dor_distance_angstrom_mean": "SER105-DOR Distance",
    "residue_min_distance_TYR188_angstrom_mean": "TYR188-DOR Distance",
    "ligand_pose_rmsd_angstrom_mean": "Ligand Pose RMSD",
    "ddg_electrostatic_mean": "Electrostatic dDG",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit a custom mechanism linear regression model on DOR controls.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("results/analysis/custom_mechanism_selected_model/tables/mechanism_panel_feature_matrix.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--ddg-summary-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/tables/mutation_ddg_summary.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/custom_mechanism_linear_regression_controls"),
    )
    parser.add_argument("--wt-reference-fold", type=float, default=1.0)
    parser.add_argument("--lasso-alphas", type=str, default="0.001,0.00178,0.00316,0.00562,0.01,0.0178,0.0316,0.0562,0.1,0.133,0.178,0.316,0.562,1,1.78,3.16,5.62,10")
    parser.add_argument("--ridge-alphas", type=str, default="0.001,0.00178,0.00316,0.00562,0.01,0.0178,0.0316,0.0562,0.1,0.133,0.178,0.316,0.562,1,1.78,3.16,5.62,10,17.8,31.6,56.2,100")
    parser.add_argument("--elastic-net-alphas", type=str, default="0.001,0.00178,0.00316,0.00562,0.01,0.0178,0.0316,0.0562,0.1,0.133,0.178,0.316,0.562,1,1.78,3.16,5.62,10")
    parser.add_argument("--elastic-net-l1-ratios", type=str, default="0.2,0.5,0.8")
    return parser.parse_args()


def _parse_float_list(text: str) -> list[float]:
    values: list[float] = []
    for token in str(text).split(","):
        token = token.strip()
        if token:
            values.append(float(token))
    if not values:
        raise ValueError(f"No floats parsed from {text!r}")
    return values


def _normalize_mutation_label(label: str) -> str:
    return str(label).strip().replace(", ", "+").replace(",", "+")


def _display_category(control_category: str, mutation: str) -> str:
    category = str(control_category).strip().lower()
    mutation_label = str(mutation).strip().upper()
    if mutation_label == "WT" or category == "wt_reference":
        return "WT"
    if category == "negative_control":
        return "Negative control"
    if category == "positive_control":
        return "Positive control"
    if category == "uncertain_limited":
        return "Test set"
    return "Other"


def _category_color(display_category: str) -> str:
    if str(display_category) == "WT":
        return "#333333"
    if str(display_category) == "Negative control":
        return CATEGORY_COLORS["Negative control"]
    if str(display_category) == "Positive control":
        return CATEGORY_COLORS["Positive control"]
    if str(display_category) == "Test set":
        return CATEGORY_COLORS["Uncertain/limited data"]
    return "#777777"


def _load_dor_susceptibilities(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path).copy()
    if df.shape[1] < 3:
        raise ValueError(f"Unexpected susceptibility workbook shape: {df.shape}")
    df = df.iloc[:, :3].copy()
    df.columns = ["mutation_sheet", "rpv_fold", "dor_fold"]
    df = df[df["mutation_sheet"].astype(str) != "Mutations"].copy()
    df["mutation"] = df["mutation_sheet"].astype(str).map(_normalize_mutation_label)
    df["dor_fold"] = pd.to_numeric(df["dor_fold"], errors="coerce")
    return df[["mutation", "dor_fold"]].dropna(subset=["dor_fold"]).drop_duplicates("mutation").reset_index(drop=True)


def _load_panel_with_targets(panel_csv: Path, susceptibility_xlsx: Path, ddg_summary_csv: Path, wt_reference_fold: float) -> pd.DataFrame:
    panel = pd.read_csv(panel_csv).copy()
    suscept = _load_dor_susceptibilities(susceptibility_xlsx)
    ddg = pd.read_csv(ddg_summary_csv, usecols=["mutation", "ddg_electrostatic_mean"]).copy()

    merged = panel.drop(columns=["target_fold_reduction"]).merge(suscept, on="mutation", how="left")
    missing = merged[merged["mutation"].astype(str).str.upper() != "WT"]
    missing = missing[missing["dor_fold"].isna()]["mutation"].astype(str).tolist()
    if missing:
        raise ValueError(f"Missing DOR folds for: {missing}")

    merged["target_fold_reduction"] = merged["dor_fold"]
    merged.loc[merged["mutation"].astype(str).str.upper() == "WT", "target_fold_reduction"] = float(wt_reference_fold)
    merged = merged.drop(columns=["dor_fold"]).merge(ddg, on="mutation", how="left")
    merged["display_category"] = [
        _display_category(control_category=row["control_category"], mutation=row["mutation"])
        for _, row in merged.iterrows()
    ]
    return merged


def _make_pipeline(model) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def _evaluate_model_cv(train_df: pd.DataFrame, feature_cols: list[str], model_factory) -> tuple[np.ndarray, dict[str, float]]:
    x = train_df[feature_cols].copy()
    y = np.log10(train_df["target_fold_reduction"].astype(float).to_numpy())
    loo = LeaveOneOut()
    preds = np.zeros(len(train_df), dtype=float)
    for tr_idx, te_idx in loo.split(x, y):
        fitted = _make_pipeline(model_factory())
        fitted.fit(x.iloc[tr_idx], y[tr_idx])
        preds[te_idx[0]] = float(fitted.predict(x.iloc[te_idx])[0])
    metrics = {
        "cv_rmse": float(mean_squared_error(y, preds) ** 0.5),
        "cv_mae": float(mean_absolute_error(y, preds)),
        "cv_r2": float(r2_score(y, preds)),
    }
    return preds, metrics


def _select_model(train_df: pd.DataFrame, args: argparse.Namespace) -> tuple[dict[str, object], np.ndarray]:
    candidates: list[tuple[str, dict[str, float], object]] = [("ols", {}, lambda: LinearRegression())]
    candidates.extend(("ridge", {"alpha": alpha}, lambda alpha=alpha: Ridge(alpha=float(alpha), random_state=0)) for alpha in _parse_float_list(args.ridge_alphas))
    candidates.extend(("lasso", {"alpha": alpha}, lambda alpha=alpha: Lasso(alpha=float(alpha), random_state=0, max_iter=200000)) for alpha in _parse_float_list(args.lasso_alphas))
    for alpha in _parse_float_list(args.elastic_net_alphas):
        for l1_ratio in _parse_float_list(args.elastic_net_l1_ratios):
            candidates.append(
                (
                    "elastic_net",
                    {"alpha": alpha, "l1_ratio": l1_ratio},
                    lambda alpha=alpha, l1_ratio=l1_ratio: ElasticNet(alpha=float(alpha), l1_ratio=float(l1_ratio), random_state=0, max_iter=200000),
                )
            )

    summary_rows: list[dict[str, object]] = []
    best_row: dict[str, object] | None = None
    best_preds: np.ndarray | None = None
    for model_name, params, factory in candidates:
        preds, metrics = _evaluate_model_cv(train_df, FEATURES, factory)
        row = {"model": model_name, "params": json.dumps(params, sort_keys=True), **metrics}
        summary_rows.append(row)
        if best_row is None or float(row["cv_rmse"]) < float(best_row["cv_rmse"]):
            best_row = row
            best_preds = preds

    assert best_row is not None and best_preds is not None
    return {"summary_df": pd.DataFrame(summary_rows).sort_values("cv_rmse", kind="stable").reset_index(drop=True), "best_row": best_row}, best_preds


def _instantiate_best_model(best_row: dict[str, object]):
    model_name = str(best_row["model"])
    params = json.loads(str(best_row["params"]))
    if model_name == "ols":
        return LinearRegression()
    if model_name == "ridge":
        return Ridge(alpha=float(params["alpha"]), random_state=0)
    if model_name == "lasso":
        return Lasso(alpha=float(params["alpha"]), random_state=0, max_iter=200000)
    if model_name == "elastic_net":
        return ElasticNet(alpha=float(params["alpha"]), l1_ratio=float(params["l1_ratio"]), random_state=0, max_iter=200000)
    raise ValueError(f"Unsupported model: {model_name}")


def _plot_predicted_vs_observed(df: pd.DataFrame, output_png: Path, title: str) -> dict[str, float]:
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    x = df["target_fold_reduction"].astype(float).to_numpy()
    y = df["predicted_fold"].astype(float).to_numpy()
    colors = [_category_color(cat) for cat in df["display_category"]]
    ax.scatter(x, y, s=56, c=colors, edgecolors="white", linewidths=0.8, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    lo = float(min(np.min(x), np.min(y)))
    hi = float(max(np.max(x), np.max(y)))
    _place_greedy_annotations(ax, x, y, df["mutation"].astype(str).tolist())
    valid_true = np.log10(df["target_fold_reduction"].astype(float).to_numpy())
    valid_pred = np.log10(df["predicted_fold"].astype(float).to_numpy())
    slope, intercept, r_value, p_value, _stderr = stats.linregress(valid_true, valid_pred)
    x_line = np.geomspace(float(np.min(x)), float(np.max(x)), 200)
    y_line = 10 ** (slope * np.log10(x_line) + intercept)
    ax.plot(x_line, y_line, color="#444444", linewidth=1.4, linestyle="-", zorder=2)
    ax.grid(alpha=0.22, linestyle=":")
    ax.set_xlabel("Observed DOR Fold-change")
    ax.set_ylabel("Predicted DOR Fold-change")
    ax.set_title(f"{title}\n$R^2$ = {r_value**2:.3f}    p = {p_value:.3g}", pad=12)
    present_categories = [cat for cat in ["Negative control", "Positive control", "Test set", "WT"] if cat in set(df["display_category"].astype(str))]
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=_category_color(category),
            markeredgecolor="white",
            markersize=8,
            label=category,
        )
        for category in present_categories
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper left", frameon=True, framealpha=0.92, facecolor="white")
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"r_squared": float(r_value**2), "p_value": float(p_value)}


def _plot_coefficients(coef_df: pd.DataFrame, output_png: Path, title: str) -> None:
    plot_df = coef_df.sort_values("coefficient", kind="stable")
    colors = ["#1d3557" if value < 0 else "#d62828" for value in plot_df["coefficient"]]
    labels = [FEATURE_LABELS.get(feature, feature) for feature in plot_df["feature"]]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.barh(labels, plot_df["coefficient"], color=colors)
    ax.axvline(0.0, color="#444444", linewidth=1.0)
    ax.set_xlabel("Coefficient")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_ranked_predictions(df: pd.DataFrame, output_png: Path, title: str) -> None:
    plot_df = df.sort_values(["predicted_fold", "target_fold_reduction"], ascending=[False, True], kind="stable").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    xs = np.arange(len(plot_df))
    colors = [_category_color(cat) for cat in plot_df["display_category"]]
    ax.bar(xs, plot_df["predicted_fold"], color=colors, edgecolor="white", linewidth=0.7)
    ax.set_yscale("log")
    ax.set_ylabel("Predicted DOR Fold-change")
    ax.set_xticks(xs)
    ax.set_xticklabels(plot_df["mutation"], rotation=45, ha="right", fontsize=8)
    ax.set_title(title)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    panel_df = _load_panel_with_targets(args.input_csv, args.susceptibility_xlsx, args.ddg_summary_csv, args.wt_reference_fold)
    panel_df.to_csv(out_tables / "linear_regression_feature_matrix.csv", index=False)

    train_df = panel_df[panel_df["control_category"].isin(["negative_control", "positive_control"])].copy().reset_index(drop=True)
    selection, cv_preds = _select_model(train_df, args)
    summary_df = selection["summary_df"]
    best_row = selection["best_row"]
    summary_df.to_csv(out_tables / "model_selection_summary.csv", index=False)

    train_cv = train_df[["mutation", "control_category", "display_category", "target_fold_reduction"]].copy()
    train_cv["observed_log10_fold"] = np.log10(train_cv["target_fold_reduction"].astype(float))
    train_cv["predicted_log10_fold"] = cv_preds
    train_cv["predicted_fold"] = np.power(10.0, cv_preds)
    train_cv.to_csv(out_tables / "cv_predictions.csv", index=False)

    estimator = _instantiate_best_model(best_row)
    fitted = _make_pipeline(estimator)
    fitted.fit(train_df[FEATURES], np.log10(train_df["target_fold_reduction"].astype(float).to_numpy()))

    all_df = panel_df.copy()
    pred_log10 = fitted.predict(all_df[FEATURES])
    all_df["predicted_log10_fold"] = pred_log10
    all_df["predicted_fold"] = np.power(10.0, pred_log10)
    all_df["prediction_source"] = np.where(all_df["mutation"].astype(str).str.upper() == "WT", "wt_reference", "full_fit")
    all_df.to_csv(out_tables / "all_mutation_predictions.csv", index=False)

    model = fitted.named_steps["model"]
    coef = np.asarray(model.coef_, dtype=float).reshape(-1)
    coef_df = pd.DataFrame(
        {
            "feature": FEATURES,
            "coefficient": coef,
            "abs_coefficient": np.abs(coef),
        }
    ).sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)
    coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)

    cv_plot_stats = _plot_predicted_vs_observed(
        train_cv,
        out_plots / "cv_predicted_vs_observed_fold_change.png",
        title=f"Control CV Regression: {str(best_row['model']).replace('_', ' ').title()}",
    )
    all_plot_df = all_df.copy().reset_index(drop=True)
    all_plot_stats = _plot_predicted_vs_observed(
        all_plot_df,
        out_plots / "all_mutation_predicted_vs_observed_fold_change.png",
        title="",
    )
    _plot_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png", title="Full-Model Regression Coefficients")
    _plot_ranked_predictions(all_df, out_plots / "all_mutation_predicted_fold_ranked.png", title="All-Mutation Predicted DOR Fold-change")

    stats_rows = [
        {"plot": "cv_predicted_vs_observed_fold_change", **cv_plot_stats},
        {"plot": "all_mutation_predicted_vs_observed_fold_change", **all_plot_stats},
    ]
    pd.DataFrame(stats_rows).to_csv(out_tables / "plot_statistics.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "input_csv": str(args.input_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "ddg_summary_csv": str(args.ddg_summary_csv),
                "output_dir": str(args.output_dir),
                "wt_reference_fold": float(args.wt_reference_fold),
                "features": FEATURES,
                "best_model": best_row,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
