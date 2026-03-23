#!/usr/bin/env python3
"""Model DOR susceptibility from ligand-pocket structural features.

This script uses ligand-pocket frame features as mutation-level features and
fits simple supervised models against experimental susceptibility from
DRM-susceptibilities.csv.xlsx.

Primary target:
  fold_reduction

Models:
  - Ridge regression
  - Random forest
  - Gradient boosting

Outputs:
  results/analysis/susceptibility_ml/{tables,plots,config}/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from ..susceptibility import load_dor_susceptibilities


R_KJ_PER_MOL_K = 0.00831446261815324
DEFAULT_TEMPERATURE_K = 300.0


def _safe_stat_corr(fn, y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2:
        return float("nan")
    try:
        value, _p = fn(y_true, y_pred)
        return float(value)
    except Exception:
        return float("nan")


def _compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
        "pearson_r": _safe_stat_corr(stats.pearsonr, y_true, y_pred),
        "spearman_rho": _safe_stat_corr(stats.spearmanr, y_true, y_pred),
    }


def _cv_predict_and_importance(
    model_name: str,
    model,
    x: pd.DataFrame,
    y: pd.Series,
    *,
    n_splits: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    splitter = KFold(n_splits=int(n_splits), shuffle=True, random_state=int(random_state))
    preds = np.full(len(y), np.nan, dtype=float)
    fold_ids = np.full(len(y), -1, dtype=int)
    importance_frames: list[pd.DataFrame] = []
    feature_names = list(x.columns)
    for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(x), start=1):
        x_train = x.iloc[train_idx]
        y_train = y.iloc[train_idx]
        x_test = x.iloc[test_idx]
        fitted = model.fit(x_train, y_train)
        preds[test_idx] = fitted.predict(x_test).astype(float)
        fold_ids[test_idx] = int(fold_idx)
        imp = _extract_importance_rows(model_name, fitted, feature_names)
        imp["fold"] = int(fold_idx)
        imp["n_train"] = int(len(train_idx))
        imp["n_test"] = int(len(test_idx))
        importance_frames.append(imp)
    return preds, fold_ids, pd.concat(importance_frames, ignore_index=True)


def _fit_ridge(feature_names: list[str]) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", Ridge(alpha=1.0, random_state=0)),
        ]
    )


def _fit_tree_models() -> dict[str, object]:
    return {
        "random_forest": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    RandomForestRegressor(
                        n_estimators=500,
                        max_depth=4,
                        min_samples_leaf=2,
                        random_state=0,
                    ),
                ),
            ]
        ),
        "gradient_boosting": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    GradientBoostingRegressor(
                        n_estimators=300,
                        learning_rate=0.03,
                        max_depth=2,
                        subsample=0.8,
                        random_state=0,
                    ),
                ),
            ]
        ),
    }


def _extract_importance_rows(model_name: str, fitted_model, feature_names: list[str]) -> pd.DataFrame:
    if model_name == "ridge":
        coef = fitted_model.named_steps["model"].coef_
        raw = np.abs(np.asarray(coef, dtype=float))
    else:
        raw = np.asarray(fitted_model.named_steps["model"].feature_importances_, dtype=float)
    total = float(np.sum(raw))
    norm = raw / total if total > 0 else np.zeros_like(raw)
    out = pd.DataFrame(
        {
            "model": model_name,
            "feature": feature_names,
            "importance_raw": raw.astype(float),
            "importance_norm": norm.astype(float),
        }
    ).sort_values("importance_norm", ascending=False, kind="stable")
    out["rank"] = np.arange(1, len(out) + 1, dtype=int)
    return out.reset_index(drop=True)


def _plot_predictions(df: pd.DataFrame, output_png: Path, target_label: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True, sharey=True)
    models = ["ridge", "random_forest", "gradient_boosting"]
    labels = {
        "ridge": "Ridge",
        "random_forest": "Random Forest",
        "gradient_boosting": "Gradient Boosting",
    }
    colors = {
        "ridge": "#1f77b4",
        "random_forest": "#2a9d8f",
        "gradient_boosting": "#d62828",
    }
    x = df["target_value"].to_numpy(dtype=float)
    lo = float(np.nanmin(x))
    hi = float(np.nanmax(x))
    pad = max(0.25, 0.08 * (hi - lo))
    for ax, model_name in zip(axes, models):
        col = f"pred_{model_name}"
        ax.scatter(
            df["target_value"],
            df[col],
            s=50,
            alpha=0.9,
            color=colors[model_name],
            linewidths=0,
        )
        for _, row in df.iterrows():
            ax.text(
                float(row["target_value"]) + 0.02,
                float(row[col]) + 0.02,
                str(row["mutation"]),
                fontsize=7,
                alpha=0.8,
            )
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#666666", linestyle="--", linewidth=1.0)
        ax.set_title(labels[model_name])
        ax.grid(alpha=0.25)
        ax.set_xlabel(f"Observed {target_label}")
    axes[0].set_ylabel(f"5-Fold CV Predicted {target_label}")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_importance(df: pd.DataFrame, output_png: Path, title: str, top_n: int = 12) -> None:
    top = df.sort_values("importance_norm", ascending=False).head(int(top_n)).copy()
    top = top.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.barh(top["feature"], top["importance_norm"], color="#457b9d")
    ax.set_xlabel("Normalized Importance")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _aggregate_importance(df: pd.DataFrame, *, model_name: str) -> pd.DataFrame:
    agg = (
        df.groupby("feature", as_index=False)
        .agg(
            importance_norm=("importance_norm", "mean"),
            importance_norm_std=("importance_norm", "std"),
            importance_raw=("importance_raw", "mean"),
            importance_raw_std=("importance_raw", "std"),
            mean_rank=("rank", "mean"),
            rank_std=("rank", "std"),
            n_folds=("fold", "nunique"),
        )
        .sort_values(["importance_norm", "mean_rank"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    agg.insert(0, "model", str(model_name))
    return agg


def _compute_feature_associations(x: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    y_values = y.to_numpy(dtype=float)
    for feature in x.columns:
        x_values = x[feature].to_numpy(dtype=float)
        pearson_r, pearson_p = stats.pearsonr(x_values, y_values)
        spearman_rho, spearman_p = stats.spearmanr(x_values, y_values)
        rows.append(
            {
                "feature": feature,
                "pearson_r": float(pearson_r),
                "pearson_pvalue": float(pearson_p),
                "pearson_abs": float(abs(pearson_r)),
                "spearman_rho": float(spearman_rho),
                "spearman_pvalue": float(spearman_p),
                "spearman_abs": float(abs(spearman_rho)),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["spearman_abs", "pearson_abs"], ascending=[False, False], kind="stable")
        .reset_index(drop=True)
    )


def _plot_associations(df: pd.DataFrame, output_png: Path, title: str, top_n: int = 12) -> None:
    top = df.head(int(top_n)).copy().iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    colors = ["#d62828" if value >= 0 else "#1d3557" for value in top["spearman_rho"]]
    ax.barh(top["feature"], top["spearman_rho"], color=colors)
    ax.set_xlabel("Spearman Rho")
    ax.set_title(title)
    ax.axvline(0.0, color="#666666", linewidth=1.0)
    ax.grid(axis="x", alpha=0.25)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_feature_target_scatter_grid(
    feat: pd.DataFrame,
    association_df: pd.DataFrame,
    *,
    target_col: str,
    target_label: str,
    output_png: Path,
    output_dir: Path,
    top_n: int = 12,
) -> None:
    top = association_df.head(int(top_n)).copy()
    if top.empty:
        return

    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_panels = int(len(top))
    n_cols = 3
    n_rows = int(np.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 4.0 * n_rows))
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])
    axes = axes.ravel()

    for ax, (_, row) in zip(axes, top.iterrows()):
        feature = str(row["feature"])
        x = pd.to_numeric(feat[feature], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(feat[target_col], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        labels = feat.loc[mask, "mutation"].astype(str).tolist()

        color = "#d62828" if float(row["spearman_rho"]) >= 0 else "#1d3557"
        ax.scatter(x, y, s=44, color=color, alpha=0.9, linewidths=0)
        grid = None
        slope = None
        intercept = None
        if len(x) >= 2:
            slope, intercept = np.polyfit(x, y, deg=1)
            grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
            ax.plot(grid, slope * grid + intercept, color="#444444", linestyle="--", linewidth=1.1)
        for xi, yi, label in zip(x, y, labels):
            ax.text(float(xi), float(yi), label, fontsize=7, alpha=0.8, ha="left", va="bottom")
        ax.set_title(
            f"{feature}\n"
            f"Spearman={float(row['spearman_rho']):.2f}, Pearson={float(row['pearson_r']):.2f}",
            fontsize=10,
        )
        ax.set_xlabel(feature)
        ax.set_ylabel(target_label)
        ax.grid(alpha=0.22)

        single_fig, single_ax = plt.subplots(figsize=(5.4, 4.4))
        single_ax.scatter(x, y, s=52, color=color, alpha=0.9, linewidths=0)
        if grid is not None and slope is not None and intercept is not None:
            single_ax.plot(grid, slope * grid + intercept, color="#444444", linestyle="--", linewidth=1.1)
        for xi, yi, label in zip(x, y, labels):
            single_ax.text(float(xi), float(yi), label, fontsize=8, alpha=0.8, ha="left", va="bottom")
        single_ax.set_title(
            f"{feature}\n"
            f"Spearman={float(row['spearman_rho']):.2f}, Pearson={float(row['pearson_r']):.2f}",
            fontsize=10,
        )
        single_ax.set_xlabel(feature)
        single_ax.set_ylabel(target_label)
        single_ax.grid(alpha=0.22)
        single_fig.tight_layout()
        single_path = output_dir / f"scatter_{feature}.png"
        single_fig.savefig(single_path, dpi=220, bbox_inches="tight")
        plt.close(single_fig)

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _mutation_feature_matrix(
    frame_df: pd.DataFrame,
    *,
    target_df: pd.DataFrame,
    temperature_k: float,
) -> pd.DataFrame:
    base_feature_cols = [
        c
        for c in frame_df.columns
        if c.startswith("pocket_") or c.startswith("ligand_") or c.startswith("contact_") or c.startswith("residue_")
    ]
    global_feature_cols = [
        c for c in base_feature_cols if not c.startswith("contact_")
    ]

    rep_mean = (
        frame_df.groupby(["mutation", "replicate"], as_index=False)[base_feature_cols]
        .mean()
    )
    mut_mean = (
        rep_mean.groupby("mutation", as_index=False)[base_feature_cols]
        .mean()
        .rename(columns={c: f"{c}_mean" for c in base_feature_cols})
    )
    mut_repstd = (
        rep_mean.groupby("mutation", as_index=False)[global_feature_cols]
        .std(ddof=1)
        .fillna(0.0)
        .rename(columns={c: f"{c}_repstd" for c in global_feature_cols})
    )

    target = target_df.copy()
    target["target_fold_reduction"] = pd.to_numeric(target["dor_fold_reduction"], errors="coerce")
    target["target_log10_fold_reduction"] = np.log10(target["target_fold_reduction"])
    target["target_ddg_exp_kj"] = float(R_KJ_PER_MOL_K * float(temperature_k)) * np.log(target["target_fold_reduction"])

    feat = target.merge(mut_mean, on="mutation", how="inner")
    feat = feat.merge(mut_repstd, on="mutation", how="left")
    feat = feat.fillna(0.0)
    return feat


def main() -> int:
    parser = argparse.ArgumentParser(description="Model DOR susceptibility from ligand-pocket structural features.")
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression/feature_screening"),
    )
    parser.add_argument("--temperature-k", type=float, default=DEFAULT_TEMPERATURE_K)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--target",
        type=str,
        default="target_fold_reduction",
        choices=["target_fold_reduction", "target_ddg_exp_kj", "target_log10_fold_reduction"],
    )
    args = parser.parse_args()

    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)
    if not args.frame_feature_csv.exists():
        raise FileNotFoundError(args.frame_feature_csv)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    target_df = load_dor_susceptibilities(args.susceptibility_xlsx)
    frame_df = pd.read_csv(args.frame_feature_csv)

    feat = _mutation_feature_matrix(
        frame_df,
        target_df=target_df,
        temperature_k=float(args.temperature_k),
    )
    feat = feat.sort_values("target_fold_reduction", ascending=True).reset_index(drop=True)
    feat.to_csv(out_tables / "mutation_feature_matrix.csv", index=False)

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
    x = feat[feature_cols].copy()
    y = feat[str(args.target)].astype(float).copy()
    n_splits = max(2, min(int(args.cv_folds), int(len(y))))

    association_df = _compute_feature_associations(x, y)
    association_df.to_csv(out_tables / "feature_target_associations.csv", index=False)
    _plot_associations(
        association_df,
        out_plots / "feature_target_associations.png",
        title="Top Feature-Target Associations",
    )

    models = {
        "ridge": _fit_ridge(feature_cols),
        **_fit_tree_models(),
    }

    prediction_df = feat[["mutation", "target_fold_reduction", "target_log10_fold_reduction", "target_ddg_exp_kj"]].copy()
    prediction_df = prediction_df.rename(columns={str(args.target): "target_value"}) if str(args.target) in prediction_df.columns else prediction_df
    prediction_df["target_value"] = y.to_numpy(dtype=float)
    cv_fold_ref: np.ndarray | None = None

    summary_rows: list[dict[str, object]] = []
    importance_frames: list[pd.DataFrame] = []
    for model_name, model in models.items():
        preds, fold_ids, fold_imp = _cv_predict_and_importance(
            model_name,
            model,
            x,
            y,
            n_splits=n_splits,
            random_state=int(args.random_state),
        )
        if cv_fold_ref is None:
            cv_fold_ref = fold_ids.copy()
        prediction_df[f"pred_{model_name}"] = preds
        metrics = _compute_metrics(y.to_numpy(dtype=float), preds)
        summary_rows.append(
            {
                "model": model_name,
                "target": str(args.target),
                "n_mutations": int(len(y)),
                "cv_folds": int(n_splits),
                **metrics,
            }
        )
        fold_imp.to_csv(out_tables / f"feature_importance_{model_name}_by_fold.csv", index=False)
        imp = _aggregate_importance(fold_imp, model_name=model_name)
        imp.to_csv(out_tables / f"feature_importance_{model_name}.csv", index=False)
        importance_frames.append(imp)
        _plot_importance(
            imp,
            out_plots / f"feature_importance_{model_name}.png",
            title=f"{model_name.replace('_', ' ').title()} Feature Importance (5-Fold Mean)",
        )

    if cv_fold_ref is not None:
        prediction_df["cv_fold"] = cv_fold_ref.astype(int)
    summary_df = pd.DataFrame(summary_rows).sort_values("mae", ascending=True).reset_index(drop=True)
    summary_df.to_csv(out_tables / "model_cv_summary.csv", index=False)
    prediction_df.to_csv(out_tables / "cv_predictions.csv", index=False)

    consensus = pd.concat(importance_frames, ignore_index=True)
    consensus = (
        consensus.groupby("feature", as_index=False)
        .agg(
            importance_norm=("importance_norm", "mean"),
            importance_norm_std=("importance_norm_std", "mean"),
            mean_rank=("mean_rank", "mean"),
        )
        .sort_values(["importance_norm", "mean_rank"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    consensus.to_csv(out_tables / "feature_importance_consensus.csv", index=False)
    _plot_importance(
        consensus,
        out_plots / "feature_importance_consensus.png",
        title="Consensus Feature Importance (5-Fold Mean)",
    )

    target_label = {
        "target_ddg_exp_kj": "Experimental DDG (kJ/mol)",
        "target_log10_fold_reduction": "log10(Fold Reduction)",
        "target_fold_reduction": "Fold Reduction",
    }[str(args.target)]
    _plot_feature_target_scatter_grid(
        feat,
        association_df,
        target_col=str(args.target),
        target_label=target_label,
        output_png=out_plots / "feature_target_scatter_grid.png",
        output_dir=out_plots / "feature_target_scatter",
    )
    _plot_predictions(prediction_df, out_plots / "cv_predictions.png", target_label=target_label)

    config = {
        "susceptibility_xlsx": str(args.susceptibility_xlsx),
        "frame_feature_csv": str(args.frame_feature_csv),
        "temperature_k": float(args.temperature_k),
        "target": str(args.target),
        "n_mutations": int(len(y)),
        "n_features": int(len(feature_cols)),
        "cv_folds": int(n_splits),
        "random_state": int(args.random_state),
        "feature_columns": feature_cols,
        "models": list(models.keys()),
        "xgboost_available": False,
    }
    (out_config / "run_config.json").write_text(json.dumps(config, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
