#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_binding_energy_summary import _place_greedy_annotations
from .plot_dor_susceptibility_bars import CATEGORY_COLORS


FEATURE_SPECS = [
    ("residue_min_distance_SER105_angstrom_mean", "SER105-DOR distance"),
    ("ligand_pose_rmsd_angstrom_mean", "DOR pose RMSD"),
]

DISPLAY_CATEGORIES = {
    "negative_control": "Negative control",
    "positive_control": "Positive control",
    "uncertain_phenotype": "Uncertain phenotype",
    "wt_reference": "WT",
}
PLOT_CATEGORY_ORDER = {"wt_reference": -1, "negative_control": 0, "uncertain_phenotype": 1, "positive_control": 2}

TRAIN_CATEGORIES = {"negative_control", "positive_control", "wt_reference"}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit a two-feature linear regression model for DOR fold-change.")
    parser.add_argument(
        "--replicate-feature-csv",
        type=Path,
        default=Path("results/analysis/new_logistic_regression/fixed_ser105_pose_rmsd/tables/replicate_level_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/new_logistic_regression/fixed_ser105_pose_rmsd_linear_regression"),
    )
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260602)
    return parser.parse_args()


def _pipeline() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LinearRegression()),
        ]
    )


def _feature_cols(rep_df: pd.DataFrame) -> list[str]:
    cols = [col for col, _label in FEATURE_SPECS if col in rep_df.columns]
    if len(cols) != 2:
        raise ValueError(f"Expected two features, found {cols}")
    return cols


def _panel_from_replicates(rep_df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    meta = (
        rep_df[["mutation", "control_category", "target_fold_change"]]
        .drop_duplicates("mutation")
        .reset_index(drop=True)
    )
    panel = rep_df.groupby("mutation", as_index=False)[features].mean().merge(meta, on="mutation", how="left")
    panel["_category_order"] = panel["control_category"].map(PLOT_CATEGORY_ORDER).fillna(99).astype(int)
    panel = panel.sort_values(["_category_order", "target_fold_change", "mutation"], kind="stable")
    return panel.drop(columns="_category_order").reset_index(drop=True)


def _train_df(panel_df: pd.DataFrame) -> pd.DataFrame:
    train = panel_df[panel_df["control_category"].isin(TRAIN_CATEGORIES)].copy().reset_index(drop=True)
    train["target_log10_fold_change"] = np.log10(train["target_fold_change"].astype(float))
    return train


def _fit_predictions(panel_df: pd.DataFrame, train_df: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = _pipeline()
    model.fit(train_df[features], train_df["target_log10_fold_change"].astype(float))
    out = panel_df[["mutation", "control_category", "target_fold_change"]].copy()
    out["observed_log10_fold_change"] = np.log10(out["target_fold_change"].astype(float))
    out["predicted_log10_fold_change"] = model.predict(panel_df[features])
    out["predicted_fold_change"] = np.power(10.0, out["predicted_log10_fold_change"])

    estimator = model.named_steps["model"]
    coef_df = pd.DataFrame(
        {
            "feature": features,
            "label": [dict(FEATURE_SPECS).get(feature, feature) for feature in features],
            "coefficient_per_standardized_feature": np.asarray(estimator.coef_, dtype=float),
        }
    )
    coef_df["abs_coefficient"] = coef_df["coefficient_per_standardized_feature"].abs()
    coef_df = coef_df.sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)
    coef_df.loc[len(coef_df)] = {
        "feature": "intercept",
        "label": "Intercept",
        "coefficient_per_standardized_feature": float(estimator.intercept_),
        "abs_coefficient": abs(float(estimator.intercept_)),
    }
    return out, coef_df


def _grouped_cv(train_rep_df: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, dict[str, float]]:
    logo = LeaveOneGroupOut()
    rows: list[dict[str, object]] = []
    x = train_rep_df[features].copy()
    y = np.log10(train_rep_df["target_fold_change"].astype(float).to_numpy())
    groups = train_rep_df["mutation"].astype(str).to_numpy()
    for train_idx, test_idx in logo.split(x, y, groups=groups):
        heldout = str(groups[test_idx][0])
        model = _pipeline()
        model.fit(x.iloc[train_idx], y[train_idx])
        rep_preds = model.predict(x.iloc[test_idx])
        pred_log10 = float(np.mean(rep_preds))
        observed_fold = float(train_rep_df.iloc[test_idx[0]]["target_fold_change"])
        rows.append(
            {
                "mutation": heldout,
                "control_category": str(train_rep_df.iloc[test_idx[0]]["control_category"]),
                "target_fold_change": observed_fold,
                "observed_log10_fold_change": float(np.log10(observed_fold)),
                "predicted_log10_fold_change": pred_log10,
                "predicted_fold_change": float(10.0**pred_log10),
                "n_replicates": int(len(test_idx)),
                "replicate_predicted_log10_values": "|".join(f"{x:.6g}" for x in rep_preds),
            }
        )
    cv = pd.DataFrame(rows).sort_values(["control_category", "target_fold_change", "mutation"], kind="stable")
    metrics = _regression_metrics(cv["observed_log10_fold_change"], cv["predicted_log10_fold_change"])
    return cv.reset_index(drop=True), {f"grouped_cv_{key}": value for key, value in metrics.items()}


def _regression_metrics(y_true: pd.Series | np.ndarray, y_pred: pd.Series | np.ndarray) -> dict[str, float]:
    y_true_arr = np.asarray(y_true, dtype=float)
    y_pred_arr = np.asarray(y_pred, dtype=float)
    return {
        "n": int(len(y_true_arr)),
        "r2": float(r2_score(y_true_arr, y_pred_arr)) if len(y_true_arr) >= 2 else np.nan,
        "rmse_log10": float(mean_squared_error(y_true_arr, y_pred_arr) ** 0.5),
        "mae_log10": float(mean_absolute_error(y_true_arr, y_pred_arr)),
    }


def _within_genotype_bootstrap(
    rep_df: pd.DataFrame,
    panel_df: pd.DataFrame,
    features: list[str],
    iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    meta = panel_df[["mutation", "control_category", "target_fold_change"]].copy()
    grouped = {mutation: group.reset_index(drop=True) for mutation, group in rep_df.groupby("mutation")}
    pred_rows: list[dict[str, object]] = []
    mutations = sorted(grouped)
    for iteration in range(int(iterations)):
        sampled_rows = []
        for mutation in mutations:
            group = grouped[mutation]
            idx = rng.integers(0, len(group), size=len(group))
            sampled_rows.append(group.iloc[idx])
        boot_rep = pd.concat(sampled_rows, ignore_index=True)
        boot_panel = boot_rep.groupby("mutation", as_index=False)[features].mean().merge(meta, on="mutation", how="left")
        boot_train = _train_df(boot_panel)
        model = _pipeline()
        model.fit(boot_train[features], boot_train["target_log10_fold_change"].astype(float))
        pred_log10 = model.predict(boot_panel[features])
        for mutation, value in zip(boot_panel["mutation"].astype(str), pred_log10):
            pred_rows.append(
                {
                    "bootstrap_iteration": int(iteration),
                    "mutation": mutation,
                    "predicted_log10_fold_change": float(value),
                    "predicted_fold_change": float(10.0**value),
                }
            )
    bootstrap_df = pd.DataFrame(pred_rows)
    quant = (
        bootstrap_df.groupby("mutation")["predicted_log10_fold_change"]
        .quantile([0.025, 0.5, 0.975])
        .unstack()
        .reset_index()
        .rename(
            columns={
                0.025: "predicted_log10_ci_lower",
                0.5: "predicted_log10_bootstrap_median",
                0.975: "predicted_log10_ci_upper",
            }
        )
    )
    for col in ["predicted_log10_ci_lower", "predicted_log10_bootstrap_median", "predicted_log10_ci_upper"]:
        quant[col.replace("log10_", "")] = np.power(10.0, quant[col].astype(float))
    full_fit, _coef = _fit_predictions(panel_df, _train_df(panel_df), features)
    summary = (
        meta.merge(full_fit[["mutation", "predicted_log10_fold_change", "predicted_fold_change"]], on="mutation", how="left")
        .merge(quant, on="mutation", how="left")
        .sort_values(["control_category", "target_fold_change", "mutation"], kind="stable")
        .reset_index(drop=True)
    )
    return bootstrap_df, summary


def _category_color(control_category: str) -> str:
    if control_category == "wt_reference":
        return "#333333"
    display = DISPLAY_CATEGORIES.get(str(control_category), str(control_category))
    return CATEGORY_COLORS.get(display, "#9aa0a6")


def _plot_predicted_vs_fold(summary_df: pd.DataFrame, output_png: Path, *, show_intervals: bool) -> dict[str, float]:
    import matplotlib.pyplot as plt

    df = summary_df[np.isfinite(summary_df["target_fold_change"])].copy()
    fig, ax = plt.subplots(figsize=(9.0, 7.2))
    order = ["wt_reference", "negative_control", "uncertain_phenotype", "positive_control"]
    for category in order:
        subset = df[df["control_category"].astype(str) == category]
        if subset.empty:
            continue
        y = subset["predicted_fold_change"].astype(float).to_numpy()
        kwargs = {
            "fmt": "o",
            "markersize": 7.5,
            "markeredgecolor": "#222222",
            "markeredgewidth": 0.6,
            "color": _category_color(category),
            "alpha": 0.95,
            "zorder": 3,
        }
        if show_intervals:
            lower = subset["predicted_ci_lower"].astype(float).to_numpy()
            upper = subset["predicted_ci_upper"].astype(float).to_numpy()
            kwargs["yerr"] = np.vstack([np.maximum(0, y - lower), np.maximum(0, upper - y)])
            kwargs["capsize"] = 3
            kwargs["elinewidth"] = 1.1
        ax.errorbar(subset["target_fold_change"].astype(float), y, **kwargs)

    fit = stats.linregress(
        np.log10(df["target_fold_change"].astype(float)),
        np.log10(df["predicted_fold_change"].astype(float)),
    )
    x_grid = np.geomspace(float(df["target_fold_change"].min()), float(df["target_fold_change"].max()), 300)
    ax.plot(x_grid, 10.0 ** (fit.slope * np.log10(x_grid) + fit.intercept), color="#222222", linewidth=1.8, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Observed fold-change", fontsize=18, fontweight="bold")
    ax.set_ylabel("Predicted fold-change", fontsize=18, fontweight="bold")
    ax.tick_params(axis="both", labelsize=13)
    ax.grid(alpha=0.24, linestyle=":")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.text(
        0.05,
        0.95,
        f"R$^2$ = {fit.rvalue**2:.2f}\np = {fit.pvalue:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )
    _place_greedy_annotations(
        ax,
        df["target_fold_change"].astype(float).to_numpy(),
        df["predicted_fold_change"].astype(float).to_numpy(),
        df["mutation"].astype(str).tolist(),
        fontsize=10,
        fixed_offsets={
            "WT": (-4, 3),
            "G190A": (-2, 2),
            "K103N+M230L": (-4, 4),
            "A98G+F227C": (-4, 5),
            "Y188L": (34, 5),
            "V106M": (3, 2),
            "V106A": (3, 2),
            "Y318F": (3, -1),
            "L100I+K103N": (8, 1),
            "V106A+F227L": (14, -14),
        },
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=2.0)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return {
        "n": int(len(df)),
        "r_squared_log10": float(fit.rvalue**2),
        "p_value_log10": float(fit.pvalue),
        "slope_log10": float(fit.slope),
        "intercept_log10": float(fit.intercept),
    }


def _plot_prediction_intervals(summary_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    order = summary_df.copy()
    order["_category_order"] = order["control_category"].map(PLOT_CATEGORY_ORDER).fillna(99).astype(int)
    order = order.sort_values(["_category_order", "target_fold_change", "mutation"], kind="stable").reset_index(drop=True)
    x = np.arange(len(order), dtype=float)
    y = order["predicted_fold_change"].astype(float).to_numpy()
    lower = order["predicted_ci_lower"].astype(float).to_numpy()
    upper = order["predicted_ci_upper"].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(max(10.0, len(order) * 0.55), 5.6))
    ax.errorbar(
        x,
        y,
        yerr=np.vstack([np.maximum(0, y - lower), np.maximum(0, upper - y)]),
        fmt="none",
        ecolor="#333333",
        elinewidth=1.1,
        capsize=3,
        zorder=1,
    )
    ax.scatter(
        x,
        y,
        s=56,
        c=[_category_color(cat) for cat in order["control_category"].astype(str)],
        edgecolor="#222222",
        linewidth=0.5,
        zorder=2,
    )
    ax.set_yscale("log")
    ax.set_ylabel("Predicted fold-change", fontsize=15, fontweight="bold")
    ax.set_xticks(x, order["mutation"].astype(str), rotation=50, ha="right", fontsize=10)
    ax.grid(axis="y", alpha=0.22, linestyle=":")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_coefficients(coef_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    plot_df = coef_df[coef_df["feature"] != "intercept"].sort_values(
        "coefficient_per_standardized_feature",
        kind="stable",
    )
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    colors = ["#4c78a8" if x < 0 else "#e45756" for x in plot_df["coefficient_per_standardized_feature"]]
    ax.barh(plot_df["label"], plot_df["coefficient_per_standardized_feature"], color=colors)
    ax.axvline(0.0, color="#444444", linewidth=1.0)
    ax.set_xlabel("Coefficient per standardized feature", fontsize=13, fontweight="bold")
    ax.grid(axis="x", alpha=0.24)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    for directory in (out_tables, out_plots, out_config):
        directory.mkdir(parents=True, exist_ok=True)

    rep_df = pd.read_csv(args.replicate_feature_csv)
    features = _feature_cols(rep_df)
    panel_df = _panel_from_replicates(rep_df, features)
    train = _train_df(panel_df)
    train_rep = rep_df[rep_df["control_category"].isin(TRAIN_CATEGORIES)].copy().reset_index(drop=True)

    full_fit, coef_df = _fit_predictions(panel_df, train, features)
    cv_df, cv_metrics = _grouped_cv(train_rep, features)
    bootstrap_df, bootstrap_summary = _within_genotype_bootstrap(
        rep_df,
        panel_df,
        features,
        int(args.bootstrap_iterations),
        int(args.bootstrap_seed),
    )

    full_metrics = _regression_metrics(full_fit["observed_log10_fold_change"], full_fit["predicted_log10_fold_change"])
    train_full = full_fit[full_fit["control_category"].isin(TRAIN_CATEGORIES)].copy()
    train_metrics = _regression_metrics(train_full["observed_log10_fold_change"], train_full["predicted_log10_fold_change"])

    panel_df.to_csv(out_tables / "linear_regression_feature_matrix.csv", index=False)
    full_fit.to_csv(out_tables / "all_mutation_predictions.csv", index=False)
    cv_df.to_csv(out_tables / "grouped_cv_predictions.csv", index=False)
    bootstrap_df.to_csv(out_tables / "within_genotype_bootstrap_predictions.csv", index=False)
    bootstrap_summary.to_csv(out_tables / "within_genotype_bootstrap_prediction_intervals.csv", index=False)
    coef_df.to_csv(out_tables / "full_model_coefficients.csv", index=False)

    plot_stats = [
        {
            "plot": "within_genotype_bootstrap_predicted_vs_fold_change",
            **_plot_predicted_vs_fold(bootstrap_summary, out_plots / "within_genotype_bootstrap_predicted_vs_fold_change.png", show_intervals=False),
        },
        {
            "plot": "within_genotype_bootstrap_predicted_vs_fold_change_with_intervals",
            **_plot_predicted_vs_fold(bootstrap_summary, out_plots / "within_genotype_bootstrap_predicted_vs_fold_change_with_intervals.png", show_intervals=True),
        },
    ]
    _plot_prediction_intervals(bootstrap_summary, out_plots / "within_genotype_bootstrap_prediction_intervals.png")
    _plot_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png")
    pd.DataFrame(plot_stats).to_csv(out_tables / "plot_statistics.csv", index=False)

    summary = {
        "model": "ordinary_least_squares",
        "target": "log10_fold_change",
        "features": features,
        "train_categories": sorted(TRAIN_CATEGORIES),
        "n_training_genotypes": int(len(train)),
        "n_training_replicates": int(len(train_rep)),
        "n_panel_genotypes": int(len(panel_df)),
        "bootstrap_iterations": int(args.bootstrap_iterations),
        "bootstrap_seed": int(args.bootstrap_seed),
        "training_full_fit": train_metrics,
        "all_panel_full_fit": full_metrics,
        **cv_metrics,
    }
    pd.DataFrame([summary]).to_csv(out_tables / "linear_model_summary.csv", index=False)
    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "replicate_feature_csv": str(args.replicate_feature_csv),
                "output_dir": str(args.output_dir),
                **summary,
            },
            indent=2,
        )
    )
    print(f"Saved {out_tables / 'linear_model_summary.csv'}")
    print(f"Saved {out_plots / 'within_genotype_bootstrap_predicted_vs_fold_change.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
