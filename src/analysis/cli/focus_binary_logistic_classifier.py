#!/usr/bin/env python3
"""Focused binary logistic resistance analysis with interpretability outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .benchmark_resistance_models import (
    _aggregate_importance,
    _binary_labels,
    _classification_metrics,
    _feature_importance_from_estimator,
    _plot_confusion_matrix,
    _selected_feature_names,
    _selection_k_grid,
)
from .plot_palm_distance_distributions import (
    _feature_slug,
    _fold_map,
    _mutation_summaries,
    _parse_triplets,
    _plot_triplet,
)

PREFERRED_TRIPLET_BASE_FEATURES = [
    "residue_min_distance_LYS101_angstrom",
    "ligand_palm_distance_angstrom",
    "residue_min_distance_PRO236_angstrom",
    "residue_min_distance_VAL108_angstrom",
    "ligand_pose_rmsd_angstrom",
]


def _binary_logistic_pipeline(random_state: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("selector", SelectKBest(score_func=f_classif)),
            (
                "model",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=int(random_state),
                ),
            ),
        ]
    )


def _binary_logistic_grid(n_features: int) -> dict[str, list[object]]:
    return {
        "selector__k": _selection_k_grid(n_features),
        "model__C": [0.1, 1.0, 10.0],
    }


def _base_frame_feature(feature_name: str, frame_feature_columns: set[str]) -> str | None:
    text = str(feature_name).strip()
    for suffix in ("_repstd", "_mean", "_median", "_std"):
        if text.endswith(suffix):
            candidate = text[: -len(suffix)]
            return candidate if candidate in frame_feature_columns else None
    return text if text in frame_feature_columns else None


def _plot_cv_probability_ranked(pred_df: pd.DataFrame, output_png: Path) -> None:
    df = pred_df.copy().sort_values(["prob_high", "target_value"], ascending=[True, True]).reset_index(drop=True)
    x = np.arange(len(df), dtype=float)
    colors = {"low": "#1d3557", "high": "#d62828"}
    markers = {"low": "o", "high": "X"}

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for idx, row in df.iterrows():
        observed = str(row["observed_class"])
        predicted = str(row["predicted_class"])
        face = colors[observed] if observed == predicted else "white"
        edge = colors[observed]
        ax.scatter(
            [idx],
            [float(row["prob_high"])],
            s=88,
            marker=markers[observed],
            facecolors=face,
            edgecolors=edge,
            linewidths=1.6,
            zorder=3,
        )
        ax.text(idx, float(row["prob_high"]) + 0.03, str(row["mutation"]), rotation=60, ha="left", va="bottom", fontsize=7)
    ax.axhline(0.5, color="#666666", linestyle="--", linewidth=1.0)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("CV Predicted Probability Of High")
    ax.set_xlabel("Mutations Sorted By Predicted Probability")
    ax.set_title("Binary Logistic CV Probabilities")
    ax.grid(axis="y", alpha=0.25)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


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
            dx_px = float(dx_pt) * fig.dpi / 72.0
            dy_px = float(dy_pt) * fig.dpi / 72.0
            ann = ax.annotate(
                label,
                xy=(x, y),
                xytext=(dx_pt, dy_pt),
                textcoords="offset points",
                fontsize=8,
                color=text_color,
                alpha=0.88,
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
            alpha=0.88,
            ha="left" if dx_pt >= 0 else "right",
            va="bottom" if dy_pt >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )
        placed_boxes.append((float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)))


def _plot_probability_vs_fold(
    pred_df: pd.DataFrame,
    output_png: Path,
    *,
    use_log10_x: bool,
) -> None:
    df = pred_df.copy()
    df["fold_reduction"] = pd.to_numeric(df["target_value"], errors="coerce")
    if use_log10_x:
        df["plot_x"] = np.log10(df["fold_reduction"])
        x_label = "log10(Fold Reduction)"
        title = "Binary Logistic CV Probability Tracks Resistance Severity"
    else:
        df["plot_x"] = df["fold_reduction"]
        x_label = "Fold Change in DOR Susceptibility"
        title = ""
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["plot_x", "prob_high"]).reset_index(drop=True)
    if df.empty:
        return

    x = df["plot_x"].to_numpy(dtype=float)
    y = df["prob_high"].to_numpy(dtype=float)
    observed = df["observed_class"].astype(str).to_numpy()
    predicted = df["predicted_class"].astype(str).to_numpy()

    pearson_r, pearson_p = stats.pearsonr(x, y)
    slope, intercept, r_value, p_value, _stderr = stats.linregress(x, y)
    x_grid = np.linspace(float(np.min(x)) - 0.05, float(np.max(x)) + 0.05, 200)
    y_grid = slope * x_grid + intercept

    colors = {"low": "#1d3557", "high": "#d62828"}
    markers = {"low": "o", "high": "X"}

    fig, ax = plt.subplots(figsize=(8.6, 5.6) if use_log10_x else (15.5, 7.2))
    for idx, row in df.iterrows():
        obs = str(row["observed_class"])
        pred = str(row["predicted_class"])
        face = colors[obs] if obs == pred else "white"
        edge = colors[obs]
        ax.scatter(
            float(row["plot_x"]),
            float(row["prob_high"]),
            s=96,
            marker=markers[obs],
            facecolors=face,
            edgecolors=edge,
            linewidths=1.7,
            alpha=0.95,
            zorder=3,
        )

    ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.6, zorder=2)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=1.0, alpha=0.9)
    if use_log10_x:
        ax.set_xlim(float(np.min(x_grid)), float(np.max(x_grid)))
    else:
        ax.set_xlim(0.0, 175.0)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Logistic Regression Probability of "High" Resistance')
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.22)
    _place_greedy_annotations(
        ax,
        x,
        y,
        df["mutation"].astype(str).tolist(),
    )
    ax.text(
        0.02,
        0.98,
        f"R^2 = {r_value**2:.3f}\np = {p_value:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.92},
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_selected_coefficients(coef_df: pd.DataFrame, output_png: Path) -> None:
    if coef_df.empty:
        return
    top = coef_df.copy().sort_values("coefficient")
    fig_h = max(4.6, 0.34 * len(top) + 1.8)
    fig, ax = plt.subplots(figsize=(8.8, fig_h))
    colors = ["#1d3557" if x < 0 else "#d62828" for x in top["coefficient"].tolist()]
    ax.barh(top["feature"], top["coefficient"], color=colors)
    ax.axvline(0.0, color="#444444", linewidth=1.0)
    ax.set_xlabel("Standardized Logistic Coefficient")
    ax.set_title("Binary Logistic Full-Model Coefficients")
    ax.grid(axis="x", alpha=0.25)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_feature_contributions(
    contrib_df: pd.DataFrame,
    mutation_summary_df: pd.DataFrame,
    *,
    false_negative_names: set[str],
    output_png: Path,
) -> None:
    if contrib_df.empty:
        return
    order_df = mutation_summary_df.copy().sort_values(["target_value", "fullfit_prob_high", "mutation"], ascending=[True, True, True])
    mutation_order = order_df["mutation"].astype(str).tolist()
    pivot = contrib_df.pivot(index="feature", columns="mutation", values="contribution").fillna(0.0)
    pivot = pivot.reindex(columns=mutation_order)
    fig_w = max(6.5, 0.9 * len(pivot.columns) + 2.2)
    fig_h = max(5.0, 0.35 * len(pivot.index) + 1.8)
    vmax = float(np.nanmax(np.abs(pivot.to_numpy(dtype=float))))
    vmax = max(vmax, 0.1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns.tolist(), rotation=45, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index.tolist())
    for idx, mutation in enumerate(pivot.columns.tolist()):
        if str(mutation) in false_negative_names:
            ax.axvline(idx - 0.5, color="#111111", linewidth=1.0, alpha=0.75)
            ax.axvline(idx + 0.5, color="#111111", linewidth=1.0, alpha=0.75)
            ax.text(idx, -1.2, "FN", ha="center", va="bottom", fontsize=8, color="#d62828", fontweight="bold", clip_on=False)
    for tick in ax.get_xticklabels():
        if tick.get_text() in false_negative_names:
            tick.set_color("#d62828")
            tick.set_fontweight("bold")
    ax.set_title("Full-Model Feature Contributions Across Mutations")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Contribution To High Logit")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _run_cv(
    feat: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    cv_folds: int,
    random_state: int,
    low_max: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = feat[feature_cols].copy()
    y = _binary_labels(feat[target_col].astype(float), low_max=low_max)
    min_class_count = int(y.value_counts().min())
    effective_folds = max(2, min(int(cv_folds), int(min_class_count)))
    outer = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(random_state))

    pos_label = "high"
    pred_rows: list[dict[str, object]] = []
    param_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []
    importance_rows: list[pd.DataFrame] = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(x, y), start=1):
        x_train = x.iloc[train_idx]
        y_train = y.iloc[train_idx]
        x_test = x.iloc[test_idx]
        y_test = y.iloc[test_idx]
        inner_splits = max(2, min(3, int(y_train.value_counts().min())))
        inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=int(random_state) + fold_idx)
        search = GridSearchCV(
            _binary_logistic_pipeline(random_state=int(random_state)),
            _binary_logistic_grid(len(feature_cols)),
            cv=inner.split(x_train, y_train),
            scoring="balanced_accuracy",
            n_jobs=1,
            refit=True,
        )
        search.fit(x_train, y_train)
        fitted = search.best_estimator_
        pred_class = fitted.predict(x_test).astype(str)
        class_order = list(fitted.named_steps["model"].classes_)
        pos_idx = int(class_order.index(pos_label))
        prob_pos = fitted.predict_proba(x_test)[:, pos_idx].astype(float)
        prob_low = 1.0 - prob_pos

        selected = _selected_feature_names(fitted, feature_cols)
        for feature in selected:
            selected_rows.append(
                {
                    "fold": int(fold_idx),
                    "feature": str(feature),
                    "selected": 1,
                }
            )
        imp = _feature_importance_from_estimator(
            "binary_logistic",
            fitted.named_steps["model"],
            feature_cols,
            selected_feature_names=selected,
        )
        if imp is not None:
            imp["fold"] = int(fold_idx)
            importance_rows.append(imp)

        param_rows.append(
            {
                "fold": int(fold_idx),
                "best_params": json.dumps(search.best_params_, sort_keys=True),
                "best_score": float(search.best_score_),
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
            }
        )

        for row_idx, mutation_idx in enumerate(test_idx):
            pred_rows.append(
                {
                    "mutation": str(feat.iloc[mutation_idx]["mutation"]),
                    "target_value": float(feat.iloc[mutation_idx][target_col]),
                    "observed_class": str(y.iloc[mutation_idx]),
                    "predicted_class": str(pred_class[row_idx]),
                    "prob_low": float(prob_low[row_idx]),
                    "prob_high": float(prob_pos[row_idx]),
                    "correct": bool(str(pred_class[row_idx]) == str(y_test.iloc[row_idx])),
                    "fold": int(fold_idx),
                }
            )

    pred_df = pd.DataFrame(pred_rows).sort_values("target_value").reset_index(drop=True)
    param_df = pd.DataFrame(param_rows).sort_values("fold").reset_index(drop=True)
    importance_df = pd.concat(importance_rows, ignore_index=True) if importance_rows else pd.DataFrame()
    return pred_df, param_df, importance_df


def _fit_full_model(
    feat: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    cv_folds: int,
    random_state: int,
    low_max: float,
) -> tuple[GridSearchCV, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = feat[feature_cols].copy()
    y = _binary_labels(feat[target_col].astype(float), low_max=low_max)
    effective_folds = max(2, min(int(cv_folds), int(y.value_counts().min())))
    splitter = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=int(random_state))
    search = GridSearchCV(
        _binary_logistic_pipeline(random_state=int(random_state)),
        _binary_logistic_grid(len(feature_cols)),
        cv=splitter.split(x, y),
        scoring="balanced_accuracy",
        n_jobs=1,
        refit=True,
    )
    search.fit(x, y)
    fitted = search.best_estimator_

    selected = _selected_feature_names(fitted, feature_cols)
    imputer = fitted.named_steps["imputer"]
    scaler = fitted.named_steps["scaler"]
    selector = fitted.named_steps["selector"]
    model = fitted.named_steps["model"]
    x_imputed = imputer.transform(x)
    x_scaled = scaler.transform(x_imputed)
    x_selected = selector.transform(x_scaled)
    coef_raw = np.asarray(model.coef_, dtype=float).reshape(-1)
    intercept_raw = float(np.asarray(model.intercept_, dtype=float).reshape(-1)[0])
    model_classes = list(model.classes_)
    if len(model_classes) != 2:
        raise ValueError(f"Expected binary logistic model, got classes={model_classes}")
    if str(model_classes[1]) == "high":
        coef = coef_raw.astype(float)
        intercept = float(intercept_raw)
    else:
        coef = (-coef_raw).astype(float)
        intercept = float(-intercept_raw)
    logits = np.asarray(x_selected, dtype=float) @ coef + intercept
    prob_pos = 1.0 / (1.0 + np.exp(-logits))
    pred_class = np.where(prob_pos >= 0.5, "high", "low")

    coef_df = pd.DataFrame(
        {
            "feature": selected,
            "coefficient": coef.astype(float),
            "abs_coefficient": np.abs(coef).astype(float),
            "direction": np.where(coef >= 0.0, "toward_high", "toward_low"),
        }
    ).sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)

    score_df = feat[["mutation", target_col]].copy().rename(columns={target_col: "target_value"})
    score_df["observed_class"] = y.astype(str).to_numpy()
    score_df["fullfit_predicted_class"] = pred_class.astype(str)
    score_df["fullfit_logit"] = logits.astype(float)
    score_df["fullfit_prob_high"] = prob_pos.astype(float)
    score_df["fullfit_prob_low"] = 1.0 - prob_pos.astype(float)

    contribution_rows: list[dict[str, object]] = []
    for row_idx, mutation in enumerate(feat["mutation"].astype(str).tolist()):
        for feat_name, coef_value, scaled_value in zip(selected, coef.tolist(), np.asarray(x_selected[row_idx], dtype=float).tolist()):
            contribution_rows.append(
                {
                    "mutation": str(mutation),
                    "feature": str(feat_name),
                    "scaled_value": float(scaled_value),
                    "coefficient": float(coef_value),
                    "contribution": float(scaled_value * coef_value),
                }
            )
    contrib_df = pd.DataFrame(contribution_rows)
    return search, coef_df, score_df, contrib_df


def _write_triplet_feature_plots(
    frame_feature_csv: Path,
    triplets_txt: Path,
    susceptibility_xlsx: Path,
    output_dir: Path,
    features: list[str],
) -> None:
    frame_df = pd.read_csv(frame_feature_csv)
    triplets = _parse_triplets(triplets_txt)
    fold_map = _fold_map(susceptibility_xlsx)

    for feature in features:
        if feature not in frame_df.columns:
            continue
        feature_slug = _feature_slug(feature)
        feature_dir = output_dir / feature_slug
        out_tables = feature_dir / "tables"
        out_plots = feature_dir / "plots"
        out_config = feature_dir / "config"
        out_tables.mkdir(parents=True, exist_ok=True)
        out_plots.mkdir(parents=True, exist_ok=True)
        out_config.mkdir(parents=True, exist_ok=True)

        rep_df, mut_df = _mutation_summaries(frame_df, feature)
        rep_df.to_csv(out_tables / f"{feature_slug}_by_replicate.csv", index=False)
        mut_df.to_csv(out_tables / f"{feature_slug}_by_mutation.csv", index=False)

        triplet_summary_rows: list[dict[str, object]] = []
        for triplet in triplets:
            wanted = set(triplet)
            sub = frame_df[frame_df["mutation"].astype(str).isin(wanted)].copy()
            if sub["mutation"].nunique() != 3:
                continue
            plot_name = f"{feature_slug}_triplet_{str(triplet[0]).replace('+', '_')}_{str(triplet[1]).replace('+', '_')}_{str(triplet[2]).replace('+', '_')}.png"
            result = _plot_triplet(
                sub,
                rep_df,
                triplet,
                feature=feature,
                fold_map=fold_map,
                output_png=out_plots / plot_name,
            )
            triplet_summary_rows.extend(result["summary_rows"])
        pd.DataFrame(triplet_summary_rows).to_csv(out_tables / f"{feature_slug}_triplet_summary.csv", index=False)
        (out_config / "run_config.json").write_text(
            json.dumps(
                {
                    "frame_feature_csv": str(frame_feature_csv),
                    "triplets_txt": str(triplets_txt),
                    "susceptibility_xlsx": str(susceptibility_xlsx),
                    "feature_column": str(feature),
                    "triplets": ["|".join(t) for t in triplets],
                },
                indent=2,
            )
        )


def _choose_triplet_features(
    importance_agg: pd.DataFrame,
    frame_feature_columns: set[str],
    *,
    top_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    by_base = (
        importance_agg.loc[importance_agg["base_frame_feature"].notna(), ["feature", "base_frame_feature", "importance_norm", "rank_mean", "n_folds"]]
        .drop_duplicates(subset=["base_frame_feature"], keep="first")
        .reset_index(drop=True)
    )
    lookup = {str(row["base_frame_feature"]): row.to_dict() for _, row in by_base.iterrows()}
    for feature in PREFERRED_TRIPLET_BASE_FEATURES:
        if feature in frame_feature_columns and feature in lookup and feature not in seen:
            rows.append(dict(lookup[feature]))
            seen.add(feature)
    for _, row in by_base.iterrows():
        base = str(row["base_frame_feature"])
        if base in seen:
            continue
        rows.append(row.to_dict())
        seen.add(base)
        if len(rows) >= int(top_n):
            break
    return pd.DataFrame(rows[: int(top_n)])


def main() -> int:
    parser = argparse.ArgumentParser(description="Focused binary logistic susceptibility analysis.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--triplets-txt",
        type=Path,
        default=Path("results/analysis/triplet_contact_story_100ns/config/triplets.txt"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression"),
    )
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--top-triplet-features", type=int, default=5)
    args = parser.parse_args()

    feat = pd.read_csv(args.feature_matrix_csv)
    target_col = "target_fold_reduction"
    non_feature_cols = {
        "drug",
        "mutation",
        "chain",
        "target_fold_reduction",
    }
    feature_cols = [c for c in feat.columns if c not in non_feature_cols]

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_triplets = args.output_dir / "feature_triplets"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)
    out_triplets.mkdir(parents=True, exist_ok=True)

    pred_df, param_df, importance_df = _run_cv(
        feat,
        feature_cols,
        target_col,
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        low_max=float(args.low_max_fold),
    )
    metrics = _classification_metrics(pred_df["observed_class"].to_numpy(dtype=str), pred_df["predicted_class"].to_numpy(dtype=str))
    cm = confusion_matrix(pred_df["observed_class"], pred_df["predicted_class"], labels=["low", "high"])
    false_neg_df = pred_df[
        (pred_df["observed_class"].astype(str) == "high")
        & (pred_df["predicted_class"].astype(str) == "low")
    ].copy().sort_values("target_value")
    false_neg_df["error_type"] = "false_negative"
    true_low_df = pred_df[
        (pred_df["observed_class"].astype(str) == "low")
        & (pred_df["predicted_class"].astype(str) == "low")
    ].copy().sort_values("target_value")
    class_prob_summary = (
        pred_df.groupby("observed_class", as_index=False)
        .agg(
            n_mutations=("mutation", "count"),
            mean_prob_high=("prob_high", "mean"),
            median_prob_high=("prob_high", "median"),
            std_prob_high=("prob_high", "std"),
        )
    )

    if importance_df.empty:
        importance_agg = pd.DataFrame(columns=["feature", "importance_norm"])
    else:
        importance_agg = _aggregate_importance(importance_df, "binary_logistic")
    frame_feature_columns = set(pd.read_csv(args.frame_feature_csv, nrows=0).columns.tolist())
    importance_agg["base_frame_feature"] = importance_agg["feature"].map(lambda x: _base_frame_feature(str(x), frame_feature_columns))
    top_frame_features = _choose_triplet_features(
        importance_agg,
        frame_feature_columns,
        top_n=int(args.top_triplet_features),
    )

    search, coef_df, score_df, contrib_df = _fit_full_model(
        feat,
        feature_cols,
        target_col,
        cv_folds=int(args.cv_folds),
        random_state=int(args.random_state),
        low_max=float(args.low_max_fold),
    )
    selected_features = coef_df["feature"].tolist()
    class_feature_summary = (
        feat.assign(observed_class=_binary_labels(feat[target_col].astype(float), low_max=float(args.low_max_fold)))
        .groupby("observed_class")[selected_features]
        .agg(["mean", "std", "min", "max"])
    )
    class_feature_summary.columns = ["_".join(col).strip("_") for col in class_feature_summary.columns.to_flat_index()]
    class_feature_summary = class_feature_summary.reset_index()

    false_neg_names = set(false_neg_df["mutation"].astype(str).tolist())

    pred_df.to_csv(out_tables / "cv_predictions.csv", index=False)
    pd.DataFrame([metrics | {"n_mutations": int(len(pred_df)), "cv_folds": int(pred_df["fold"].nunique())}]).to_csv(out_tables / "cv_summary.csv", index=False)
    param_df.to_csv(out_tables / "cv_best_params_by_fold.csv", index=False)
    importance_df.to_csv(out_tables / "feature_importance_by_fold.csv", index=False)
    importance_agg.to_csv(out_tables / "feature_importance.csv", index=False)
    top_frame_features.to_csv(out_tables / "top_triplet_features.csv", index=False)
    false_neg_df.to_csv(out_tables / "false_negative_cases.csv", index=False)
    true_low_df.to_csv(out_tables / "true_low_cases.csv", index=False)
    class_prob_summary.to_csv(out_tables / "class_probability_summary.csv", index=False)
    coef_df.to_csv(out_tables / "full_model_feature_coefficients.csv", index=False)
    score_df.to_csv(out_tables / "full_model_mutation_scores.csv", index=False)
    contrib_df.to_csv(out_tables / "full_model_feature_contributions.csv", index=False)
    class_feature_summary.to_csv(out_tables / "selected_feature_class_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "best_params": json.dumps(search.best_params_, sort_keys=True),
                "best_score": float(search.best_score_),
            }
        ]
    ).to_csv(out_tables / "full_model_best_params.csv", index=False)
    pd.DataFrame(cm, index=["obs_low", "obs_high"], columns=["pred_low", "pred_high"]).reset_index().rename(columns={"index": "observed"}).to_csv(
        out_tables / "confusion_matrix.csv",
        index=False,
    )

    _plot_confusion_matrix(
        cm,
        ["low", "high"],
        "Binary Logistic Confusion Matrix",
        out_plots / "confusion_matrix.png",
    )
    _plot_cv_probability_ranked(pred_df, out_plots / "cv_probability_ranked.png")
    _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_log10_fold.png", use_log10_x=True)
    _plot_probability_vs_fold(pred_df, out_plots / "cv_probability_vs_fold.png", use_log10_x=False)
    _plot_selected_coefficients(coef_df, out_plots / "full_model_feature_coefficients.png")
    _plot_feature_contributions(
        contrib_df,
        score_df,
        false_negative_names=false_neg_names,
        output_png=out_plots / "feature_contributions.png",
    )

    _write_triplet_feature_plots(
        args.frame_feature_csv,
        args.triplets_txt,
        args.susceptibility_xlsx,
        out_triplets,
        top_frame_features["base_frame_feature"].dropna().astype(str).tolist(),
    )

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "frame_feature_csv": str(args.frame_feature_csv),
                "triplets_txt": str(args.triplets_txt),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "target_col": str(target_col),
                "cv_folds_requested": int(args.cv_folds),
                "random_state": int(args.random_state),
                "low_max_fold": float(args.low_max_fold),
                "top_triplet_features": int(args.top_triplet_features),
                "n_mutations": int(len(feat)),
                "n_features": int(len(feature_cols)),
                "feature_columns": feature_cols,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
