#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from ..susceptibility import load_dor_susceptibilities
from .model_susceptibility_from_state_features import _mutation_feature_matrix


META_COLUMNS = {"drug", "mutation", "chain", "target_fold_reduction", "target_binary_class"}
FEATURE_COLORS = ["#1d3557", "#d62828", "#2a9d8f", "#f4a261", "#6d597a"]


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
        (30, 18),
        (-30, 18),
        (30, -18),
        (-30, -18),
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
                fontsize=9,
                color=text_color,
                alpha=0.92,
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
            fontsize=9,
            color=text_color,
            alpha=0.92,
            ha="left" if dx_pt >= 0 else "right",
            va="bottom" if dy_pt >= 0 else "top",
            arrowprops={"arrowstyle": "-", "color": "#999999", "linewidth": 0.8, "alpha": 0.7},
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )
        placed_boxes.append((float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)))


def _feature_slug(feature_name: str) -> str:
    return str(feature_name).replace("_angstrom", "").replace("_", "_")


def _pretty_feature_label(feature_name: str) -> str:
    text = str(feature_name)
    stat_suffix = ""
    if text.endswith("_repstd"):
        text = text[: -len("_repstd")]
        stat_suffix = " (Replicate SD)"
    elif text.endswith("_sd"):
        text = text[: -len("_sd")]
        stat_suffix = " (SD)"
    elif text.endswith("_mean"):
        text = text[: -len("_mean")]
        stat_suffix = " (Mean)"
    elif text.endswith("_std"):
        text = text[: -len("_std")]
        stat_suffix = " (SD)"

    residue_prefix = "residue_min_distance_"
    if text.startswith(residue_prefix):
        residue = text[len(residue_prefix) :]
        residue = residue.replace("_angstrom", "")
        label = f"{residue[:3].title()}{residue[3:]} Minimum Distance"
        return f"{label}{stat_suffix}"

    replacements = {
        "ligand_palm_distance_angstrom": "Ligand Palm Distance",
        "ligand_pose_rmsd_angstrom": "Ligand Pose RMSD",
        "ligand_entrance_distance_angstrom": "Ligand Entrance Distance",
        "ligand_pocket_center_distance_angstrom": "Ligand Pocket Center Distance",
        "ligand_palm_depth_projection_angstrom": "Ligand Palm Depth Projection",
        "pocket_ca_rmsd_angstrom": "Pocket CA RMSD",
    }
    if text in replacements:
        return f"{replacements[text]}{stat_suffix}"

    label = text.replace("_angstrom", "").replace("_", " ").title()
    label = label.replace(" Ca ", " CA ").replace(" Rmsd", " RMSD").replace("Sd", "SD")
    return f"{label}{stat_suffix}"


def _feature_specs(feature_cols: list[str]) -> list[tuple[str, str, str]]:
    specs: list[tuple[str, str, str]] = []
    for idx, feature in enumerate(feature_cols):
        specs.append((feature, _pretty_feature_label(feature), FEATURE_COLORS[idx % len(FEATURE_COLORS)]))
    return specs


def _augment_with_wt_row(
    feat: pd.DataFrame,
    *,
    frame_feature_csv: Path,
    mmgbsa_replicate_csv: Path,
    susceptibility_xlsx: Path,
) -> pd.DataFrame:
    if "WT" in feat["mutation"].astype(str).tolist():
        return feat.copy()
    frame_df = pd.read_csv(frame_feature_csv)
    mmgbsa_df = pd.read_csv(mmgbsa_replicate_csv) if mmgbsa_replicate_csv.exists() else None
    target_df = load_dor_susceptibilities(susceptibility_xlsx)
    wt_row = pd.DataFrame(
        [{"drug": "DOR", "mutation": "WT", "chain": "A", "dor_fold_reduction": 1.0, "order": -1}]
    )
    target_aug = pd.concat([target_df, wt_row], ignore_index=True)
    full_aug = _mutation_feature_matrix(
        frame_df,
        target_df=target_aug,
        temperature_k=300.0,
        dispersion_mode="replicate_sd",
        mmgbsa_df=mmgbsa_df,
    )
    if "WT" not in full_aug["mutation"].astype(str).tolist():
        return feat.copy()
    wt_aug = full_aug[full_aug["mutation"].astype(str) == "WT"].copy()
    keep_cols = [c for c in feat.columns if c in wt_aug.columns]
    wt_aug = wt_aug[keep_cols]
    out = pd.concat([feat, wt_aug], ignore_index=True)
    if "target_fold_reduction" in out.columns:
        out = out.sort_values("target_fold_reduction", ascending=True, kind="stable").reset_index(drop=True)
    return out


def _feature_base_and_kind(feature_name: str) -> tuple[str, str]:
    text = str(feature_name)
    if text.endswith("_repstd"):
        return text[: -len("_repstd")], "repstd"
    if text.endswith("_mean"):
        return text[: -len("_mean")], "mean"
    if text.endswith("_sd"):
        return text[: -len("_sd")], "sd"
    if text.endswith("_std"):
        return text[: -len("_std")], "std"
    return text, "raw"


def _sem(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size <= 1:
        return 0.0
    return float(np.std(arr, ddof=1) / np.sqrt(arr.size))


def _sd_chi_square_ci(
    values: np.ndarray,
    *,
    ci_level: float,
) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size <= 1:
        point = 0.0
        return point, point
    sample_sd = float(np.std(arr, ddof=1))
    if sample_sd <= 0.0:
        return 0.0, 0.0
    df = int(arr.size - 1)
    alpha = 1.0 - float(ci_level)
    chi2_low = float(stats.chi2.ppf(alpha / 2.0, df))
    chi2_high = float(stats.chi2.ppf(1.0 - alpha / 2.0, df))
    if chi2_low <= 0.0 or chi2_high <= 0.0:
        return 0.0, np.nan
    var = sample_sd**2
    ci_low = np.sqrt(df * var / chi2_high)
    ci_high = np.sqrt(df * var / chi2_low)
    alpha = 1.0 - float(ci_level)
    del alpha
    return float(ci_low), float(ci_high)


def _prepare_replicate_sources(
    frame_df: pd.DataFrame,
    mmgbsa_df: pd.DataFrame | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    structural_cols = [
        c
        for c in frame_df.columns
        if c.startswith("pocket_") or c.startswith("ligand_") or c.startswith("contact_") or c.startswith("residue_")
    ]
    structural_rep = frame_df.groupby(["mutation", "replicate"], as_index=False)[structural_cols].mean()
    if mmgbsa_df is None or mmgbsa_df.empty:
        return structural_rep, pd.DataFrame(columns=["mutation", "replicate"])
    energy_cols = [
        c
        for c in mmgbsa_df.columns
        if c.startswith("binding_dg") and not c.endswith("_std") and not c.endswith("_sem")
    ]
    mmgbsa_rep = mmgbsa_df[["mutation", "replicate", *energy_cols]].copy()
    return structural_rep, mmgbsa_rep


def _replicate_feature_values(
    feature_name: str,
    *,
    structural_rep: pd.DataFrame,
    mmgbsa_rep: pd.DataFrame,
) -> pd.DataFrame:
    base_name, _kind = _feature_base_and_kind(feature_name)
    if base_name.startswith("binding_dg"):
        if mmgbsa_rep.empty or base_name not in mmgbsa_rep.columns:
            return pd.DataFrame(columns=["mutation", "replicate", "replicate_value"])
        return mmgbsa_rep[["mutation", "replicate", base_name]].rename(columns={base_name: "replicate_value"}).copy()
    if base_name not in structural_rep.columns:
        return pd.DataFrame(columns=["mutation", "replicate", "replicate_value"])
    return structural_rep[["mutation", "replicate", base_name]].rename(columns={base_name: "replicate_value"}).copy()


def _feature_uncertainty_summary(
    feat: pd.DataFrame,
    feature_cols: list[str],
    *,
    structural_rep: pd.DataFrame,
    mmgbsa_rep: pd.DataFrame,
    ci_level: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for feature in feature_cols:
        base_name, kind = _feature_base_and_kind(feature)
        rep_df = _replicate_feature_values(feature, structural_rep=structural_rep, mmgbsa_rep=mmgbsa_rep)
        if rep_df.empty:
            continue
        for mutation, sub in rep_df.groupby("mutation", sort=True):
            values = pd.to_numeric(sub["replicate_value"], errors="coerce").dropna().to_numpy(dtype=float)
            if values.size == 0:
                continue
            matrix_match = feat.loc[feat["mutation"].astype(str) == str(mutation), feature]
            feature_value = float(pd.to_numeric(matrix_match, errors="coerce").iloc[0]) if not matrix_match.empty else np.nan
            n_reps = int(values.size)
            if kind == "mean":
                sem = _sem(values)
                rows.append(
                    {
                        "feature": str(feature),
                        "base_feature": str(base_name),
                        "mutation": str(mutation),
                        "n_replicates": n_reps,
                        "feature_value": float(feature_value),
                        "error_method": "replicate_sem",
                        "sem": float(sem),
                        "ci_level": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "yerr_lower": float(sem),
                        "yerr_upper": float(sem),
                    }
                )
            elif kind == "repstd":
                point = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
                ci_low, ci_high = _sd_chi_square_ci(
                    values,
                    ci_level=float(ci_level),
                )
                rows.append(
                    {
                        "feature": str(feature),
                        "base_feature": str(base_name),
                        "mutation": str(mutation),
                        "n_replicates": n_reps,
                        "feature_value": float(feature_value if np.isfinite(feature_value) else point),
                        "error_method": f"replicate_sd_chi_square_{int(round(float(ci_level) * 100.0))}ci",
                        "sem": np.nan,
                        "ci_level": float(ci_level),
                        "ci_low": float(ci_low),
                        "ci_high": float(ci_high),
                        "yerr_lower": float(max(0.0, point - ci_low)),
                        "yerr_upper": float(max(0.0, ci_high - point)),
                    }
                )
            else:
                rows.append(
                    {
                        "feature": str(feature),
                        "base_feature": str(base_name),
                        "mutation": str(mutation),
                        "n_replicates": n_reps,
                        "feature_value": float(feature_value),
                        "error_method": "none",
                        "sem": np.nan,
                        "ci_level": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "yerr_lower": 0.0,
                        "yerr_upper": 0.0,
                    }
                )
    return pd.DataFrame(rows)


def _feature_fold_stats(feat: pd.DataFrame, specs: list[tuple[str, str, str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    df = feat.copy()
    df["fold_reduction"] = pd.to_numeric(df["target_fold_reduction"], errors="coerce")
    df["log10_fold_reduction"] = np.log10(df["fold_reduction"])
    df["is_combo"] = df["mutation"].astype(str).str.contains(r"\+")
    for feature, label, _color in specs:
        sub = df[["mutation", "fold_reduction", "log10_fold_reduction", feature, "is_combo"]].dropna().copy()
        if len(sub) >= 3:
            pearson_r, pearson_p = stats.pearsonr(sub["log10_fold_reduction"], sub[feature])
            spearman_rho, spearman_p = stats.spearmanr(sub["log10_fold_reduction"], sub[feature])
            slope, intercept, r_value, p_value, _stderr = stats.linregress(sub["log10_fold_reduction"], sub[feature])
        else:
            pearson_r = pearson_p = spearman_rho = spearman_p = slope = intercept = r_value = p_value = np.nan
        rows.append(
            {
                "feature": feature,
                "label": label,
                "n_mutations": int(len(sub)),
                "pearson_r": float(pearson_r) if np.isfinite(pearson_r) else np.nan,
                "pearson_pvalue": float(pearson_p) if np.isfinite(pearson_p) else np.nan,
                "spearman_rho": float(spearman_rho) if np.isfinite(spearman_rho) else np.nan,
                "spearman_pvalue": float(spearman_p) if np.isfinite(spearman_p) else np.nan,
                "r_squared": float(r_value**2) if np.isfinite(r_value) else np.nan,
                "slope_per_log10_fold": float(slope) if np.isfinite(slope) else np.nan,
                "intercept": float(intercept) if np.isfinite(intercept) else np.nan,
                "n_singles": int((~sub["is_combo"]).sum()),
                "n_combos": int(sub["is_combo"].sum()),
            }
        )
    return pd.DataFrame(rows)


def _feature_replicate_points(
    feat: pd.DataFrame,
    feature_cols: list[str],
    *,
    structural_rep: pd.DataFrame,
    mmgbsa_rep: pd.DataFrame,
) -> pd.DataFrame:
    base = feat[["mutation", "target_fold_reduction"]].copy()
    base["fold_reduction"] = pd.to_numeric(base["target_fold_reduction"], errors="coerce").astype(float)
    base["is_combo"] = base["mutation"].astype(str).str.contains(r"\+")
    rows: list[pd.DataFrame] = []
    for feature in feature_cols:
        base_name, kind = _feature_base_and_kind(feature)
        if kind != "mean":
            continue
        rep_df = _replicate_feature_values(feature, structural_rep=structural_rep, mmgbsa_rep=mmgbsa_rep)
        if rep_df.empty:
            continue
        merged = rep_df.merge(base[["mutation", "fold_reduction", "is_combo"]], on="mutation", how="left")
        merged["feature"] = str(feature)
        merged["base_feature"] = str(base_name)
        rows.append(merged[["feature", "base_feature", "mutation", "replicate", "fold_reduction", "is_combo", "replicate_value"]].copy())
    if not rows:
        return pd.DataFrame(columns=["feature", "base_feature", "mutation", "replicate", "fold_reduction", "is_combo", "replicate_value"])
    return pd.concat(rows, ignore_index=True)


def _replicate_x_positions(fold_values: np.ndarray, replicates: np.ndarray) -> np.ndarray:
    x = np.asarray(fold_values, dtype=float).copy()
    reps = np.asarray(replicates)
    unique_reps = sorted(pd.unique(reps).tolist())
    if len(unique_reps) <= 1:
        return x
    offsets = np.linspace(-0.035, 0.035, num=len(unique_reps))
    offset_map = {rep: float(offset) for rep, offset in zip(unique_reps, offsets)}
    factors = np.array([1.0 + offset_map[rep] for rep in reps], dtype=float)
    return x * factors


def _plot_feature_vs_fold_change(
    feat: pd.DataFrame,
    uncertainty_df: pd.DataFrame,
    replicate_points_df: pd.DataFrame,
    feature: str,
    label: str,
    color: str,
    output_png: Path,
) -> None:
    df = feat.copy()
    df["fold_reduction"] = pd.to_numeric(df["target_fold_reduction"], errors="coerce")
    df["log10_fold_reduction"] = np.log10(df["fold_reduction"])
    df["is_combo"] = df["mutation"].astype(str).str.contains(r"\+")
    feature_unc = uncertainty_df[uncertainty_df["feature"].astype(str) == str(feature)].copy()
    if not feature_unc.empty:
        df = df.merge(
            feature_unc[
                [
                    "mutation",
                    "error_method",
                    "sem",
                    "ci_low",
                    "ci_high",
                    "yerr_lower",
                    "yerr_upper",
                    "n_replicates",
                ]
            ],
            on="mutation",
            how="left",
        )
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log10_fold_reduction", feature]).reset_index(drop=True)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14.2, 9.2))
    _base_name, kind = _feature_base_and_kind(feature)
    singles = df[~df["is_combo"]].copy()
    combos = df[df["is_combo"]].copy()
    for subset, marker, legend_label in [
        (singles, "o", "Single DRM"),
        (combos, "s", "Combination DRM"),
    ]:
        if subset.empty:
            continue
        if kind == "mean":
            err_subset = subset.dropna(subset=["yerr_lower", "yerr_upper"]).copy()
            for row in err_subset.itertuples(index=False):
                ax.errorbar(
                    float(row.fold_reduction),
                    float(getattr(row, feature)),
                    yerr=np.array([[float(row.yerr_lower)], [float(row.yerr_upper)]]),
                    fmt="none",
                    ecolor=color,
                    elinewidth=1.2,
                    capsize=3.0,
                    alpha=0.75,
                    zorder=2,
                )
        ax.scatter(
            subset["fold_reduction"],
            subset[feature],
            s=78,
            color=color,
            marker=marker,
            edgecolors="white",
            linewidths=0.8,
            alpha=0.95,
            label=legend_label,
            zorder=3,
        )

    valid = df[["log10_fold_reduction", "fold_reduction", feature]].dropna()
    if len(valid) >= 3:
        slope, intercept, r_value, p_value, _stderr = stats.linregress(valid["log10_fold_reduction"], valid[feature])
        pearson_r, pearson_p = stats.pearsonr(valid["log10_fold_reduction"], valid[feature])
        spearman_rho, spearman_p = stats.spearmanr(valid["log10_fold_reduction"], valid[feature])
        x_grid = np.geomspace(float(df["fold_reduction"].min()) * 0.9, float(df["fold_reduction"].max()) * 1.1, 300)
        y_grid = slope * np.log10(x_grid) + intercept
        ax.plot(x_grid, y_grid, color="#444444", linestyle="--", linewidth=1.5, zorder=2)
        annotation = (
            f"Pearson r = {pearson_r:.3f}\n"
            f"Spearman rho = {spearman_rho:.3f}\n"
            f"R^2 = {r_value**2:.3f}\n"
            f"p = {p_value:.3f}"
        )
    else:
        annotation = "Pearson r = NA\nSpearman rho = NA\nR^2 = NA\np = NA"
    if kind == "mean":
        annotation += "\nError bars: replicate SEM"
    elif kind == "repstd":
        annotation += "\nNo error bars: summary replicate SD"

    ax.set_xscale("log")
    ax.grid(alpha=0.25)
    ax.set_title(f"{label} vs Fold-change", fontsize=14, fontweight="bold")
    ax.set_xlabel("Fold-change")
    ax.set_ylabel(label)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, fontsize=10, frameon=False)
    _place_greedy_annotations(
        ax,
        df["fold_reduction"].to_numpy(dtype=float),
        df[feature].to_numpy(dtype=float),
        df["mutation"].astype(str).tolist(),
    )
    ax.text(
        0.02,
        0.98,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot structural logistic features vs susceptibility fold change.")
    parser.add_argument(
        "--feature-matrix-csv",
        type=Path,
        default=Path("results/analysis/logistic_regression_contrib_ge_1p5/feature_screening/tables/mutation_feature_matrix.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/logistic_regression_contrib_ge_1p5/feature_vs_fold"),
    )
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--mmgbsa-replicate-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/last20frames/mmgbsa_replicate_metrics_last20frames.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument("--include-wt-reference", action="store_true")
    parser.add_argument("--ci-level", type=float, default=0.95)
    args = parser.parse_args()

    if not args.feature_matrix_csv.exists():
        raise FileNotFoundError(args.feature_matrix_csv)
    if not args.frame_feature_csv.exists():
        raise FileNotFoundError(args.frame_feature_csv)
    if args.include_wt_reference and not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feat = pd.read_csv(args.feature_matrix_csv)
    if args.include_wt_reference:
        feat = _augment_with_wt_row(
            feat,
            frame_feature_csv=args.frame_feature_csv,
            mmgbsa_replicate_csv=args.mmgbsa_replicate_csv,
            susceptibility_xlsx=args.susceptibility_xlsx,
        )
    frame_df = pd.read_csv(args.frame_feature_csv)
    mmgbsa_df = pd.read_csv(args.mmgbsa_replicate_csv) if args.mmgbsa_replicate_csv.exists() else pd.DataFrame()
    feature_cols = [c for c in feat.columns if c not in META_COLUMNS]
    specs = _feature_specs(feature_cols)
    structural_rep, mmgbsa_rep = _prepare_replicate_sources(frame_df, mmgbsa_df)
    uncertainty_df = _feature_uncertainty_summary(
        feat,
        feature_cols,
        structural_rep=structural_rep,
        mmgbsa_rep=mmgbsa_rep,
        ci_level=float(args.ci_level),
    ).sort_values(["feature", "mutation"], kind="stable").reset_index(drop=True)
    replicate_points_df = _feature_replicate_points(
        feat,
        feature_cols,
        structural_rep=structural_rep,
        mmgbsa_rep=mmgbsa_rep,
    ).sort_values(["feature", "mutation", "replicate"], kind="stable").reset_index(drop=True)

    stats_df = _feature_fold_stats(feat, specs).sort_values("r_squared", ascending=False, kind="stable").reset_index(drop=True)
    stats_df.to_csv(out_tables / "feature_vs_fold_change_stats.csv", index=False)
    uncertainty_df.to_csv(out_tables / "feature_uncertainty_by_mutation.csv", index=False)
    replicate_points_df.to_csv(out_tables / "feature_replicate_points.csv", index=False)
    feat.to_csv(out_tables / "mutation_feature_matrix.csv", index=False)

    for feature, label, color in specs:
        _plot_feature_vs_fold_change(
            feat,
            uncertainty_df,
            replicate_points_df,
            feature,
            label,
            color,
            out_plots / f"{_feature_slug(feature)}_vs_fold_change.png",
        )

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "feature_matrix_csv": str(args.feature_matrix_csv),
                "output_dir": str(args.output_dir),
                "frame_feature_csv": str(args.frame_feature_csv),
                "mmgbsa_replicate_csv": str(args.mmgbsa_replicate_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "include_wt_reference": bool(args.include_wt_reference),
                "repstd_ci_method": "exact_chi_square",
                "ci_level": float(args.ci_level),
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
