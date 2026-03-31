#!/usr/bin/env python3
"""Engineer event-like mutation features and screen them on curated controls."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .search_feature_combo_logistic_controls import NEGATIVE_CONTROLS, POSITIVE_CONTROLS, UNCERTAIN_LIMITED
from ..susceptibility import load_dor_susceptibilities


META_COLUMNS = {"drug", "mutation", "chain", "target_fold_reduction", "control_category"}


def _category_for_mutation(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation in NEGATIVE_CONTROLS:
        return "negative_control"
    if mutation in POSITIVE_CONTROLS:
        return "positive_control"
    if mutation in UNCERTAIN_LIMITED:
        return "uncertain_limited"
    if mutation == "WT":
        return "wt_reference"
    raise ValueError(f"Mutation is not assigned to a control category: {label}")


def _load_targets(path: Path) -> pd.DataFrame:
    target_df = load_dor_susceptibilities(path).copy()
    target_df["target_fold_reduction"] = pd.to_numeric(target_df["dor_fold_reduction"], errors="coerce").astype(float)
    target_df["control_category"] = target_df["mutation"].astype(str).map(_category_for_mutation)
    keep_cols = [c for c in ("drug", "mutation", "chain", "target_fold_reduction", "control_category") if c in target_df.columns]
    return target_df[keep_cols].copy()


def _build_frame_event_features(frame_df: pd.DataFrame) -> pd.DataFrame:
    base_cols = [
        c
        for c in frame_df.columns
        if c.startswith("pocket_") or c.startswith("ligand_") or c.startswith("residue_")
    ]
    rows: list[dict[str, object]] = []
    for (mutation, replicate), sub in frame_df.groupby(["mutation", "replicate"], sort=True):
        row: dict[str, object] = {"mutation": str(mutation), "replicate": int(replicate)}
        for col in base_cols:
            values = pd.to_numeric(sub[col], errors="coerce").dropna().astype(float)
            if values.empty:
                continue
            row[f"{col}_mean"] = float(values.mean())
            row[f"{col}_median"] = float(values.median())
            row[f"{col}_q10"] = float(values.quantile(0.10))
            row[f"{col}_q90"] = float(values.quantile(0.90))
            if col.startswith("residue_min_distance_"):
                row[f"{col}_occ_lt_3p5"] = float((values < 3.5).mean())
                row[f"{col}_occ_lt_4p0"] = float((values < 4.0).mean())
            elif col == "ligand_pose_rmsd_angstrom":
                row[f"{col}_occ_gt_1p5"] = float((values > 1.5).mean())
                row[f"{col}_occ_gt_2p0"] = float((values > 2.0).mean())
            elif col == "pocket_ca_rmsd_angstrom":
                row[f"{col}_occ_gt_0p75"] = float((values > 0.75).mean())
                row[f"{col}_occ_gt_1p0"] = float((values > 1.0).mean())
        rows.append(row)
    rep_df = pd.DataFrame(rows)
    feature_cols = [c for c in rep_df.columns if c not in {"mutation", "replicate"}]
    mut_df = rep_df.groupby("mutation", as_index=False)[feature_cols].mean()
    return mut_df


def _build_mmgbsa_features(mmgbsa_df: pd.DataFrame) -> pd.DataFrame:
    if mmgbsa_df.empty:
        return pd.DataFrame(columns=["mutation"])
    energy_cols = [
        "binding_dg",
        "binding_dg_vdw",
        "binding_dg_electrostatic",
        "binding_dg_gb",
        "binding_dg_sa",
    ]
    avail = [c for c in energy_cols if c in mmgbsa_df.columns]
    if not avail:
        return pd.DataFrame(columns=["mutation"])
    rep = mmgbsa_df[["mutation", "replicate", *avail]].copy()
    rows = []
    for mutation, sub in rep.groupby("mutation", sort=True):
        row = {"mutation": str(mutation)}
        for col in avail:
            values = pd.to_numeric(sub[col], errors="coerce").dropna().astype(float)
            if values.empty:
                continue
            row[f"{col}_mean"] = float(values.mean())
            row[f"{col}_median"] = float(values.median())
        rows.append(row)
    return pd.DataFrame(rows)


def _cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    total = 0
    for xi in x:
        total += np.sum(xi > y) - np.sum(xi < y)
    return float(total / (len(x) * len(y)))


def _safe_corr(x: pd.Series, y: pd.Series) -> tuple[float, float, float, float]:
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return (float("nan"),) * 4
    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_rho, spearman_p = stats.spearmanr(x, y)
    return float(pearson_r), float(pearson_p), float(spearman_rho), float(spearman_p)


def _single_feature_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="l2",
                    C=1.0,
                    solver="lbfgs",
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=0,
                ),
            ),
        ]
    )


def _single_feature_loocv(df: pd.DataFrame, *, feature: str) -> dict[str, float]:
    x = pd.to_numeric(df[feature], errors="coerce").astype(float).to_numpy(dtype=float).reshape(-1, 1)
    y = df["control_category"].astype(str).map(lambda value: 1 if value == "positive_control" else 0).to_numpy(dtype=int)
    loo = LeaveOneOut()
    prob = np.full(len(df), np.nan, dtype=float)
    pred = np.full(len(df), -1, dtype=int)
    for train_idx, test_idx in loo.split(x, y):
        pipe = _single_feature_pipeline()
        pipe.fit(x[train_idx], y[train_idx])
        prob[test_idx[0]] = float(pipe.predict_proba(x[test_idx])[:, 1][0])
        pred[test_idx[0]] = int(prob[test_idx[0]] >= 0.5)
    prob_clip = np.clip(prob, 1e-6, 1 - 1e-6)
    logit = np.log(prob_clip / (1.0 - prob_clip))
    y_fold = pd.to_numeric(df["target_fold_reduction"], errors="coerce").astype(float)
    raw_pr, raw_pp, raw_sr, raw_sp = _safe_corr(pd.Series(logit), y_fold)
    log_pr, log_pp, log_sr, log_sp = _safe_corr(pd.Series(logit), y_fold.map(math.log10))
    return {
        "single_feature_loocv_accuracy": float(accuracy_score(y, pred)),
        "single_feature_loocv_balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "single_feature_loocv_roc_auc": float(roc_auc_score(y, prob)),
        "single_feature_loocv_logit_raw_pearson_r": raw_pr,
        "single_feature_loocv_logit_raw_pearson_pvalue": raw_pp,
        "single_feature_loocv_logit_raw_spearman_rho": raw_sr,
        "single_feature_loocv_logit_raw_spearman_pvalue": raw_sp,
        "single_feature_loocv_logit_log10_pearson_r": log_pr,
        "single_feature_loocv_logit_log10_pearson_pvalue": log_pp,
        "single_feature_loocv_logit_log10_spearman_rho": log_sr,
        "single_feature_loocv_logit_log10_spearman_pvalue": log_sp,
    }


def _screen_features(mut_feat: pd.DataFrame) -> pd.DataFrame:
    controls = mut_feat[mut_feat["control_category"].isin({"negative_control", "positive_control"})].copy()
    feature_cols = [c for c in controls.columns if c not in META_COLUMNS]
    rows: list[dict[str, object]] = []
    for feature in feature_cols:
        values = pd.to_numeric(controls[feature], errors="coerce").astype(float)
        work = controls[["mutation", "control_category", "target_fold_reduction"]].copy()
        work["feature_value"] = values
        work = work.dropna().copy()
        if len(work) < 6:
            continue
        neg = work[work["control_category"] == "negative_control"]["feature_value"].to_numpy(dtype=float)
        pos = work[work["control_category"] == "positive_control"]["feature_value"].to_numpy(dtype=float)
        if len(np.unique(work["feature_value"])) < 2:
            continue
        auc = roc_auc_score((work["control_category"] == "positive_control").astype(int), work["feature_value"])
        oriented_sign = 1.0 if auc >= 0.5 else -1.0
        auc_abs = float(max(auc, 1.0 - auc))
        mw = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        cliff = _cliffs_delta(pos, neg)
        if oriented_sign > 0:
            class_gap = float(np.min(pos) - np.max(neg))
        else:
            class_gap = float(np.min(neg) - np.max(pos))
        neg_mean = float(np.mean(neg))
        pos_mean = float(np.mean(pos))
        pooled_sd = float(np.sqrt((((len(neg) - 1) * np.var(neg, ddof=1)) + ((len(pos) - 1) * np.var(pos, ddof=1))) / max(1, len(neg) + len(pos) - 2)))
        cohen_d = float((pos_mean - neg_mean) / pooled_sd) if pooled_sd > 0 else float("nan")
        raw_pr, raw_pp, raw_sr, raw_sp = _safe_corr(work["feature_value"], work["target_fold_reduction"])
        log_pr, log_pp, log_sr, log_sp = _safe_corr(work["feature_value"], work["target_fold_reduction"].map(math.log10))
        loocv = _single_feature_loocv(controls[["mutation", "control_category", "target_fold_reduction", feature]].copy(), feature=feature)
        rows.append(
            {
                "feature": str(feature),
                "n_controls": int(len(work)),
                "orientation": "higher_in_positive" if oriented_sign > 0 else "lower_in_positive",
                "negative_mean": neg_mean,
                "positive_mean": pos_mean,
                "class_gap": class_gap,
                "cohen_d": cohen_d,
                "cliffs_delta": cliff,
                "mannwhitney_u_pvalue": float(mw.pvalue),
                "single_feature_abs_auroc": auc_abs,
                "raw_fold_pearson_r": raw_pr,
                "raw_fold_pearson_pvalue": raw_pp,
                "raw_fold_spearman_rho": raw_sr,
                "raw_fold_spearman_pvalue": raw_sp,
                "log10_fold_pearson_r": log_pr,
                "log10_fold_pearson_pvalue": log_pp,
                "log10_fold_spearman_rho": log_sr,
                "log10_fold_spearman_pvalue": log_sp,
                **loocv,
            }
        )
    out = pd.DataFrame(rows)
    out["screen_rank_score"] = (
        0.45 * out["single_feature_abs_auroc"].astype(float)
        + 0.35 * out["single_feature_loocv_logit_log10_pearson_r"].astype(float).clip(lower=-1, upper=1).fillna(0.0)
        + 0.20 * out["single_feature_loocv_balanced_accuracy"].astype(float)
    )
    return out.sort_values(
        [
            "screen_rank_score",
            "single_feature_abs_auroc",
            "single_feature_loocv_logit_log10_pearson_r",
            "class_gap",
        ],
        ascending=[False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)


def _feature_value_by_mutation(mut_feat: pd.DataFrame, feature: str) -> pd.DataFrame:
    cols = ["mutation", "control_category", "target_fold_reduction", feature]
    out = mut_feat[cols].copy().rename(columns={feature: "feature_value"})
    return out.sort_values(["control_category", "target_fold_reduction", "mutation"], kind="stable").reset_index(drop=True)


def _plot_top_features(mut_feat: pd.DataFrame, ranking_df: pd.DataFrame, output_dir: Path, top_n: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    shortlist = ranking_df.head(int(top_n))
    for record in shortlist.itertuples(index=False):
        feature = str(record.feature)
        df = _feature_value_by_mutation(mut_feat, feature)
        controls = df[df["control_category"].isin({"negative_control", "positive_control"})].copy()
        controls["x"] = controls["control_category"].map({"negative_control": 0, "positive_control": 1}).astype(float)
        jitter = np.linspace(-0.08, 0.08, len(controls))
        controls = controls.sort_values(["control_category", "feature_value", "mutation"], kind="stable").reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(6.6, 4.4))
        for xpos, label, color in [(0, "Negative", "#1f77b4"), (1, "Positive", "#d62728")]:
            sub = controls[controls["x"] == xpos].copy()
            if sub.empty:
                continue
            ax.scatter(
                np.full(len(sub), xpos) + jitter[: len(sub)],
                sub["feature_value"],
                s=52,
                color=color,
                alpha=0.9,
                edgecolor="white",
                linewidth=0.6,
                label=label,
            )
            q1 = float(sub["feature_value"].quantile(0.25))
            med = float(sub["feature_value"].median())
            q3 = float(sub["feature_value"].quantile(0.75))
            ax.vlines(xpos, q1, q3, color=color, linewidth=5, alpha=0.35)
            ax.hlines(med, xpos - 0.16, xpos + 0.16, color=color, linewidth=2.4)
        ax.set_xticks([0, 1], ["Negative", "Positive"])
        ax.set_ylabel(feature)
        ax.set_title(
            f"{feature}\nAUROC={record.single_feature_abs_auroc:.3f} | "
            f"LOOCV logit-log10 r={record.single_feature_loocv_logit_log10_pearson_r:.3f}"
        )
        ax.grid(axis="y", alpha=0.22)
        fig.tight_layout()
        out = output_dir / f"{feature}.png"
        fig.savefig(out, dpi=240, bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Screen event-like mutation features on curated controls.")
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
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/control_event_feature_screen"),
    )
    parser.add_argument("--top-n-plots", type=int, default=12)
    args = parser.parse_args()

    frame_df = pd.read_csv(args.frame_feature_csv)
    mmgbsa_df = pd.read_csv(args.mmgbsa_replicate_csv) if args.mmgbsa_replicate_csv.exists() else pd.DataFrame()
    target_df = _load_targets(args.susceptibility_xlsx)

    frame_feat = _build_frame_event_features(frame_df)
    energy_feat = _build_mmgbsa_features(mmgbsa_df)
    mut_feat = target_df.merge(frame_feat, on="mutation", how="left")
    if not energy_feat.empty:
        mut_feat = mut_feat.merge(energy_feat, on="mutation", how="left")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    ranking = _screen_features(mut_feat)
    shortlist = ranking[
        (ranking["single_feature_abs_auroc"] >= 0.90)
        & (ranking["single_feature_loocv_logit_log10_pearson_r"] >= 0.55)
    ].copy()
    ranking.to_csv(out_tables / "feature_screen_ranking.csv", index=False)
    shortlist.to_csv(out_tables / "feature_screen_shortlist.csv", index=False)
    mut_feat.to_csv(out_tables / "engineered_feature_matrix.csv", index=False)
    _plot_top_features(mut_feat, ranking, out_plots, top_n=int(args.top_n_plots))

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "frame_feature_csv": str(args.frame_feature_csv),
                "mmgbsa_replicate_csv": str(args.mmgbsa_replicate_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "top_n_plots": int(args.top_n_plots),
                "feature_families": {
                    "generic_stats": ["mean", "median", "q10", "q90"],
                    "residue_contact_occupancy_thresholds_angstrom": [3.5, 4.0],
                    "ligand_pose_rmsd_thresholds_angstrom": [1.5, 2.0],
                    "pocket_ca_rmsd_thresholds_angstrom": [0.75, 1.0],
                },
                "shortlist_rule": {
                    "single_feature_abs_auroc_min": 0.90,
                    "single_feature_loocv_logit_log10_pearson_r_min": 0.55,
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
