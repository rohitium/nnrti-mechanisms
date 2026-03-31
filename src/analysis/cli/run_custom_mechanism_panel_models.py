#!/usr/bin/env python3
"""Build and evaluate a custom mechanism-driven feature panel for DOR susceptibility."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import pearsonr, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Lasso, LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_dor_susceptibility_bars import NEGATIVE_CONTROLS, POSITIVE_CONTROLS, UNCERTAIN_LIMITED
from .plot_triplet_contact_story import _load_replicate_meta


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run custom mechanism-panel models for DOR susceptibility.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument(
        "--frame-features-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/custom_mechanism_panel_models"),
    )
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=5)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    return parser.parse_args()


def _control_category(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation == "WT":
        return "wt_reference"
    if mutation in NEGATIVE_CONTROLS:
        return "negative_control"
    if mutation in POSITIVE_CONTROLS:
        return "positive_control"
    if mutation in UNCERTAIN_LIMITED:
        return "uncertain_limited"
    return "other"


def _binary_class(fold: float, low_max_fold: float) -> str:
    return "low" if float(fold) < float(low_max_fold) else "high"


def _sidechain_atoms(residue):
    sc = residue.atoms.select_atoms("not name N CA C O OXT and not name H*")
    if sc.n_atoms == 0:
        sc = residue.atoms.select_atoms("not name H*")
    return sc


def _compute_custom_replicate_means(
    metas,
    ligand_resname: str,
    resid_offset: int,
    frame_stride: int,
) -> pd.DataFrame:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    traj_map = {
        "SER105": int(105 + resid_offset),
        "VAL106": int(106 + resid_offset),
        "VAL179": int(179 + resid_offset),
        "TYR188": int(188 + resid_offset),
        "PRO225": int(225 + resid_offset),
        "PHE227": int(227 + resid_offset),
    }

    rows: list[dict[str, object]] = []
    for meta in metas:
        u = mda.Universe(str(meta.topology_pdb), str(meta.analysis_dcd), format="DCD")
        lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if lig.n_atoms == 0:
            continue

        residue_atoms: dict[str, object] = {}
        for name, traj_resid in traj_map.items():
            sel = u.select_atoms(f"protein and resid {traj_resid} and not name H*")
            if sel.n_atoms == 0:
                residue_atoms[name] = None
                continue
            residue_atoms[name] = sel

        v106_res = u.select_atoms(f"protein and resid {traj_map['VAL106']}")
        f227_res = u.select_atoms(f"protein and resid {traj_map['PHE227']}")
        v106_sc = _sidechain_atoms(v106_res.residues[0]) if v106_res.n_atoms else None
        f227_sc = _sidechain_atoms(f227_res.residues[0]) if f227_res.n_atoms else None

        vals = {
            "ser105_dor_distance_angstrom": [],
            "v106_f227_sidechain_distance_angstrom": [],
        }
        for ts in u.trajectory[:: max(1, int(frame_stride))]:
            box = u.dimensions
            ser = residue_atoms["SER105"]
            if ser is not None and ser.n_atoms > 0:
                vals["ser105_dor_distance_angstrom"].append(
                    float(distance_array(ser.positions, lig.positions, box=box).min())
                )
            if v106_sc is not None and f227_sc is not None and v106_sc.n_atoms > 0 and f227_sc.n_atoms > 0:
                vals["v106_f227_sidechain_distance_angstrom"].append(
                    float(distance_array(v106_sc.positions, f227_sc.positions, box=box).min())
                )

        rows.append(
            {
                "mutation": str(meta.mutation),
                "replicate": int(meta.replicate),
                "ser105_dor_distance_angstrom_mean": float(np.nanmean(vals["ser105_dor_distance_angstrom"])),
                "v106_f227_sidechain_distance_angstrom_mean": float(np.nanmean(vals["v106_f227_sidechain_distance_angstrom"])),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No custom replicate metrics were computed.")
    return out


def _aggregate_frame_feature_means(frame_features_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(frame_features_csv).copy()
    wanted = [
        "residue_min_distance_VAL106_angstrom",
        "residue_min_distance_VAL179_angstrom",
        "residue_min_distance_TYR188_angstrom",
        "residue_min_distance_PRO225_angstrom",
        "residue_min_distance_PHE227_angstrom",
        "ligand_pose_rmsd_angstrom",
        "ligand_palm_distance_angstrom",
        "ligand_entrance_distance_angstrom",
    ]
    rep = (
        df.groupby(["mutation", "replicate"], as_index=False)[wanted]
        .mean()
        .rename(columns={c: f"{c}_repmean" for c in wanted})
    )
    agg = rep.groupby("mutation", as_index=False).agg(
        **{
            f"{c}_mean": (f"{c}_repmean", "mean")
            for c in wanted
        }
    )
    return agg


def _logistic_pipeline(penalty: str, c_value: float) -> Pipeline:
    solver = "liblinear" if penalty == "l1" else "lbfgs"
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(penalty=penalty, C=float(c_value), solver=solver, max_iter=5000, random_state=0)),
        ]
    )


def _regression_pipeline(model) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def main() -> int:
    args = _parse_args()
    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)
    if not args.frame_features_csv.exists():
        raise FileNotFoundError(args.frame_features_csv)

    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    needed = set(NEGATIVE_CONTROLS) | set(POSITIVE_CONTROLS) | set(UNCERTAIN_LIMITED) | {"WT"}
    metas = _load_replicate_meta(args.manifest, needed_mutations=needed)
    rep_custom = _compute_custom_replicate_means(
        metas=metas,
        ligand_resname=str(args.ligand_resname),
        resid_offset=int(args.resid_offset),
        frame_stride=int(args.frame_stride),
    )
    rep_custom.to_csv(out_tables / "custom_replicate_feature_means.csv", index=False)

    custom_mut = rep_custom.groupby("mutation", as_index=False).agg(
        ser105_dor_distance_angstrom_mean=("ser105_dor_distance_angstrom_mean", "mean"),
        v106_f227_sidechain_distance_angstrom_mean=("v106_f227_sidechain_distance_angstrom_mean", "mean"),
    )
    custom_mut.to_csv(out_tables / "custom_mutation_feature_means.csv", index=False)

    frame_mut = _aggregate_frame_feature_means(args.frame_features_csv)
    panel_df = frame_mut.merge(custom_mut, on="mutation", how="outer")
    panel_df["control_category"] = panel_df["mutation"].map(_control_category)

    fold_map = (
        pd.read_csv("results/analysis/contact_occupancy_feature_screen/tables/occupancy_mean_feature_matrix.csv")[["mutation", "dor_fold_reduction"]]
        .drop_duplicates()
        .rename(columns={"dor_fold_reduction": "target_fold_reduction"})
    )
    panel_df = panel_df.merge(fold_map, on="mutation", how="left")
    panel_df["target_binary_class"] = panel_df["target_fold_reduction"].map(
        lambda x: _binary_class(x, float(args.low_max_fold)) if np.isfinite(x) else np.nan
    )

    panel_features = [
        "ser105_dor_distance_angstrom_mean",
        "v106_f227_sidechain_distance_angstrom_mean",
        "residue_min_distance_VAL106_angstrom_mean",
        "residue_min_distance_VAL179_angstrom_mean",
        "residue_min_distance_TYR188_angstrom_mean",
        "residue_min_distance_PRO225_angstrom_mean",
        "residue_min_distance_PHE227_angstrom_mean",
        "ligand_pose_rmsd_angstrom_mean",
        "ligand_palm_distance_angstrom_mean",
        "ligand_entrance_distance_angstrom_mean",
    ]
    panel_df[["mutation", "control_category", "target_fold_reduction", "target_binary_class"] + panel_features].to_csv(
        out_tables / "mechanism_panel_feature_matrix.csv", index=False
    )

    train = panel_df[panel_df["control_category"].isin({"negative_control", "positive_control"})].copy().reset_index(drop=True)
    hold = panel_df[panel_df["control_category"].isin({"uncertain_limited", "wt_reference"})].copy().reset_index(drop=True)
    x_train = train[panel_features].copy()
    y_train = train["target_binary_class"].astype(str).copy()

    loo = LeaveOneOut()
    log_rows: list[dict[str, object]] = []
    for penalty, c_values in [("l2", [0.03, 0.1, 0.3, 1.0, 3.0, 10.0]), ("l1", [0.03, 0.1, 0.3, 1.0, 3.0, 10.0])]:
        for c_value in c_values:
            probs = np.zeros(len(train), dtype=float)
            for tr_idx, te_idx in loo.split(x_train, y_train):
                fitted = _logistic_pipeline(penalty, c_value)
                fitted.fit(x_train.iloc[tr_idx], y_train.iloc[tr_idx])
                pos_idx = list(fitted.named_steps["model"].classes_).index("high")
                probs[te_idx[0]] = float(fitted.predict_proba(x_train.iloc[te_idx])[:, pos_idx][0])
            pred = np.where(probs >= 0.5, "high", "low")
            logits = logit(np.clip(probs, 1e-6, 1.0 - 1e-6))
            full = _logistic_pipeline(penalty, c_value)
            full.fit(x_train, y_train)
            coef = full.named_steps["model"].coef_[0]
            nz = [panel_features[i] for i, v in enumerate(coef) if abs(v) > 1e-10]
            wt_prob = np.nan
            wt_row = hold[hold["mutation"] == "WT"]
            if not wt_row.empty:
                pos_idx = list(full.named_steps["model"].classes_).index("high")
                wt_prob = float(full.predict_proba(wt_row[panel_features])[:, pos_idx][0])
            log_rows.append(
                {
                    "penalty": penalty,
                    "c_value": float(c_value),
                    "n_nonzero": int(len(nz)),
                    "features_nonzero": "|".join(nz),
                    "accuracy": float(accuracy_score(y_train, pred)),
                    "balanced_accuracy": float(balanced_accuracy_score(y_train, pred)),
                    "macro_f1": float(f1_score(y_train, pred, average="macro")),
                    "roc_auc": float(roc_auc_score((y_train == "high").astype(int), probs)),
                    "pearson_r_logit_vs_fold": float(pearsonr(logits, train["target_fold_reduction"].to_numpy(dtype=float)).statistic),
                    "pearson_r_logit_vs_log10_fold": float(
                        pearsonr(logits, np.log10(train["target_fold_reduction"].to_numpy(dtype=float))).statistic
                    ),
                    "wt_prob_high": wt_prob,
                }
            )
    log_df = pd.DataFrame(log_rows).sort_values(
        ["balanced_accuracy", "roc_auc", "macro_f1", "accuracy"],
        ascending=[False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)
    log_df.to_csv(out_tables / "control_logistic_sweep.csv", index=False)

    best = log_df.iloc[0]
    best_log = _logistic_pipeline(str(best["penalty"]), float(best["c_value"]))
    best_log.fit(x_train, y_train)
    coef = best_log.named_steps["model"].coef_[0]
    pd.DataFrame(
        {"feature": panel_features, "coefficient": coef, "abs_coefficient": np.abs(coef)}
    ).sort_values("abs_coefficient", ascending=False, kind="stable").to_csv(
        out_tables / "best_control_logistic_coefficients.csv", index=False
    )
    pos_idx = list(best_log.named_steps["model"].classes_).index("high")
    hold_out = hold[["mutation", "control_category", "target_fold_reduction"]].copy()
    hold_out["prob_high"] = best_log.predict_proba(hold[panel_features])[:, pos_idx]
    hold_out["predicted_class"] = np.where(hold_out["prob_high"] >= 0.5, "high", "low")
    hold_out.to_csv(out_tables / "best_control_logistic_holdout_predictions.csv", index=False)

    reg_df = panel_df[panel_df["mutation"] != "WT"].copy().reset_index(drop=True)
    x_reg = reg_df[panel_features].copy()
    y_reg = np.log10(reg_df["target_fold_reduction"].to_numpy(dtype=float))
    reg_rows: list[dict[str, object]] = []
    reg_specs = [
        ("ridge_0p3", Ridge(alpha=0.3)),
        ("ridge_1", Ridge(alpha=1.0)),
        ("ridge_3", Ridge(alpha=3.0)),
        ("lasso_0p01", Lasso(alpha=0.01, max_iter=20000)),
        ("lasso_0p03", Lasso(alpha=0.03, max_iter=20000)),
    ]
    for name, model in reg_specs:
        pipe = _regression_pipeline(model)
        pred = cross_val_predict(pipe, x_reg, y_reg, cv=LeaveOneOut())
        reg_rows.append(
            {
                "model": name,
                "r2": float(r2_score(y_reg, pred)),
                "mae": float(mean_absolute_error(y_reg, pred)),
                "rmse": float(np.sqrt(mean_squared_error(y_reg, pred))),
                "pearson_r": float(pearsonr(y_reg, pred).statistic),
                "spearman_rho": float(spearmanr(y_reg, pred).statistic),
            }
        )
        pd.DataFrame(
            {
                "mutation": reg_df["mutation"],
                "observed_log10_fold": y_reg,
                "predicted_log10_fold": pred,
                "residual": pred - y_reg,
            }
        ).sort_values("observed_log10_fold", kind="stable").to_csv(
            out_tables / f"{name}_log10_fold_cv_predictions.csv", index=False
        )
    reg_out = pd.DataFrame(reg_rows).sort_values(["r2", "pearson_r"], ascending=[False, False], kind="stable").reset_index(drop=True)
    reg_out.to_csv(out_tables / "all_label_regression_sweep.csv", index=False)

    wt_pred = pd.DataFrame()
    wt_row = panel_df[panel_df["mutation"] == "WT"]
    if not wt_row.empty:
        best_reg_name = str(reg_out.iloc[0]["model"])
        best_reg_model = dict(reg_specs)[best_reg_name]
        best_reg = _regression_pipeline(best_reg_model)
        best_reg.fit(x_reg, y_reg)
        wt_log10 = float(best_reg.predict(wt_row[panel_features])[0])
        wt_pred = pd.DataFrame(
            {
                "mutation": ["WT"],
                "predicted_log10_fold": [wt_log10],
                "predicted_fold": [float(10 ** wt_log10)],
            }
        )
        wt_pred.to_csv(out_tables / "wt_regression_prediction.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "frame_features_csv": str(args.frame_features_csv),
                "output_dir": str(args.output_dir),
                "ligand_resname": str(args.ligand_resname),
                "resid_offset": int(args.resid_offset),
                "frame_stride": int(args.frame_stride),
                "panel_features": panel_features,
                "training_design": "control_only_logistic_plus_all_label_regression",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
