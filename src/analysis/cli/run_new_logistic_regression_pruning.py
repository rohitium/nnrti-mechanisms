#!/usr/bin/env python3
"""Greedy backward-pruned logistic model over WT-contacted residue-DOR distances."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import linregress
from sklearn.impute import SimpleImputer
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_curve,
    roc_auc_score,
)
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .plot_dor_susceptibility_bars import NEGATIVE_CONTROLS, POSITIVE_CONTROLS, UNCERTAIN_PHENOTYPE
from .plot_triplet_contact_story import _load_replicate_meta
from ..susceptibility import load_dor_susceptibilities

warnings.filterwarnings("ignore", category=FutureWarning, module=r"sklearn\.linear_model\._logistic")
warnings.filterwarnings("ignore", category=UserWarning, module=r"sklearn\.linear_model\._logistic")
warnings.filterwarnings("ignore", category=ConvergenceWarning)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the 3x-threshold residue-distance logistic pruning analysis.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument(
        "--display-residue-csv",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/"
            "all_mutation_wt_referenced_occupancy_heatmap_wt_contacted_residues_by_region_excluding_f227c_display_residues.csv"
        ),
    )
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument("--frame-features-csv", type=Path, default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"))
    parser.add_argument(
        "--binding-ddg-summary-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/last20frames/tables/mutation_ddg_summary.csv"),
    )
    parser.add_argument(
        "--binding-ddg-full-csv",
        type=Path,
        default=Path("results/analysis/binding_energy/last20frames/tables/ddg_full.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/new_logistic_regression"))
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=5)
    parser.add_argument("--high-threshold-fold", type=float, default=3.0)
    parser.add_argument("--c-values", type=str, default="0.01,0.03,0.1,0.3,1,3,10,30")
    parser.add_argument("--penalties", type=str, default="l2,l1")
    parser.add_argument("--class-weights", type=str, default="balanced,none")
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--parallel-backend", choices=["loky", "threading"], default="loky")
    parser.add_argument("--parallel-verbose", type=int, default=0)
    parser.add_argument("--include-wt-as-negative-control", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-ligand-pose-rmsd", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--replicate-level-grouped-cv", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--within-genotype-bootstrap", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260518)
    parser.add_argument(
        "--exhaustive-max-features",
        type=int,
        default=0,
        help="Also evaluate every feature subset up to this size; 0 disables this broader search.",
    )
    parser.add_argument(
        "--reuse-distance-features",
        action="store_true",
        help="Reuse the existing replicate distance matrix in the output directory when present.",
    )
    parser.add_argument(
        "--reuse-distance-features-from",
        type=Path,
        default=None,
        help="Read the replicate distance matrix from this CSV instead of the output directory cache.",
    )
    parser.add_argument(
        "--fixed-features",
        type=str,
        default="",
        help="Comma-separated feature list to evaluate without pruning. Accepts exact columns or labels such as SER105, TYR188, ligand_pose_rmsd.",
    )
    return parser.parse_args()


def _parse_float_list(text: str) -> list[float]:
    vals = [float(token.strip()) for token in str(text).split(",") if token.strip()]
    if not vals:
        raise ValueError("Expected at least one C value.")
    return vals


def _parse_str_list(text: str) -> list[str]:
    vals = [token.strip().lower() for token in str(text).split(",") if token.strip()]
    if not vals:
        raise ValueError("Expected at least one value.")
    return vals


def _hyperparameter_grid(c_values: list[float], penalties: list[str], class_weights: list[str]) -> list[dict[str, object]]:
    grid: list[dict[str, object]] = []
    for c_value in c_values:
        for penalty in penalties:
            if penalty not in {"l1", "l2"}:
                raise ValueError(f"Unsupported penalty: {penalty}")
            for class_weight_label in class_weights:
                if class_weight_label not in {"balanced", "none"}:
                    raise ValueError(f"Unsupported class weight: {class_weight_label}")
                grid.append(
                    {
                        "c_value": float(c_value),
                        "penalty": penalty,
                        "solver": "liblinear" if penalty == "l1" else "lbfgs",
                        "class_weight": None if class_weight_label == "none" else "balanced",
                        "class_weight_label": class_weight_label,
                    }
                )
    return grid


def _model_params_from_row(row: pd.Series | dict[str, object]) -> dict[str, object]:
    class_weight_label = str(row["class_weight"])
    return {
        "c_value": float(row["c_value"]),
        "penalty": str(row["penalty"]),
        "solver": str(row["solver"]),
        "class_weight": None if class_weight_label == "none" else "balanced",
        "class_weight_label": class_weight_label,
    }


def _feature_name(label: str) -> str:
    return f"residue_min_distance_{label}_angstrom_mean"


def _resolve_feature_tokens(text: str) -> list[str]:
    features: list[str] = []
    for raw_token in str(text).split(","):
        token = raw_token.strip()
        if not token:
            continue
        normalized = token.strip().upper().replace(" ", "_").replace("-", "_")
        if normalized in {"DOR_POSE_RMSD", "LIGAND_POSE_RMSD", "POSE_RMSD"}:
            features.append("ligand_pose_rmsd_angstrom_mean")
        elif normalized in {"DDG_ELECTROSTATIC", "ELECTROSTATIC_DDG", "ELECTROSTATIC_BINDING_ENERGY", "DDG_ELEC"}:
            features.append("ddg_electrostatic_mean")
        elif normalized in {"DDG_VDW", "VDW_DDG", "VAN_DER_WAALS_DDG", "VDW_BINDING_ENERGY"}:
            features.append("ddg_vdw_mean")
        elif token.endswith("_mean") or token.startswith("residue_min_distance_"):
            features.append(token)
        else:
            features.append(_feature_name(normalized))
    return features


def _category(mutation: str) -> str:
    mutation = str(mutation).strip().upper()
    if mutation == "WT":
        return "wt_reference"
    if mutation in NEGATIVE_CONTROLS:
        return "negative_control"
    if mutation in POSITIVE_CONTROLS:
        return "positive_control"
    if mutation in UNCERTAIN_PHENOTYPE:
        return "uncertain_phenotype"
    return "other"


def _binary_class(fold_change: float, threshold: float) -> str:
    return "high" if float(fold_change) >= float(threshold) else "low"


def _pipeline(model_params: dict[str, object]) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    C=float(model_params["c_value"]),
                    penalty=str(model_params["penalty"]),
                    solver=str(model_params["solver"]),
                    class_weight=model_params["class_weight"],
                    max_iter=5000,
                    random_state=0,
                ),
            ),
        ]
    )


def _compute_replicate_distance_means(
    *,
    manifest: Path,
    residue_df: pd.DataFrame,
    ligand_resname: str,
    resid_offset: int,
    frame_stride: int,
) -> pd.DataFrame:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    needed = set(NEGATIVE_CONTROLS) | set(POSITIVE_CONTROLS) | set(UNCERTAIN_PHENOTYPE) | {"WT"}
    metas = _load_replicate_meta(manifest, needed_mutations=needed)
    residue_specs = [
        (str(row.label), int(row.auth_resid), int(row.auth_resid) + int(resid_offset))
        for row in residue_df.itertuples(index=False)
    ]

    rows: list[dict[str, object]] = []
    for meta in metas:
        u = mda.Universe(str(meta.topology_pdb), str(meta.analysis_dcd), format="DCD")
        ligand = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if ligand.n_atoms == 0:
            raise ValueError(f"No ligand atoms found for {meta.mutation} rep{meta.replicate}")
        residue_atoms = {
            label: u.select_atoms(f"protein and resid {traj_resid} and not name H*")
            for label, _auth_resid, traj_resid in residue_specs
        }
        vals = {label: [] for label, _auth_resid, _traj_resid in residue_specs}
        n_frames = 0
        for _ts in u.trajectory[:: max(1, int(frame_stride))]:
            n_frames += 1
            box = u.dimensions
            for label, _auth_resid, _traj_resid in residue_specs:
                atoms = residue_atoms[label]
                if atoms.n_atoms == 0:
                    vals[label].append(np.nan)
                else:
                    vals[label].append(float(distance_array(atoms.positions, ligand.positions, box=box).min()))
        row = {"mutation": str(meta.mutation), "replicate": int(meta.replicate), "n_frames": int(n_frames)}
        for label, auth_resid, traj_resid in residue_specs:
            row.update(
                {
                    f"residue_min_distance_{label}_angstrom_mean": float(np.nanmean(vals[label])),
                    f"residue_min_distance_{label}_auth_resid": int(auth_resid),
                    f"residue_min_distance_{label}_traj_resid": int(traj_resid),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)


def _mutation_feature_matrix(rep_df: pd.DataFrame, residue_df: pd.DataFrame) -> pd.DataFrame:
    feature_cols = [_feature_name(str(label)) for label in residue_df["label"].astype(str)]
    agg = rep_df.groupby("mutation", as_index=False)[feature_cols].mean()
    return agg


def _ligand_pose_rmsd_feature(frame_features_csv: Path) -> pd.DataFrame:
    frame_df = pd.read_csv(frame_features_csv, usecols=["mutation", "ligand_pose_rmsd_angstrom"])
    return (
        frame_df.groupby("mutation", as_index=False)["ligand_pose_rmsd_angstrom"]
        .mean()
        .rename(columns={"ligand_pose_rmsd_angstrom": "ligand_pose_rmsd_angstrom_mean"})
    )


def _ligand_pose_rmsd_replicate_feature(frame_features_csv: Path) -> pd.DataFrame:
    frame_df = pd.read_csv(frame_features_csv, usecols=["mutation", "replicate", "ligand_pose_rmsd_angstrom"])
    return (
        frame_df.groupby(["mutation", "replicate"], as_index=False)["ligand_pose_rmsd_angstrom"]
        .mean()
        .rename(columns={"ligand_pose_rmsd_angstrom": "ligand_pose_rmsd_angstrom_mean"})
    )


def _ddg_component_feature(binding_ddg_summary_csv: Path, component_col: str) -> pd.DataFrame:
    ddg_df = pd.read_csv(binding_ddg_summary_csv, usecols=["mutation", component_col])
    wt = pd.DataFrame([{"mutation": "WT", component_col: 0.0}])
    return pd.concat([ddg_df, wt], ignore_index=True).drop_duplicates("mutation", keep="last")


def _ddg_component_replicate_feature(binding_ddg_full_csv: Path, raw_col: str, feature_col: str) -> pd.DataFrame:
    ddg_df = pd.read_csv(binding_ddg_full_csv, usecols=["mutation", "replicate", raw_col])
    return ddg_df.rename(columns={raw_col: feature_col})


def _metadata_columns(panel_df: pd.DataFrame) -> pd.DataFrame:
    return panel_df[["mutation", "control_category", "target_fold_change", "target_binary_class"]].drop_duplicates(
        "mutation", keep="last"
    )


def _build_replicate_feature_matrix(
    rep_df: pd.DataFrame,
    panel_df: pd.DataFrame,
    feature_cols: list[str],
    *,
    frame_features_csv: Path,
    binding_ddg_full_csv: Path,
    include_ligand_pose_rmsd: bool,
    include_ddg_electrostatic: bool,
    include_ddg_vdw: bool,
) -> pd.DataFrame:
    base_cols = ["mutation", "replicate", *[c for c in feature_cols if c in rep_df.columns]]
    out = rep_df[base_cols].copy()
    if include_ligand_pose_rmsd:
        out = out.merge(_ligand_pose_rmsd_replicate_feature(frame_features_csv), on=["mutation", "replicate"], how="left")
    if include_ddg_electrostatic:
        out = out.merge(
            _ddg_component_replicate_feature(binding_ddg_full_csv, "ddg_electrostatic", "ddg_electrostatic_mean"),
            on=["mutation", "replicate"],
            how="left",
        )
    if include_ddg_vdw:
        out = out.merge(
            _ddg_component_replicate_feature(binding_ddg_full_csv, "ddg_vdw", "ddg_vdw_mean"),
            on=["mutation", "replicate"],
            how="left",
        )
    out = out.merge(_metadata_columns(panel_df), on="mutation", how="left")
    return out[["mutation", "replicate", "control_category", "target_fold_change", "target_binary_class", *feature_cols]].sort_values(
        ["mutation", "replicate"], kind="stable"
    )


def _metrics(y_true: np.ndarray, y_pred: np.ndarray, prob_high: np.ndarray) -> dict[str, float]:
    y_bin = (y_true == "high").astype(int)
    out = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
    }
    if len(np.unique(y_bin)) == 2:
        out["roc_auc"] = float(roc_auc_score(y_bin, prob_high))
    else:
        out["roc_auc"] = np.nan
    cm = confusion_matrix(y_true, y_pred, labels=["low", "high"])
    out.update({"tn": int(cm[0, 0]), "fp": int(cm[0, 1]), "fn": int(cm[1, 0]), "tp": int(cm[1, 1])})
    return out


def _evaluate_features(train_df: pd.DataFrame, features: list[str], hyperparams: list[dict[str, object]]) -> dict[str, object]:
    x = train_df[features].copy()
    y = train_df["target_binary_class"].astype(str).to_numpy()
    loo = LeaveOneOut()

    best: dict[str, object] | None = None
    for params in hyperparams:
        pred = np.empty(len(train_df), dtype=object)
        prob_high = np.full(len(train_df), np.nan, dtype=float)
        for train_idx, test_idx in loo.split(x, y):
            model = _pipeline(params)
            model.fit(x.iloc[train_idx], y[train_idx])
            classes = list(model.named_steps["model"].classes_)
            high_idx = int(classes.index("high"))
            prob = float(model.predict_proba(x.iloc[test_idx])[:, high_idx][0])
            prob_high[test_idx[0]] = prob
            pred[test_idx[0]] = "high" if prob >= 0.5 else "low"
        row = {
            "c_value": float(params["c_value"]),
            "penalty": str(params["penalty"]),
            "solver": str(params["solver"]),
            "class_weight": str(params["class_weight_label"]),
            "n_features": int(len(features)),
            "features": "|".join(features),
            **_metrics(y, pred.astype(str), prob_high),
        }
        if best is None or (
            row["accuracy"],
            row["balanced_accuracy"],
            row["macro_f1"],
            row["roc_auc"],
            1 if row["penalty"] == "l1" else 0,
            1 if row["class_weight"] == "balanced" else 0,
            -row["c_value"],
        ) > (
            best["accuracy"],
            best["balanced_accuracy"],
            best["macro_f1"],
            best["roc_auc"],
            1 if best["penalty"] == "l1" else 0,
            1 if best["class_weight"] == "balanced" else 0,
            -best["c_value"],
        ):
            best = row
    assert best is not None
    return best


def _heldout_probability(
    *,
    panel_df: pd.DataFrame,
    train_df: pd.DataFrame,
    features: list[str],
    model_params: dict[str, object],
    mutation: str,
) -> float:
    heldout = panel_df[panel_df["mutation"] == mutation]
    if heldout.empty:
        return np.nan
    model = _pipeline(model_params)
    model.fit(train_df[features], train_df["target_binary_class"].astype(str))
    classes = list(model.named_steps["model"].classes_)
    high_idx = int(classes.index("high"))
    return float(model.predict_proba(heldout[features])[:, high_idx][0])


def _loo_predictions(train_df: pd.DataFrame, features: list[str], model_params: dict[str, object]) -> pd.DataFrame:
    x = train_df[features].copy()
    y = train_df["target_binary_class"].astype(str).to_numpy()
    loo = LeaveOneOut()
    rows: list[dict[str, object]] = []
    for fold, (train_idx, test_idx) in enumerate(loo.split(x, y), start=1):
        model = _pipeline(model_params)
        model.fit(x.iloc[train_idx], y[train_idx])
        classes = list(model.named_steps["model"].classes_)
        high_idx = int(classes.index("high"))
        prob = float(model.predict_proba(x.iloc[test_idx])[:, high_idx][0])
        sample = train_df.iloc[test_idx[0]]
        rows.append(
            {
                "fold": int(fold),
                "mutation": str(sample["mutation"]),
                "control_category": str(sample["control_category"]),
                "target_fold_change": float(sample["target_fold_change"]),
                "observed_class": str(sample["target_binary_class"]),
                "prob_high": prob,
                "predicted_class": "high" if prob >= 0.5 else "low",
            }
        )
    return pd.DataFrame(rows)


def _evaluate_replicate_grouped_cv(
    train_rep_df: pd.DataFrame,
    features: list[str],
    hyperparams: list[dict[str, object]],
) -> tuple[dict[str, object], pd.DataFrame]:
    mutations = train_rep_df["mutation"].drop_duplicates().tolist()
    best: dict[str, object] | None = None
    best_pred_df = pd.DataFrame()

    for params in hyperparams:
        rows: list[dict[str, object]] = []
        for mutation in mutations:
            tr = train_rep_df[train_rep_df["mutation"] != mutation]
            te = train_rep_df[train_rep_df["mutation"] == mutation]
            model = _pipeline(params)
            model.fit(tr[features], tr["target_binary_class"].astype(str))
            classes = list(model.named_steps["model"].classes_)
            high_idx = int(classes.index("high"))
            rep_probs = model.predict_proba(te[features])[:, high_idx]
            prob = float(np.mean(rep_probs))
            sample = te.iloc[0]
            rows.append(
                {
                    "mutation": str(mutation),
                    "control_category": str(sample["control_category"]),
                    "target_fold_change": float(sample["target_fold_change"]),
                    "observed_class": str(sample["target_binary_class"]),
                    "prob_high": prob,
                    "predicted_class": "high" if prob >= 0.5 else "low",
                    "n_replicates": int(len(te)),
                    "replicate_prob_high_values": "|".join(f"{p:.6g}" for p in rep_probs),
                }
            )
        pred_df = pd.DataFrame(rows)
        row = {
            "c_value": float(params["c_value"]),
            "penalty": str(params["penalty"]),
            "solver": str(params["solver"]),
            "class_weight": str(params["class_weight_label"]),
            "n_features": int(len(features)),
            "features": "|".join(features),
            **_metrics(
                pred_df["observed_class"].astype(str).to_numpy(),
                pred_df["predicted_class"].astype(str).to_numpy(),
                pred_df["prob_high"].astype(float).to_numpy(),
            ),
        }
        if best is None or (
            row["accuracy"],
            row["balanced_accuracy"],
            row["macro_f1"],
            row["roc_auc"],
            1 if row["penalty"] == "l1" else 0,
            1 if row["class_weight"] == "balanced" else 0,
            -row["c_value"],
        ) > (
            best["accuracy"],
            best["balanced_accuracy"],
            best["macro_f1"],
            best["roc_auc"],
            1 if best["penalty"] == "l1" else 0,
            1 if best["class_weight"] == "balanced" else 0,
            -best["c_value"],
        ):
            best = row
            best_pred_df = pred_df
    assert best is not None
    return best, best_pred_df


def _fit_all_predictions(panel_df: pd.DataFrame, train_df: pd.DataFrame, features: list[str], model_params: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = _pipeline(model_params)
    model.fit(train_df[features], train_df["target_binary_class"].astype(str))
    classes = list(model.named_steps["model"].classes_)
    high_idx = int(classes.index("high"))
    out = panel_df[["mutation", "control_category", "target_fold_change", "target_binary_class"]].copy()
    out["prob_high"] = model.predict_proba(panel_df[features])[:, high_idx]
    out["predicted_class"] = np.where(out["prob_high"] >= 0.5, "high", "low")
    coefs = pd.DataFrame(
        {
            "feature": features,
            "coefficient": model.named_steps["model"].coef_.reshape(-1),
        }
    )
    coefs["abs_coefficient"] = coefs["coefficient"].abs()
    coefs = coefs.sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)
    return out, coefs


def _fit_replicate_level_predictions(
    panel_df: pd.DataFrame,
    rep_df: pd.DataFrame,
    train_rep_df: pd.DataFrame,
    features: list[str],
    model_params: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = _pipeline(model_params)
    model.fit(train_rep_df[features], train_rep_df["target_binary_class"].astype(str))
    classes = list(model.named_steps["model"].classes_)
    high_idx = int(classes.index("high"))
    scored = rep_df[["mutation", "replicate", *features]].copy()
    scored["replicate_prob_high"] = model.predict_proba(rep_df[features])[:, high_idx]
    out = (
        scored.groupby("mutation", as_index=False)["replicate_prob_high"]
        .mean()
        .rename(columns={"replicate_prob_high": "prob_high"})
    )
    out = _metadata_columns(panel_df).merge(out, on="mutation", how="left")
    out["predicted_class"] = np.where(out["prob_high"] >= 0.5, "high", "low")
    coefs = pd.DataFrame(
        {
            "feature": features,
            "coefficient": model.named_steps["model"].coef_.reshape(-1),
        }
    )
    coefs["abs_coefficient"] = coefs["coefficient"].abs()
    return out, coefs.sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)


def _bootstrap_probability_intervals(
    panel_df: pd.DataFrame,
    train_df: pd.DataFrame,
    features: list[str],
    model_params: dict[str, object],
    n_bootstrap: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(int(seed))
    high_idx = train_df.index[train_df["target_binary_class"].astype(str) == "high"].to_numpy()
    low_idx = train_df.index[train_df["target_binary_class"].astype(str) == "low"].to_numpy()
    if len(high_idx) == 0 or len(low_idx) == 0:
        raise ValueError("Bootstrap requires at least one low and one high training sample.")

    pred_rows: list[dict[str, object]] = []
    for iteration in range(1, int(n_bootstrap) + 1):
        boot_index = np.concatenate(
            [
                rng.choice(low_idx, size=len(low_idx), replace=True),
                rng.choice(high_idx, size=len(high_idx), replace=True),
            ]
        )
        boot_train = train_df.loc[boot_index].reset_index(drop=True)
        model = _pipeline(model_params)
        model.fit(boot_train[features], boot_train["target_binary_class"].astype(str))
        classes = list(model.named_steps["model"].classes_)
        high_class_idx = int(classes.index("high"))
        probs = model.predict_proba(panel_df[features])[:, high_class_idx]
        for mutation, prob in zip(panel_df["mutation"], probs):
            pred_rows.append({"bootstrap_iteration": iteration, "mutation": str(mutation), "prob_high": float(prob)})

    bootstrap_df = pd.DataFrame(pred_rows)
    summary = (
        bootstrap_df.groupby("mutation")["prob_high"]
        .quantile([0.025, 0.5, 0.975])
        .unstack()
        .rename(columns={0.025: "prob_high_ci_lower", 0.5: "prob_high_bootstrap_median", 0.975: "prob_high_ci_upper"})
        .reset_index()
    )
    metadata_cols = ["mutation", "control_category", "target_fold_change", "target_binary_class"]
    full_fit, _coefs = _fit_all_predictions(panel_df, train_df, features, model_params)
    summary = (
        panel_df[metadata_cols]
        .merge(full_fit[["mutation", "prob_high", "predicted_class"]], on="mutation", how="left")
        .merge(summary, on="mutation", how="left")
    )
    return bootstrap_df, summary


def _within_genotype_bootstrap_probability_intervals(
    rep_df: pd.DataFrame,
    panel_df: pd.DataFrame,
    train_categories: set[str],
    features: list[str],
    model_params: dict[str, object],
    n_bootstrap: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(int(seed))
    mutations = rep_df["mutation"].drop_duplicates().tolist()
    pred_rows: list[dict[str, object]] = []
    meta = _metadata_columns(panel_df)

    for iteration in range(1, int(n_bootstrap) + 1):
        boot_rows: list[pd.DataFrame] = []
        for mutation in mutations:
            sub = rep_df[rep_df["mutation"] == mutation]
            sampled_idx = rng.choice(sub.index.to_numpy(), size=len(sub), replace=True)
            boot_rows.append(sub.loc[sampled_idx])
        boot_rep = pd.concat(boot_rows, ignore_index=True)
        boot_panel = boot_rep.groupby("mutation", as_index=False)[features].mean().merge(meta, on="mutation", how="left")
        boot_train = boot_panel[boot_panel["control_category"].isin(train_categories)].copy()
        model = _pipeline(model_params)
        model.fit(boot_train[features], boot_train["target_binary_class"].astype(str))
        classes = list(model.named_steps["model"].classes_)
        high_idx = int(classes.index("high"))
        probs = model.predict_proba(boot_panel[features])[:, high_idx]
        for mutation, prob in zip(boot_panel["mutation"], probs):
            pred_rows.append({"bootstrap_iteration": int(iteration), "mutation": str(mutation), "prob_high": float(prob)})

    bootstrap_df = pd.DataFrame(pred_rows)
    summary = (
        bootstrap_df.groupby("mutation")["prob_high"]
        .quantile([0.025, 0.5, 0.975])
        .unstack()
        .rename(columns={0.025: "prob_high_ci_lower", 0.5: "prob_high_bootstrap_median", 0.975: "prob_high_ci_upper"})
        .reset_index()
    )
    full_panel = rep_df.groupby("mutation", as_index=False)[features].mean().merge(meta, on="mutation", how="left")
    full_train = full_panel[full_panel["control_category"].isin(train_categories)].copy()
    full_fit, _coefs = _fit_all_predictions(full_panel, full_train, features, model_params)
    return bootstrap_df, full_fit.merge(summary, on="mutation", how="left")


def _rank_models(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(
        ["accuracy", "balanced_accuracy", "macro_f1", "roc_auc", "wt_prob_high", "n_features"],
        ascending=[False, False, False, False, True, True],
        kind="stable",
    ).reset_index(drop=True)


def _evaluate_subset_search(
    *,
    panel_df: pd.DataFrame,
    train_df: pd.DataFrame,
    feature_cols: list[str],
    hyperparams: list[dict[str, object]],
    max_features: int,
    n_jobs: int,
    backend: str,
    verbose: int,
) -> pd.DataFrame:
    max_features = min(int(max_features), len(feature_cols))
    subsets = [
        combo
        for size in range(1, max_features + 1)
        for combo in combinations(feature_cols, size)
    ]

    def _eval_combo(combo: tuple[str, ...]) -> dict[str, object]:
        features = list(combo)
        row = _evaluate_features(train_df, features, hyperparams)
        row["wt_prob_high"] = _heldout_probability(
            panel_df=panel_df,
            train_df=train_df,
            features=features,
            model_params=_model_params_from_row(row),
            mutation="WT",
        )
        return row

    rows = Parallel(n_jobs=int(n_jobs), backend=str(backend), verbose=int(verbose), batch_size=32)(
        delayed(_eval_combo)(combo) for combo in subsets
    )
    return _rank_models(pd.DataFrame(rows))


def _plot_pruning_path(path_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    ax.plot(path_df["n_features"], path_df["accuracy"], marker="o", linewidth=2.4, label="Accuracy")
    ax.plot(path_df["n_features"], path_df["balanced_accuracy"], marker="s", linewidth=2.0, label="Balanced accuracy")
    ax.invert_xaxis()
    ax.set_xlabel("Number of residue-distance features", fontsize=16, fontweight="bold")
    ax.set_ylabel("LOO CV performance", fontsize=16, fontweight="bold")
    ax.tick_params(axis="both", labelsize=13)
    ax.set_ylim(0.0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linestyle=":")
    ax.legend(frameon=False, fontsize=12)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)


def _plot_probabilities(pred_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    order = pred_df.sort_values(["control_category", "target_fold_change", "mutation"], kind="stable").reset_index(drop=True)
    colors = {
        "negative_control": "#4c78a8",
        "positive_control": "#e45756",
        "uncertain_phenotype": "#9aa0a6",
        "wt_reference": "#333333",
        "other": "#cccccc",
    }
    fig, ax = plt.subplots(figsize=(13.2, 6.2))
    x = np.arange(len(order))
    ax.bar(x, order["prob_high"], color=[colors.get(c, "#cccccc") for c in order["control_category"]], edgecolor="#333333", linewidth=0.5)
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=1.0)
    ax.set_xticks(x, order["mutation"].tolist(), rotation=45, ha="right")
    ax.set_ylabel("Model-predicted probability\nof >=3x DOR fold-change", fontsize=16, fontweight="bold")
    ax.set_xlabel("")
    ax.tick_params(axis="both", labelsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.22, linestyle=":")
    fig.subplots_adjust(left=0.14, bottom=0.34, right=0.99, top=0.97)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320)
    plt.close(fig)


def _plot_probability_intervals(summary_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    order = summary_df.sort_values(["control_category", "target_fold_change", "mutation"], kind="stable").reset_index(drop=True)
    colors = {
        "negative_control": "#4c78a8",
        "positive_control": "#e45756",
        "uncertain_phenotype": "#9aa0a6",
        "wt_reference": "#333333",
        "other": "#cccccc",
    }
    fig, ax = plt.subplots(figsize=(13.2, 6.2))
    x = np.arange(len(order))
    y = order["prob_high"].to_numpy(dtype=float)
    lower = order["prob_high_ci_lower"].to_numpy(dtype=float)
    upper = order["prob_high_ci_upper"].to_numpy(dtype=float)
    yerr = np.vstack([np.maximum(0, y - lower), np.maximum(0, upper - y)])
    ax.bar(x, y, color=[colors.get(c, "#cccccc") for c in order["control_category"]], edgecolor="#333333", linewidth=0.5)
    ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#222222", elinewidth=1.1, capsize=2.8, capthick=1.1)
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=1.0)
    ax.set_xticks(x, order["mutation"].tolist(), rotation=45, ha="right")
    ax.set_ylabel("Predicted probability\n(DOR fold-change >=3x)", fontsize=14, fontweight="bold", labelpad=10)
    ax.tick_params(axis="both", labelsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.22, linestyle=":")
    fig.subplots_adjust(left=0.18, bottom=0.34, right=0.99, top=0.97)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320)
    plt.close(fig)


def _plot_probability_vs_fold_with_intervals(
    summary_df: pd.DataFrame,
    output_png: Path,
    *,
    show_intervals: bool = True,
) -> None:
    import matplotlib.pyplot as plt

    df = summary_df[np.isfinite(summary_df["target_fold_change"])].copy()
    df["x"] = np.log10(df["target_fold_change"].astype(float).clip(lower=1e-6))
    colors = {
        "negative_control": "#4c78a8",
        "positive_control": "#e45756",
        "uncertain_phenotype": "#9aa0a6",
        "wt_reference": "#333333",
        "other": "#cccccc",
    }
    labels = {
        "negative_control": "Negative control",
        "positive_control": "Positive control",
        "uncertain_phenotype": "Uncertain phenotype",
        "wt_reference": "WT",
        "other": "Other",
    }
    fig, ax = plt.subplots(figsize=(7.8, 6.4))
    for category, sub in df.groupby("control_category", sort=False):
        y = sub["prob_high"].to_numpy(dtype=float)
        plot_kwargs = {
            "fmt": "o",
            "markersize": 7.5,
            "markeredgecolor": "#222222",
            "markeredgewidth": 0.6,
            "color": colors.get(category, "#cccccc"),
            "label": labels.get(category, str(category)),
            "alpha": 0.95,
        }
        if show_intervals:
            lower = sub["prob_high_ci_lower"].to_numpy(dtype=float)
            upper = sub["prob_high_ci_upper"].to_numpy(dtype=float)
            plot_kwargs["yerr"] = np.vstack([np.maximum(0, y - lower), np.maximum(0, upper - y)])
            plot_kwargs["capsize"] = 3
            plot_kwargs["elinewidth"] = 1.1
        ax.errorbar(sub["x"], y, **plot_kwargs)
    fit_df = df[np.isfinite(df["x"]) & np.isfinite(df["prob_high"])].copy()
    if len(fit_df) >= 3:
        fit = linregress(fit_df["x"].to_numpy(dtype=float), fit_df["prob_high"].to_numpy(dtype=float))
        x_line = np.linspace(float(fit_df["x"].min()), float(fit_df["x"].max()), 100)
        y_line = fit.intercept + fit.slope * x_line
        ax.plot(x_line, y_line, color="#222222", linewidth=1.8)
        ax.text(
            0.05,
            0.95,
            f"R$^2$ = {fit.rvalue ** 2:.2f}\np = {fit.pvalue:.3g}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
        )
    label_offsets = {
        "WT": (-4, 3, "right"),
        "G190A": (-4, 4, "right"),
        "K103N+M230L": (-4, 4, "right"),
        "A98G+F227C": (-4, 5, "right"),
        "Y188L": (-5, 6, "right"),
        "V106A+P225H": (4, -5, "left"),
        "V106A+L234I": (4, -8, "left"),
        "V106I+F227C": (-4, 5, "right"),
        "K103N": (4, -7, "left"),
        "V106I": (4, 5, "left"),
        "V106M": (4, 5, "left"),
        "L100I+K103N": (4, 5, "left"),
        "K103N+P225H": (4, 5, "left"),
        "G190E": (4, 5, "left"),
        "G190S": (4, 5, "left"),
        "V106A": (4, 5, "left"),
        "Y318F": (4, 5, "left"),
        "V106A+F227L": (4, 5, "left"),
    }
    for row in df.itertuples(index=False):
        mutation = str(row.mutation)
        xoff, yoff, ha = label_offsets.get(mutation, (4, 3, "left"))
        ax.annotate(
            mutation,
            (float(row.x), float(row.prob_high)),
            xytext=(xoff, yoff),
            textcoords="offset points",
            fontsize=8,
            ha=ha,
            va="bottom",
        )
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=1.0)
    ax.axvline(np.log10(3.0), color="#777777", linestyle=":", linewidth=1.0)
    ax.set_xlabel("Fold-change", fontsize=18, fontweight="bold")
    ax.set_ylabel("Model-predicted probability\nof >=3x DOR fold-change", fontsize=18, fontweight="bold")
    ax.set_ylim(-0.04, 1.04)
    ticks = [1, 3, 10, 30, 100, 300]
    ax.set_xticks(np.log10(ticks), [str(t) for t in ticks])
    ax.set_xlim(float(df["x"].min()) - 0.1, float(df["x"].max()) + 0.32)
    ax.tick_params(axis="both", labelsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.22, linestyle=":")
    ax.legend(frameon=False, fontsize=11, loc="lower right")
    fig.subplots_adjust(left=0.16, bottom=0.14, right=0.98, top=0.98)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320)
    plt.close(fig)


def _plot_confusion_matrix_png(cm: np.ndarray, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    cm_display = cm[[1, 0], :]
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(cm_display, cmap="Blues", vmin=0)
    ax.set_xticks([0, 1], ["Low", "High"], fontsize=13)
    ax.set_yticks([0, 1], ["High", "Low"], fontsize=13)
    ax.set_xlabel("Predicted", fontsize=15, fontweight="bold")
    ax.set_ylabel("Observed", fontsize=15, fontweight="bold")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(int(cm_display[i, j])), ha="center", va="center", fontsize=18, fontweight="bold", color="#111111")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.subplots_adjust(left=0.18, bottom=0.16, right=0.90, top=0.96)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def _plot_roc_curve_png(df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    y_true = (df["target_binary_class"].astype(str) == "high").astype(int).to_numpy()
    y_score = df["prob_high"].astype(float).to_numpy()
    fpr, tpr, _thresholds = roc_curve(y_true, y_score)
    auc_val = roc_auc_score(y_true, y_score)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot(fpr, tpr, color="#333333", linewidth=2.4, label=f"AUC = {auc_val:.3f}")
    ax.plot([0, 1], [0, 1], color="#9aa0a6", linewidth=1.2, linestyle="--")
    ax.set_xlabel("False positive rate", fontsize=15, fontweight="bold")
    ax.set_ylabel("True positive rate", fontsize=15, fontweight="bold")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.tick_params(axis="both", labelsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.22, linestyle=":")
    ax.legend(frameon=False, fontsize=12, loc="lower right")
    fig.subplots_adjust(left=0.16, bottom=0.15, right=0.98, top=0.98)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320)
    plt.close(fig)


def _plot_coefficients_png(coef_df: pd.DataFrame, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    label_map = {
        "residue_min_distance_SER105_angstrom_mean": "SER105-DOR distance",
        "residue_min_distance_TYR188_angstrom_mean": "TYR188-DOR distance",
        "ligand_pose_rmsd_angstrom_mean": "DOR pose RMSD",
    }
    df = coef_df.copy().sort_values("coefficient", ascending=True, kind="stable")
    labels = [label_map.get(f, str(f).replace("residue_min_distance_", "").replace("_angstrom_mean", "")) for f in df["feature"]]
    colors = np.where(df["coefficient"].to_numpy(dtype=float) >= 0, "#e45756", "#4c78a8")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.barh(np.arange(len(df)), df["coefficient"], color=colors, edgecolor="#222222", linewidth=0.5)
    ax.axvline(0, color="#333333", linewidth=1.0)
    ax.set_yticks(np.arange(len(df)), labels, fontsize=13)
    ax.set_xlabel("Standardized logistic coefficient", fontsize=15, fontweight="bold")
    ax.tick_params(axis="x", labelsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.22, linestyle=":")
    fig.subplots_adjust(left=0.34, bottom=0.18, right=0.98, top=0.96)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320)
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    c_values = _parse_float_list(args.c_values)
    penalties = _parse_str_list(args.penalties)
    class_weights = _parse_str_list(args.class_weights)
    hyperparams = _hyperparameter_grid(c_values, penalties, class_weights)
    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    for d in (out_tables, out_plots, out_config):
        d.mkdir(parents=True, exist_ok=True)

    residue_df = pd.read_csv(args.display_residue_csv)
    feature_cols = [_feature_name(str(label)) for label in residue_df["label"].astype(str)]
    fixed_features = _resolve_feature_tokens(args.fixed_features)
    include_ligand_pose_rmsd = bool(args.include_ligand_pose_rmsd) or "ligand_pose_rmsd_angstrom_mean" in fixed_features
    include_ddg_electrostatic = "ddg_electrostatic_mean" in fixed_features
    include_ddg_vdw = "ddg_vdw_mean" in fixed_features
    if include_ligand_pose_rmsd:
        feature_cols.append("ligand_pose_rmsd_angstrom_mean")
    if include_ddg_electrostatic:
        feature_cols.append("ddg_electrostatic_mean")
    if include_ddg_vdw:
        feature_cols.append("ddg_vdw_mean")
    if fixed_features:
        missing_fixed = sorted(set(fixed_features) - set(feature_cols))
        if missing_fixed:
            raise ValueError(f"Fixed features are not available: {missing_fixed}")
        feature_cols = fixed_features

    replicate_matrix_path = out_tables / "replicate_residue_distance_feature_matrix.csv"
    if args.reuse_distance_features_from is not None:
        rep_df = pd.read_csv(args.reuse_distance_features_from)
    elif args.reuse_distance_features and replicate_matrix_path.exists():
        rep_df = pd.read_csv(replicate_matrix_path)
    else:
        rep_df = _compute_replicate_distance_means(
            manifest=args.manifest,
            residue_df=residue_df,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            frame_stride=int(args.frame_stride),
        )
    rep_df.to_csv(out_tables / "replicate_residue_distance_feature_matrix.csv", index=False)

    panel_df = _mutation_feature_matrix(rep_df, residue_df)
    if include_ligand_pose_rmsd:
        panel_df = panel_df.merge(_ligand_pose_rmsd_feature(args.frame_features_csv), on="mutation", how="left")
    if include_ddg_electrostatic:
        panel_df = panel_df.merge(
            _ddg_component_feature(args.binding_ddg_summary_csv, "ddg_electrostatic_mean"),
            on="mutation",
            how="left",
        )
    if include_ddg_vdw:
        panel_df = panel_df.merge(
            _ddg_component_feature(args.binding_ddg_summary_csv, "ddg_vdw_mean"),
            on="mutation",
            how="left",
        )
    fold_df = load_dor_susceptibilities(args.susceptibility_xlsx)[["mutation", "dor_fold_reduction"]].rename(
        columns={"dor_fold_reduction": "target_fold_change"}
    )
    fold_df = pd.concat([fold_df, pd.DataFrame([{"mutation": "WT", "target_fold_change": 1.0}])], ignore_index=True)
    panel_df = panel_df.merge(fold_df.drop_duplicates("mutation", keep="last"), on="mutation", how="left")
    panel_df["control_category"] = panel_df["mutation"].map(_category)
    panel_df["target_binary_class"] = panel_df["target_fold_change"].map(
        lambda x: _binary_class(float(x), float(args.high_threshold_fold)) if np.isfinite(x) else np.nan
    )
    panel_df = panel_df[["mutation", "control_category", "target_fold_change", "target_binary_class", *feature_cols]].sort_values(
        ["control_category", "target_fold_change", "mutation"],
        kind="stable",
    )
    panel_df.to_csv(out_tables / "mutation_residue_distance_feature_matrix.csv", index=False)

    train_categories = {"negative_control", "positive_control"}
    if bool(args.include_wt_as_negative_control):
        train_categories.add("wt_reference")
    train_df = panel_df[panel_df["control_category"].isin(train_categories)].copy().reset_index(drop=True)
    if set(train_df["target_binary_class"]) != {"low", "high"}:
        raise ValueError("Training controls do not contain both high and low classes.")

    rep_panel_df = pd.DataFrame()
    train_rep_df = pd.DataFrame()
    grouped_best_model: dict[str, object] | None = None
    grouped_best_params: dict[str, object] | None = None
    grouped_cv_df = pd.DataFrame()
    grouped_all_pred_df = pd.DataFrame()
    grouped_coef_df = pd.DataFrame()
    within_bootstrap_df = pd.DataFrame()
    within_bootstrap_summary_df = pd.DataFrame()
    if bool(args.replicate_level_grouped_cv) or bool(args.within_genotype_bootstrap):
        rep_panel_df = _build_replicate_feature_matrix(
            rep_df=rep_df,
            panel_df=panel_df,
            feature_cols=feature_cols,
            frame_features_csv=args.frame_features_csv,
            binding_ddg_full_csv=args.binding_ddg_full_csv,
            include_ligand_pose_rmsd=include_ligand_pose_rmsd,
            include_ddg_electrostatic=include_ddg_electrostatic,
            include_ddg_vdw=include_ddg_vdw,
        )
        train_rep_df = rep_panel_df[rep_panel_df["control_category"].isin(train_categories)].copy().reset_index(drop=True)

    active_features = list(feature_cols)
    all_round_rows: list[dict[str, object]] = []
    path_rows: list[dict[str, object]] = []
    round_idx = 0
    current_eval = _evaluate_features(train_df, active_features, hyperparams)
    current_eval["wt_prob_high"] = _heldout_probability(
        panel_df=panel_df,
        train_df=train_df,
        features=active_features,
        model_params=_model_params_from_row(current_eval),
        mutation="WT",
    )
    current_eval.update({"round": round_idx, "dropped_feature": "", "feature_set_source": "initial_all_features"})
    path_rows.append(current_eval)

    while len(active_features) > 1 and not fixed_features:
        round_idx += 1

        def _drop_eval(drop_feature: str) -> dict[str, object]:
            candidate = [feature for feature in active_features if feature != drop_feature]
            row = _evaluate_features(train_df, candidate, hyperparams)
            row["wt_prob_high"] = _heldout_probability(
                panel_df=panel_df,
                train_df=train_df,
                features=candidate,
                model_params=_model_params_from_row(row),
                mutation="WT",
            )
            row.update({"round": round_idx, "dropped_feature": drop_feature, "feature_set_source": "candidate_drop"})
            return row

        candidates = Parallel(
            n_jobs=int(args.n_jobs),
            backend=str(args.parallel_backend),
            verbose=int(args.parallel_verbose),
        )(
            delayed(_drop_eval)(feature) for feature in active_features
        )
        all_round_rows.extend(candidates)
        best_candidate = sorted(
            candidates,
            key=lambda row: (
                row["accuracy"],
                row["balanced_accuracy"],
                row["macro_f1"],
                row["roc_auc"],
                -row["wt_prob_high"],
                -row["n_features"],
                row["dropped_feature"],
            ),
            reverse=True,
        )[0]
        path_rows.append({**best_candidate, "feature_set_source": "selected_drop"})
        active_features = str(best_candidate["features"]).split("|")

    path_df = pd.DataFrame(path_rows)
    all_round_df = pd.DataFrame(all_round_rows)
    best_model = path_df.sort_values(
        ["accuracy", "balanced_accuracy", "macro_f1", "roc_auc", "wt_prob_high", "n_features"],
        ascending=[False, False, False, False, True, True],
        kind="stable",
    ).iloc[0]
    best_features = str(best_model["features"]).split("|")
    best_params = _model_params_from_row(best_model)

    loo_df = _loo_predictions(train_df, best_features, best_params)
    all_pred_df, coef_df = _fit_all_predictions(panel_df.reset_index(drop=True), train_df, best_features, best_params)
    bootstrap_df, bootstrap_summary_df = _bootstrap_probability_intervals(
        panel_df.reset_index(drop=True),
        train_df,
        best_features,
        best_params,
        int(args.bootstrap_iterations),
        int(args.bootstrap_seed),
    )

    if bool(args.replicate_level_grouped_cv):
        grouped_best_model, grouped_cv_df = _evaluate_replicate_grouped_cv(train_rep_df, best_features, hyperparams)
        grouped_best_params = _model_params_from_row(grouped_best_model)
        grouped_all_pred_df, grouped_coef_df = _fit_replicate_level_predictions(
            panel_df.reset_index(drop=True),
            rep_panel_df.reset_index(drop=True),
            train_rep_df,
            best_features,
            grouped_best_params,
        )

    if bool(args.within_genotype_bootstrap):
        bootstrap_params = grouped_best_params if grouped_best_params is not None else best_params
        within_bootstrap_df, within_bootstrap_summary_df = _within_genotype_bootstrap_probability_intervals(
            rep_panel_df.reset_index(drop=True),
            panel_df.reset_index(drop=True),
            train_categories,
            best_features,
            bootstrap_params,
            int(args.bootstrap_iterations),
            int(args.bootstrap_seed),
        )

    exhaustive_df = pd.DataFrame()
    exhaustive_best = None
    exhaustive_pred_df = pd.DataFrame()
    exhaustive_coef_df = pd.DataFrame()
    if int(args.exhaustive_max_features) > 0:
        exhaustive_df = _evaluate_subset_search(
            panel_df=panel_df,
            train_df=train_df,
            feature_cols=feature_cols,
            hyperparams=hyperparams,
            max_features=int(args.exhaustive_max_features),
            n_jobs=int(args.n_jobs),
            backend=str(args.parallel_backend),
            verbose=int(args.parallel_verbose),
        )
        exhaustive_best = exhaustive_df.iloc[0]
        exhaustive_best_features = str(exhaustive_best["features"]).split("|")
        exhaustive_pred_df, exhaustive_coef_df = _fit_all_predictions(
            panel_df.reset_index(drop=True),
            train_df,
            exhaustive_best_features,
            _model_params_from_row(exhaustive_best),
        )

    path_df.to_csv(out_tables / "greedy_pruning_path.csv", index=False)
    all_round_df.to_csv(out_tables / "greedy_pruning_all_candidate_drops.csv", index=False)
    pd.DataFrame([best_model.to_dict()]).to_csv(out_tables / "best_pruned_model_summary.csv", index=False)
    loo_df.to_csv(out_tables / "best_pruned_model_loo_predictions.csv", index=False)
    all_pred_df.to_csv(out_tables / "best_pruned_model_all_mutation_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "best_pruned_model_coefficients.csv", index=False)
    bootstrap_df.to_csv(out_tables / "best_pruned_model_bootstrap_probabilities.csv", index=False)
    bootstrap_summary_df.to_csv(out_tables / "best_pruned_model_bootstrap_probability_intervals.csv", index=False)
    all_cm = confusion_matrix(all_pred_df["target_binary_class"], all_pred_df["predicted_class"], labels=["low", "high"])
    pd.DataFrame(all_cm, index=["observed_low", "observed_high"], columns=["predicted_low", "predicted_high"]).to_csv(
        out_tables / "all_mutation_full_fit_confusion_matrix.csv"
    )
    if not rep_panel_df.empty:
        rep_panel_df.to_csv(out_tables / "replicate_level_feature_matrix.csv", index=False)
    if grouped_best_model is not None:
        pd.DataFrame([grouped_best_model]).to_csv(out_tables / "replicate_level_grouped_cv_model_summary.csv", index=False)
        grouped_cv_df.to_csv(out_tables / "replicate_level_grouped_cv_predictions.csv", index=False)
        grouped_all_pred_df.to_csv(out_tables / "replicate_level_full_fit_all_mutation_predictions.csv", index=False)
        grouped_coef_df.to_csv(out_tables / "replicate_level_full_fit_coefficients.csv", index=False)
    if not within_bootstrap_df.empty:
        within_bootstrap_df.to_csv(out_tables / "within_genotype_bootstrap_probabilities.csv", index=False)
        within_bootstrap_summary_df.to_csv(out_tables / "within_genotype_bootstrap_probability_intervals.csv", index=False)
        within_bootstrap_cm = confusion_matrix(
            within_bootstrap_summary_df["target_binary_class"],
            within_bootstrap_summary_df["predicted_class"],
            labels=["low", "high"],
        )
        pd.DataFrame(
            within_bootstrap_cm,
            index=["observed_low", "observed_high"],
            columns=["predicted_low", "predicted_high"],
        ).to_csv(out_tables / "within_genotype_bootstrap_confusion_matrix.csv")
    if not exhaustive_df.empty:
        exhaustive_df.to_csv(out_tables / f"exhaustive_subset_search_up_to_{int(args.exhaustive_max_features)}_features.csv", index=False)
        pd.DataFrame([exhaustive_best.to_dict()]).to_csv(out_tables / "best_wt_sanity_model_summary.csv", index=False)
        exhaustive_pred_df.to_csv(out_tables / "best_wt_sanity_model_all_mutation_predictions.csv", index=False)
        exhaustive_coef_df.to_csv(out_tables / "best_wt_sanity_model_coefficients.csv", index=False)

    _plot_pruning_path(path_df, out_plots / "greedy_pruning_path.png")
    _plot_probabilities(all_pred_df, out_plots / "best_pruned_model_probabilities.png")
    _plot_probability_intervals(bootstrap_summary_df, out_plots / "best_pruned_model_probability_intervals.png")
    _plot_probability_vs_fold_with_intervals(bootstrap_summary_df, out_plots / "all_mutation_probability_vs_fold_change.png")
    _plot_confusion_matrix_png(all_cm, out_plots / "all_mutation_full_fit_confusion_matrix.png")
    _plot_roc_curve_png(all_pred_df, out_plots / "all_mutation_full_fit_roc_curve.png")
    _plot_coefficients_png(coef_df, out_plots / "full_model_feature_coefficients.png")
    if not grouped_all_pred_df.empty:
        grouped_cm = confusion_matrix(grouped_all_pred_df["target_binary_class"], grouped_all_pred_df["predicted_class"], labels=["low", "high"])
        _plot_confusion_matrix_png(grouped_cm, out_plots / "replicate_level_full_fit_confusion_matrix.png")
        _plot_roc_curve_png(grouped_all_pred_df, out_plots / "replicate_level_full_fit_roc_curve.png")
        _plot_coefficients_png(grouped_coef_df, out_plots / "replicate_level_full_fit_coefficients.png")
    if not within_bootstrap_summary_df.empty:
        _plot_probability_vs_fold_with_intervals(
            within_bootstrap_summary_df,
            out_plots / "within_genotype_bootstrap_probability_vs_fold_change.png",
            show_intervals=False,
        )
        _plot_confusion_matrix_png(within_bootstrap_cm, out_plots / "replicate_level_full_fit_confusion_matrix.png")
        _plot_confusion_matrix_png(within_bootstrap_cm, out_plots / "within_genotype_bootstrap_confusion_matrix.png")
        _plot_probability_intervals(
            within_bootstrap_summary_df,
            out_plots / "within_genotype_bootstrap_probability_intervals.png",
        )
    if not exhaustive_pred_df.empty:
        _plot_probabilities(exhaustive_pred_df, out_plots / "best_wt_sanity_model_probabilities.png")

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "analysis": "new_logistic_regression",
                "design": "greedy_backward_pruning_from_27_wt_contacted_residue_dor_minimum_heavy_atom_distances",
                "manifest": str(args.manifest),
                "display_residue_csv": str(args.display_residue_csv),
                "frame_features_csv": str(args.frame_features_csv),
                "binding_ddg_summary_csv": str(args.binding_ddg_summary_csv),
                "binding_ddg_full_csv": str(args.binding_ddg_full_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "ligand_resname": str(args.ligand_resname),
                "resid_offset": int(args.resid_offset),
                "frame_stride": int(args.frame_stride),
                "high_threshold_fold": float(args.high_threshold_fold),
                "c_values": c_values,
                "penalties": penalties,
                "class_weights": class_weights,
                "n_hyperparameter_combinations": int(len(hyperparams)),
                "n_jobs": int(args.n_jobs),
                "parallel_backend": str(args.parallel_backend),
                "parallel_verbose": int(args.parallel_verbose),
                "exhaustive_max_features": int(args.exhaustive_max_features),
                "reuse_distance_features": bool(args.reuse_distance_features),
                "reuse_distance_features_from": str(args.reuse_distance_features_from) if args.reuse_distance_features_from else None,
                "fixed_features": fixed_features,
                "n_features_start": int(len(feature_cols)),
                "feature_columns": feature_cols,
                "n_training_controls": int(len(train_df)),
                "train_categories": sorted(train_categories),
                "include_wt_as_negative_control": bool(args.include_wt_as_negative_control),
                "include_ligand_pose_rmsd": bool(include_ligand_pose_rmsd),
                "include_ddg_electrostatic": bool(include_ddg_electrostatic),
                "include_ddg_vdw": bool(include_ddg_vdw),
                "replicate_level_grouped_cv": bool(args.replicate_level_grouped_cv),
                "within_genotype_bootstrap": bool(args.within_genotype_bootstrap),
                "bootstrap_iterations": int(args.bootstrap_iterations),
                "bootstrap_seed": int(args.bootstrap_seed),
                "negative_controls": sorted(NEGATIVE_CONTROLS),
                "positive_controls": sorted(POSITIVE_CONTROLS),
                "uncertain_phenotype": sorted(UNCERTAIN_PHENOTYPE),
                "best_model": best_model.to_dict(),
                "replicate_level_grouped_cv_best_model": grouped_best_model,
            },
            indent=2,
        )
    )
    print(f"Saved {out_tables / 'best_pruned_model_summary.csv'}")
    print(f"Saved {out_plots / 'greedy_pruning_path.png'}")
    print(f"Saved {out_plots / 'best_pruned_model_probabilities.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
