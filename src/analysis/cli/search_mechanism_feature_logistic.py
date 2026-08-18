#!/usr/bin/env python3
"""Search small logistic models over a hand-picked mechanism-driven feature panel."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

from .model_susceptibility_from_state_features import _mutation_feature_matrix
from .plot_best_combo_logistic import main as plot_best_combo_logistic_main
from .search_feature_combo_logistic import (
    META_COLUMNS,
    _evaluate_combo,
    _feature_target_correlations,
    _parse_float_list,
    _parse_int_list,
)
from ..susceptibility import load_dor_susceptibilities


MECHANISM_FEATURES = [
    "binding_dg_electrostatic_mean",
    "ligand_pose_rmsd_angstrom_mean",
    "ligand_palm_distance_angstrom_mean",
    "residue_min_distance_LYS103_angstrom_mean",
    "residue_min_distance_VAL106_angstrom_mean",
    "residue_min_distance_VAL179_angstrom_mean",
    "residue_min_distance_TYR188_angstrom_mean",
    "residue_min_distance_VAL189_angstrom_mean",
    "residue_min_distance_PHE227_angstrom_mean",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Mechanism-driven logistic combo search for DOR susceptibility.")
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
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/mechanism_feature_logistic"),
    )
    parser.add_argument("--target-col", type=str, default="target_fold_reduction")
    parser.add_argument("--combo-sizes", type=str, default="1,2,3,4")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--low-max-fold", type=float, default=10.0)
    parser.add_argument("--inner-scoring", type=str, default="balanced_accuracy")
    parser.add_argument("--l2-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--l1-c-values", type=str, default="0.01,0.1,1,10")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--parallel-backend", type=str, default="loky", choices=["threading", "loky"])
    parser.add_argument("--top-n-diagnostics", type=int, default=4)
    args = parser.parse_args()

    if not args.frame_feature_csv.exists():
        raise FileNotFoundError(args.frame_feature_csv)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    combo_sizes = _parse_int_list(args.combo_sizes)
    l2_c_values = _parse_float_list(args.l2_c_values)
    l1_c_values = _parse_float_list(args.l1_c_values)

    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    frame_df = pd.read_csv(args.frame_feature_csv)
    mmgbsa_df = pd.read_csv(args.mmgbsa_replicate_csv) if args.mmgbsa_replicate_csv.exists() else None
    target_df = load_dor_susceptibilities(args.susceptibility_xlsx)
    feat_all = _mutation_feature_matrix(
        frame_df,
        target_df=target_df,
        temperature_k=300.0,
        dispersion_mode="replicate_sd",
        mmgbsa_df=mmgbsa_df,
    )

    missing = [feature for feature in MECHANISM_FEATURES if feature not in feat_all.columns]
    if missing:
        raise ValueError(f"Missing mechanism features: {missing}")

    feat = feat_all[[column for column in feat_all.columns if column in META_COLUMNS] + MECHANISM_FEATURES].copy()
    feat = feat.sort_values(["target_fold_reduction", "mutation"], ascending=[True, True], kind="stable").reset_index(drop=True)
    feat.to_csv(out_tables / "mutation_feature_matrix.csv", index=False)

    corr_df = _feature_target_correlations(
        feat,
        feature_cols=list(MECHANISM_FEATURES),
        target_col=str(args.target_col),
    )
    filtered_df = corr_df.copy().reset_index(drop=True)
    filtered_df.to_csv(out_tables / "filtered_features.csv", index=False)
    corr_df.to_csv(out_tables / "feature_target_correlations.csv", index=False)

    tasks: list[tuple[tuple[str, ...], int]] = []
    for combo_size in sorted(set(int(size) for size in combo_sizes)):
        if combo_size <= 0 or combo_size > len(MECHANISM_FEATURES):
            continue
        for combo in itertools.combinations(MECHANISM_FEATURES, combo_size):
            tasks.append((combo, combo_size))

    penalty_c_map = {
        "l2": l2_c_values,
        "l1": l1_c_values,
    }

    parallel_kwargs = {"n_jobs": int(args.n_jobs), "verbose": 10}
    if str(args.parallel_backend) == "threading":
        parallel_kwargs["prefer"] = "threads"
    results = Parallel(**parallel_kwargs)(
        delayed(_evaluate_combo)(
            feat,
            combo=combo,
            combo_size=int(combo_size),
            feature_count_filtered=int(len(MECHANISM_FEATURES)),
            target_col=str(args.target_col),
            low_max_fold=float(args.low_max_fold),
            cv_folds=int(args.cv_folds),
            random_state=int(args.random_state),
            inner_scoring=str(args.inner_scoring),
            penalty_c_map={key: list(values) for key, values in penalty_c_map.items()},
        )
        for combo, combo_size in tasks
    )

    summary_rows = [row for combo_summary, _pred, _coef, _cm in results for row in combo_summary]
    pred_rows = [row for _summary, combo_pred, _coef, _cm in results for row in combo_pred]
    coef_rows = [row for _summary, _pred, combo_coef, _cm in results for row in combo_coef]
    cm_rows = [row for _summary, _pred, _coef, combo_cm in results for row in combo_cm]

    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["balanced_accuracy", "macro_f1", "roc_auc", "average_precision", "accuracy"],
        ascending=[False, False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)
    pred_df = pd.DataFrame(pred_rows)
    coef_df = pd.DataFrame(coef_rows).sort_values(
        ["penalty", "combo_size", "feature_combo", "abs_coefficient"],
        ascending=[True, True, True, False],
        kind="stable",
    ).reset_index(drop=True)
    cm_df = pd.DataFrame(cm_rows)

    summary_df.to_csv(out_tables / "combo_model_summary.csv", index=False)
    pred_df.to_csv(out_tables / "combo_cv_predictions.csv", index=False)
    coef_df.to_csv(out_tables / "combo_fullfit_coefficients.csv", index=False)
    cm_df.to_csv(out_tables / "combo_confusion_matrices.csv", index=False)
    summary_df.head(25).to_csv(out_tables / "top_combo_models.csv", index=False)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "frame_feature_csv": str(args.frame_feature_csv),
                "mmgbsa_replicate_csv": str(args.mmgbsa_replicate_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "target_col": str(args.target_col),
                "low_max_fold": float(args.low_max_fold),
                "cv_folds": int(args.cv_folds),
                "random_state": int(args.random_state),
                "inner_scoring": str(args.inner_scoring),
                "combo_sizes": combo_sizes,
                "l2_c_values": l2_c_values,
                "l1_c_values": l1_c_values,
                "mechanism_features": list(MECHANISM_FEATURES),
                "feature_design": "hand-picked mechanism-driven mean features only",
            },
            indent=2,
        )
    )

    plot_best_combo_logistic_main_args = [
        "--input-dir",
        str(args.output_dir),
        "--output-dir",
        str(args.output_dir / "diagnostics"),
        "--top-n",
        str(int(args.top_n_diagnostics)),
    ]
    import sys

    argv_prev = sys.argv[:]
    try:
        sys.argv = ["plot_best_combo_logistic.py", *plot_best_combo_logistic_main_args]
        plot_best_combo_logistic_main()
    finally:
        sys.argv = argv_prev
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
