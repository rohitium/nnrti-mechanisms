#!/usr/bin/env python3
"""Incremental analysis with checkpointing to prevent data loss on crashes.

This script separates data collection, MM/GBSA computation, and plotting into
independent steps that can be run separately or resumed from checkpoints.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd


def checkpoint_exists(checkpoint_path: Path) -> bool:
    """Check if a checkpoint file exists and is non-empty."""
    return checkpoint_path.exists() and checkpoint_path.stat().st_size > 0


def _dataset_score(df: pd.DataFrame) -> tuple[int, int, int]:
    """Score a dataset by completeness for downstream plotting."""
    if df.empty:
        return (0, 0, 0)
    if "mutation" not in df.columns:
        return (0, 0, len(df))
    mutations = pd.Series(df["mutation"]).dropna().astype(str)
    has_wt = int((mutations == "WT").any())
    n_mutations = int(mutations.nunique())
    return (has_wt, n_mutations, len(df))


def _load_best_dataset(
    candidates: list[Path],
    *,
    required_columns: set[str] | None = None,
) -> tuple[pd.DataFrame, Path | None]:
    """Load the most complete available dataset from candidate CSV files."""
    best_df = pd.DataFrame()
    best_path: Path | None = None
    best_score = (0, 0, 0)

    for path in candidates:
        if not checkpoint_exists(path):
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if required_columns and not required_columns.issubset(set(df.columns)):
            continue
        score = _dataset_score(df)
        if score > best_score:
            best_df = df
            best_path = path
            best_score = score

    return best_df, best_path


def _with_fold_change_alias(df: pd.DataFrame) -> pd.DataFrame:
    """Add a fold_change output alias for backward-compatible FR->FC transition."""
    if "fold_reduction" in df.columns and "fold_change" not in df.columns:
        out = df.copy()
        insert_at = list(out.columns).index("fold_reduction") + 1
        out.insert(insert_at, "fold_change", out["fold_reduction"])
        return out
    return df


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Incremental analysis with checkpointing")
    parser.add_argument("--step", choices=["collect", "mmgbsa", "metrics", "plots", "all"], default="all")
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--force", action="store_true", help="Force recomputation even if checkpoint exists")
    parser.add_argument("--mmgbsa-snapshots", type=int, default=100)
    parser.add_argument("--mmgbsa-discard-fraction", type=float, default=0.25)
    parser.add_argument("--metric-frame-stride", type=int, default=5)
    parser.add_argument("--metric-max-frames", type=int, default=200)
    parser.add_argument(
        "--profile-workers",
        type=int,
        default=4,
        help="Worker threads for trajectory profile collection (RMSD/COM/pocket volume)",
    )
    args = parser.parse_args()

    results_dir = args.results_dir
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    ckpt_dir = results_dir / ".checkpoints"
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    # Checkpoint paths
    ckpt_md_metadata = ckpt_dir / ".checkpoint_md_metadata.csv"
    ckpt_mmgbsa = ckpt_dir / ".checkpoint_mmgbsa_replicate_metrics.csv"
    ckpt_rmsd = ckpt_dir / ".checkpoint_rmsd_ca_profiles.csv"
    ckpt_com = ckpt_dir / ".checkpoint_com_distance_profiles.csv"
    ckpt_pocket = ckpt_dir / ".checkpoint_pocket_volume_profiles.csv"
    ckpt_boundness = ckpt_dir / ".checkpoint_boundness_qc.csv"
    ckpt_structural = ckpt_dir / ".checkpoint_structural_metrics.csv"

    # === STEP 1: Collect MD metadata ===
    if args.step in ("collect", "all"):
        if args.force or not checkpoint_exists(ckpt_md_metadata):
            logging.info("STEP 1: Collecting MD run metadata")
            from src.analysis.result_collector import collect_md_results

            md_df = collect_md_results(args.manifest, args.results_dir)
            md_df.to_csv(ckpt_md_metadata, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_md_metadata}")
        else:
            logging.info(f"STEP 1: Using cached MD metadata from {ckpt_md_metadata}")
            md_df = pd.read_csv(ckpt_md_metadata)

    # === STEP 2: MM/GBSA (incremental, per-replicate) ===
    if args.step in ("mmgbsa", "all"):
        if not checkpoint_exists(ckpt_md_metadata):
            logging.error("MD metadata checkpoint not found. Run --step collect first.")
            return 1

        md_df = pd.read_csv(ckpt_md_metadata)

        if args.force or not checkpoint_exists(ckpt_mmgbsa):
            logging.info("STEP 2: Computing MM/GBSA metrics (incremental)")
            from src.analysis.result_collector import compute_mmgbsa_metrics

            mmgbsa_df = compute_mmgbsa_metrics(
                md_df,
                ligand_resname="2KW",
                n_snapshots=args.mmgbsa_snapshots,
                discard_fraction=args.mmgbsa_discard_fraction,
            )
            mmgbsa_df.to_csv(ckpt_mmgbsa, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_mmgbsa}")
        else:
            logging.info(f"STEP 2: Using cached MM/GBSA metrics from {ckpt_mmgbsa}")

    # === STEP 3: Compute structural metrics ===
    if args.step in ("metrics", "all"):
        if not checkpoint_exists(ckpt_md_metadata):
            logging.error("MD metadata checkpoint not found. Run --step collect first.")
            return 1

        md_df = pd.read_csv(ckpt_md_metadata)
        logging.info(f"STEP 3: profile workers={max(1, args.profile_workers)}")

        # RMSD
        if args.force or not checkpoint_exists(ckpt_rmsd):
            logging.info("STEP 3a: Computing Cα RMSD profiles")
            from src.analysis.result_collector import collect_ca_rmsd_profiles

            rmsd_df = collect_ca_rmsd_profiles(
                md_df,
                frame_stride=args.metric_frame_stride,
                max_frames=args.metric_max_frames,
                workers=args.profile_workers,
            )
            rmsd_df.to_csv(ckpt_rmsd, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_rmsd}")
        else:
            logging.info(f"STEP 3a: Using cached RMSD from {ckpt_rmsd}")

        # COM distance
        if args.force or not checkpoint_exists(ckpt_com):
            logging.info("STEP 3b: Computing COM distance profiles")
            from src.analysis.result_collector import collect_com_distance_profiles

            com_df = collect_com_distance_profiles(
                md_df,
                ligand_resname="2KW",
                frame_stride=args.metric_frame_stride,
                max_frames=args.metric_max_frames,
                workers=args.profile_workers,
            )
            com_df.to_csv(ckpt_com, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_com}")
        else:
            logging.info(f"STEP 3b: Using cached COM distance from {ckpt_com}")

        # Pocket volume time series
        if args.force or not checkpoint_exists(ckpt_pocket):
            logging.info("STEP 3c: Computing pocket-volume profiles")
            from src.analysis.result_collector import collect_pocket_volume_profiles

            pocket_df = collect_pocket_volume_profiles(
                md_df,
                ligand_resname="2KW",
                frame_stride=args.metric_frame_stride,
                max_frames=args.metric_max_frames,
                workers=args.profile_workers,
            )
            pocket_df.to_csv(ckpt_pocket, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_pocket}")
        else:
            logging.info(f"STEP 3c: Using cached pocket-volume profiles from {ckpt_pocket}")

        # Boundness QC
        if args.force or not checkpoint_exists(ckpt_boundness):
            logging.info("STEP 3d: Computing boundness QC")
            from src.analysis.result_collector import compute_boundness_qc

            boundness_df = compute_boundness_qc(md_df, ligand_resname="2KW")
            boundness_df.to_csv(ckpt_boundness, index=False)
            logging.info(f"  Saved checkpoint: {ckpt_boundness}")
        else:
            logging.info(f"STEP 3d: Using cached boundness QC from {ckpt_boundness}")

        # Structural metrics (use parallel script for speed)
        if args.force or not checkpoint_exists(ckpt_structural):
            logging.info("STEP 3e: Computing structural metrics (parallel)")
            import subprocess
            force_flag = ["--force"] if args.force else []
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "src.analysis.cli.compute_structural_metrics_parallel",
                    "--manifest", str(args.manifest),
                    "--results-dir", str(args.results_dir),
                    "--output", str(ckpt_structural),
                    "--frame-stride", str(args.metric_frame_stride),
                    "--max-frames", str(args.metric_max_frames),
                ] + force_flag,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                logging.error(f"Structural metrics computation failed: {result.stderr}")
                return 1
            logging.info(f"  Saved checkpoint: {ckpt_structural}")
        else:
            logging.info(f"STEP 3e: Using cached structural metrics from {ckpt_structural}")

    # === STEP 4: Generate plots (individual, independent) ===
    if args.step in ("plots", "all"):
        logging.info("STEP 4: Generating plots")

        # Load all checkpointed data
        logging.info("  Loading checkpointed data...")
        try:
            md_df = pd.read_csv(ckpt_md_metadata) if ckpt_md_metadata.exists() else pd.DataFrame()
            mmgbsa_df, mmgbsa_source = _load_best_dataset(
                [
                    results_dir / "mmgbsa_replicate_metrics.csv",
                    ckpt_mmgbsa,
                ],
                required_columns={"mutation", "replicate", "binding_dg"},
            )
            if mmgbsa_source is not None:
                logging.info(f"  Using MM/GBSA source: {mmgbsa_source}")
            else:
                logging.warning("  No usable MM/GBSA source found.")
            rmsd_df = pd.read_csv(ckpt_rmsd) if ckpt_rmsd.exists() else pd.DataFrame()
            com_df = pd.read_csv(ckpt_com) if ckpt_com.exists() else pd.DataFrame()
            pocket_df = pd.read_csv(ckpt_pocket) if ckpt_pocket.exists() else pd.DataFrame()
            boundness_df = pd.read_csv(ckpt_boundness) if ckpt_boundness.exists() else pd.DataFrame()
            structural_df, structural_source = _load_best_dataset(
                [
                    results_dir / "structural_metrics.csv",
                    ckpt_structural,
                ],
                required_columns={"mutation", "replicate"},
            )
            if structural_source is not None:
                logging.info(f"  Using structural source: {structural_source}")
            else:
                logging.warning("  No usable structural metrics source found.")
        except Exception as exc:
            logging.error(f"Failed to load checkpointed data: {exc}")
            return 1

        # Merge data for final outputs
        if not mmgbsa_df.empty and not structural_df.empty:
            from src.analysis.result_collector import (
                compute_binding_ddg,
                merge_with_structural_metrics,
            )

            ddg_df = compute_binding_ddg(mmgbsa_df)
            ddg_df = merge_with_structural_metrics(ddg_df, structural_df)
            ddg_df = _with_fold_change_alias(ddg_df)
            ddg_df.to_csv(results_dir / "ddg_full.csv", index=False)
            logging.info(f"  Saved: {results_dir / 'ddg_full.csv'}")
        else:
            logging.warning("Skipping ddG merge - missing MM/GBSA or structural metrics")
            ddg_df = pd.DataFrame()

        # Save other outputs
        if not rmsd_df.empty:
            rmsd_df.to_csv(results_dir / "rmsd_ca_profiles.csv", index=False)
        if not com_df.empty:
            com_df.to_csv(results_dir / "com_distance_profiles.csv", index=False)
        if not pocket_df.empty:
            pocket_df.to_csv(results_dir / "pocket_volume_profiles.csv", index=False)
        if not boundness_df.empty:
            boundness_df.to_csv(results_dir / "boundness_qc.csv", index=False)
        if not structural_df.empty:
            _with_fold_change_alias(structural_df).to_csv(results_dir / "structural_metrics.csv", index=False)
        if not mmgbsa_df.empty:
            _with_fold_change_alias(mmgbsa_df).to_csv(results_dir / "mmgbsa_replicate_metrics.csv", index=False)

        # Generate individual plots
        from src.analysis.plotting import (
            plot_all_metrics_vs_fold_reduction,
            plot_all_mutation_pocket_volume_timeseries,
            plot_boundness_qc,
            plot_simulation_convergence,
        )
        from src.utils import project_paths

        paths = project_paths(results_dir.parent)

        try:
            logging.info("  Plot 1/5: Simulation convergence (RMSD + COM distance)")
            if not rmsd_df.empty or not com_df.empty:
                plot_simulation_convergence(rmsd_df, com_df, paths)
                logging.info(f"    Saved: {plots_dir / 'rmsd_convergence.png'}")
                logging.info(f"    Saved: {plots_dir / 'com_distance_convergence.png'}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        try:
            logging.info("  Plot 2/6: Boundness QC")
            if not boundness_df.empty:
                plot_boundness_qc(boundness_df, paths)
                logging.info(f"    Saved: {plots_dir / 'boundness_qc_min_distance.png'}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        try:
            logging.info("  Plot 3/6: All metrics vs fold change")
            if not ddg_df.empty:
                plot_all_metrics_vs_fold_reduction(ddg_df, paths)
                logging.info(f"    Saved: {plots_dir / 'all_metrics_vs_fold_reduction.png'}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        try:
            logging.info("  Plot 4/6: Pocket volume time series (WT vs mutant)")
            if not pocket_df.empty:
                plot_all_mutation_pocket_volume_timeseries(pocket_df, paths)
                logging.info(f"    Saved: {plots_dir / 'pocket_volume_timeseries'}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        # Resistance heatmap
        try:
            logging.info("  Plot 5/6: Resistance heatmap")
            import subprocess
            result = subprocess.run(
                [sys.executable, "-m", "src.analysis.cli.plot_resistance_heatmap"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                logging.info(f"    Saved: {plots_dir / 'resistance_heatmap.png'}")
            else:
                logging.error(f"    Failed: {result.stderr}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        # MM/GBSA plots
        try:
            logging.info("  Plot 6/6: MM/GBSA component plots")
            if not mmgbsa_df.empty:
                import subprocess
                result = subprocess.run(
                    [sys.executable, "-m", "src.analysis.cli.plot_mmgbsa_tables"],
                    capture_output=True,
                    text=True,
                )
                if result.returncode == 0:
                    logging.info(f"    Saved: {plots_dir / 'mmgbsa_*.png'}")
                else:
                    logging.error(f"    Failed: {result.stderr}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        # DRM distance traces
        try:
            logging.info("  Plot 6/7: DRM sidechain distance traces")
            import subprocess
            result = subprocess.run(
                [sys.executable, "-m", "src.analysis.cli.plot_all_mutation_drm_distances"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                logging.info(f"    Saved: {plots_dir / 'drm_distances'}")
            else:
                logging.error(f"    Failed: {result.stderr}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

        # Crystal-derived DOR contact distance traces
        try:
            logging.info("  Plot 7/7: Crystal-derived DOR contact distance traces")
            import subprocess
            result = subprocess.run(
                [sys.executable, "-m", "src.analysis.cli.plot_all_mutation_dor_key_contacts"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                logging.info(f"    Saved: {plots_dir / 'dor_key_contacts'}")
            else:
                logging.error(f"    Failed: {result.stderr}")
        except Exception as exc:
            logging.error(f"    Failed: {exc}")

    logging.info("Analysis complete!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
