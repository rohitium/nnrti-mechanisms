#!/usr/bin/env python3
"""Parallel structural metrics computation with per-replicate checkpointing.

This script computes ensemble metrics (contacts, H-bonds, pocket volume) in parallel
across multiple processes, saving results after each replicate to prevent data loss.
"""
from __future__ import annotations

import argparse
import gc
import logging
import multiprocessing as mp
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parents[3]


def _nonempty_path(value: object) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    return Path(text)


def _remap_to_local_workspace(candidate: Path | None) -> Path | None:
    """Map stale absolute paths (e.g., /scratch/.../nnrti-mechanisms/...) to local repo."""
    if candidate is None:
        return None
    if candidate.exists():
        return candidate

    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker not in text:
        return candidate

    rel = text.split(marker, 1)[1]
    mapped = project_root / rel
    if mapped.exists():
        return mapped
    return candidate


def _resolve_local_path(candidate: Path | None, fallback: Path | None = None) -> Path | None:
    candidate = _remap_to_local_workspace(candidate)
    fallback = _remap_to_local_workspace(fallback)
    if candidate is not None and candidate.exists():
        return candidate
    if fallback is not None and fallback.exists():
        return fallback
    return candidate or fallback


def _infer_rep_dir(row: pd.Series) -> Path:
    for key in ("analysis_dcd", "trajectory_dcd", "prepared_topology_pdb", "minimized_pdb"):
        val = str(row.get(key) or "").strip()
        if val:
            return Path(val).parent
    return Path(".")


def _normalize_sample_window_ns(value: float | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except Exception:
        return None
    return parsed if parsed > 0.0 else None


def _compute_one_replicate(args_tuple):
    """Worker function for parallel processing."""
    (row_dict, ligand_resname, frame_stride, max_frames, sample_window_ns) = args_tuple

    # Import here to avoid issues with multiprocessing
    from src.analysis.metrics import compute_ensemble_metrics

    row = pd.Series(row_dict)
    rep_dir = _infer_rep_dir(row)
    safe = str(row["safe_label"])
    rep = int(row["replicate"])
    mutation = row["mutation"]

    topo = _resolve_local_path(
        _nonempty_path(row.get("analysis_topology_pdb")),
        rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
    )
    dcd = _resolve_local_path(
        _nonempty_path(row.get("analysis_dcd")),
        rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
    )

    if topo is None or dcd is None or not topo.exists() or not dcd.exists():
        return {
            "mutation": mutation,
            "replicate": rep,
            "error": "missing_files",
        }

    try:
        ens = compute_ensemble_metrics(
            topology_pdb_path=topo,
            trajectory_dcd_path=dcd,
            ligand_resname=ligand_resname,
            frame_stride=frame_stride,
            max_frames=max_frames,
            sample_window_ns=sample_window_ns,
        )
        return {
            "structure": row["structure"],
            "mutation": mutation,
            "safe_label": safe,
            "replicate": rep,
            "contact_count": ens.contact_count_mean,
            "contact_count_std": ens.contact_count_std,
            "hbond_count": ens.hbond_count_mean,
            "hbond_count_std": ens.hbond_count_std,
            "pocket_volume_proxy": ens.pocket_volume_proxy_mean,
            "pocket_volume_proxy_std": ens.pocket_volume_proxy_std,
            "metric_n_frames": ens.n_frames,
            "metric_source": "trajectory",
            "metric_sample_window_ns": (
                float(sample_window_ns)
                if sample_window_ns is not None and float(sample_window_ns) > 0.0
                else float("nan")
            ),
            "fold_reduction": row["fold_reduction"],
            "error": None,
        }
    except Exception as exc:
        return {
            "mutation": mutation,
            "replicate": rep,
            "error": str(exc),
        }


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Parallel structural metrics computation")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, help="Output CSV path")
    parser.add_argument("--frame-stride", type=int, default=5)
    parser.add_argument("--max-frames", type=int, default=200)
    parser.add_argument(
        "--sample-window-ns",
        type=float,
        default=0.0,
        help="If > 0, sample only the last N ns. Default 0 uses all sampled frames.",
    )
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--workers", type=int, default=None, help="Number of parallel workers (default: CPU count)")
    parser.add_argument("--force", action="store_true", help="Recompute even if checkpoint exists")
    args = parser.parse_args()

    sample_window_ns = _normalize_sample_window_ns(args.sample_window_ns)

    ckpt_dir = args.results_dir / ".checkpoints"
    output_path = args.output or (ckpt_dir / ".checkpoint_structural_metrics.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load MD metadata
    from src.analysis.result_collector import collect_md_results

    logging.info(f"Loading MD metadata from {args.manifest}")
    md_df = collect_md_results(args.manifest, args.results_dir)

    if md_df.empty:
        logging.error("No MD results found")
        return 1

    # Load existing checkpoint if available
    existing_results = set()
    existing_df = pd.DataFrame()
    if output_path.exists() and not args.force:
        try:
            existing_df = pd.read_csv(output_path)
            # Filter out error entries
            valid_df = existing_df[existing_df["error"].isna()]
            existing_results = set(zip(valid_df["mutation"], valid_df["replicate"]))
            logging.info(f"Loaded {len(existing_results)} existing structural metrics from {output_path}")
        except Exception as exc:
            logging.warning(f"Could not load existing results: {exc}")

    # Filter to only replicates that need processing
    to_process = []
    for idx, row in md_df.iterrows():
        mutation = str(row["mutation"])
        replicate = int(row["replicate"])
        if (mutation, replicate) not in existing_results:
            to_process.append(
                (
                    row.to_dict(),
                    args.ligand_resname,
                    args.frame_stride,
                    args.max_frames,
                    sample_window_ns,
                )
            )

    if not to_process:
        logging.info("All structural metrics already computed!")
        return 0

    logging.info(f"Processing {len(to_process)} replicates using {args.workers or mp.cpu_count()} workers")

    # Process in parallel
    results = []
    if args.workers == 1:
        # Serial processing for debugging
        for i, task_args in enumerate(to_process, 1):
            logging.info(f"[{i}/{len(to_process)}] Processing {task_args[0]['mutation']} rep{task_args[0]['replicate']}")
            result = _compute_one_replicate(task_args)
            if result.get("error") is None:
                results.append(result)
                logging.info(f"  ✓ Success")
            else:
                logging.warning(f"  ✗ Failed: {result['error']}")

            # Save checkpoint after each success (in serial mode)
            if results:
                new_df = pd.DataFrame(results)
                if not existing_df.empty:
                    combined_df = pd.concat([existing_df, new_df], ignore_index=True)
                else:
                    combined_df = new_df
                combined_df.to_csv(output_path, index=False)
                results.clear()
    else:
        # Parallel processing
        with mp.Pool(processes=args.workers) as pool:
            for i, result in enumerate(pool.imap_unordered(_compute_one_replicate, to_process), 1):
                if result.get("error") is None:
                    results.append(result)
                    logging.info(f"[{i}/{len(to_process)}] ✓ {result['mutation']} rep{result['replicate']}")
                else:
                    logging.warning(f"[{i}/{len(to_process)}] ✗ {result.get('mutation', '?')} rep{result.get('replicate', '?')}: {result['error']}")

                # Save checkpoint every 5 successful results
                if len(results) >= 5:
                    new_df = pd.DataFrame(results)
                    if not existing_df.empty:
                        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
                    else:
                        combined_df = new_df
                    combined_df.to_csv(output_path, index=False)
                    logging.info(f"  Checkpoint saved ({len(results)} new results)")
                    existing_df = combined_df  # Update for next batch
                    results.clear()

    # Save any remaining results
    if results:
        new_df = pd.DataFrame(results)
        if not existing_df.empty:
            combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        else:
            combined_df = new_df
        combined_df.to_csv(output_path, index=False)
        logging.info(f"Final checkpoint saved")

    logging.info(f"Structural metrics computation complete. Results saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
