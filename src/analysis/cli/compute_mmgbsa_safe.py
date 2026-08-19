#!/usr/bin/env python3
"""Safe, incremental MM/GBSA computation with per-replicate checkpointing.

This script processes MM/GBSA calculations one replicate at a time,
saving results after each successful computation to prevent data loss.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gc
import json
import logging
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parents[3]


def _nonempty_path(value: object) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    return Path(text)


def _remap_to_local_workspace(candidate: Path | None) -> Path | None:
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


def _append_dedup_rows(output_path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    new_df = pd.DataFrame(rows)
    if output_path.exists():
        existing_df = pd.read_csv(output_path)
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
    else:
        combined_df = new_df
    combined_df = combined_df.drop_duplicates(subset=["mutation", "replicate"], keep="last")
    combined_df = combined_df.sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)
    combined_df.to_csv(output_path, index=False)


def _steps_per_ns(timestep_fs: float) -> int:
    return int(round((1000.0 * 1000.0) / float(timestep_fs)))


def _normalize_sample_window_ns(value: float | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except Exception:
        return None
    return parsed if parsed > 0.0 else None


def _infer_total_steps_from_state_csv(rep_dir: Path, safe_label: str, replicate: int) -> int | None:
    state_csv = rep_dir / f"{safe_label}_rep{replicate:02d}_md_state.csv"
    if not state_csv.exists():
        return None
    try:
        df = pd.read_csv(state_csv)
    except Exception:
        return None


def _infer_total_steps(row: pd.Series, rep_dir: Path, safe_label: str, replicate: int) -> int | None:
    step_candidates: list[int] = []

    output_json = _resolve_local_path(
        _nonempty_path(row.get("output_json")),
        rep_dir / f"{safe_label}_rep{replicate:02d}.json",
    )
    if output_json is not None and output_json.exists():
        try:
            payload = json.loads(output_json.read_text())
            steps = int(
                payload.get("md_production_steps_completed")
                or payload.get("md_production_steps")
                or 0
            )
            if steps > 0:
                step_candidates.append(steps)
        except Exception:
            pass

    state_steps = _infer_total_steps_from_state_csv(rep_dir, safe_label, replicate)
    if state_steps is not None and state_steps > 0:
        step_candidates.append(int(state_steps))

    if not step_candidates:
        return None
    return max(step_candidates)
    col = None
    for candidate in ('#"Step"', "Step"):
        if candidate in df.columns:
            col = candidate
            break
    if col is None or df.empty:
        return None
    try:
        return int(pd.to_numeric(df[col], errors="coerce").dropna().max())
    except Exception:
        return None


def _compute_one_task(task: dict) -> tuple[bool, dict | None, str]:
    from src.md.openmm.mmgbsa import compute_mmgbsa_from_trajectory

    mutation = task["mutation"]
    replicate = int(task["replicate"])

    try:
        mm = compute_mmgbsa_from_trajectory(
            minimized_pdb_path=Path(task["min_pdb"]),
            trajectory_dcd_path=Path(task["dcd"]),
            ligand_resname=task["ligand_resname"],
            ligand_sdf=Path(task["ligand_sdf"]),
            n_snapshots=int(task["snapshots"]),
            discard_fraction=float(task["discard_fraction"]),
            sample_window_ns=task["sample_window_ns"],
            total_time_ns=task["total_time_ns"],
            sample_last_frames=task["sample_last_frames"],
            analysis_topology_pdb_path=Path(task["analysis_topo"]),
            allowed_frames=task.get("allowed_frames"),
            snapshot_relaxation=task["snapshot_relaxation"],
            relaxation_iterations=int(task["relaxation_iterations"]),
        )
        row = {
            "structure": task["structure"],
            "mutation": mutation,
            "safe_label": task["safe_label"],
            "replicate": replicate,
            "fold_reduction": task["fold_reduction"],
            "binding_dg": mm.binding_dg_mean,
            "binding_dg_std": mm.binding_dg_std,
            "binding_dg_sem": mm.binding_dg_sem,
            "binding_dg_vdw": mm.delta_e_vdw_mean,
            "binding_dg_vdw_std": mm.delta_e_vdw_std,
            "binding_dg_vdw_sem": mm.delta_e_vdw_sem,
            "binding_dg_electrostatic": mm.delta_e_elec_mean,
            "binding_dg_electrostatic_std": mm.delta_e_elec_std,
            "binding_dg_electrostatic_sem": mm.delta_e_elec_sem,
            "binding_dg_gb": mm.delta_g_gb_mean,
            "binding_dg_gb_std": mm.delta_g_gb_std,
            "binding_dg_gb_sem": mm.delta_g_gb_sem,
            "binding_dg_sa": mm.delta_g_sa_mean,
            "binding_dg_sa_std": mm.delta_g_sa_std,
            "binding_dg_sa_sem": mm.delta_g_sa_sem,
            "mmgbsa_snapshots": mm.n_snapshots,
            "mmgbsa_discard_fraction": float(task["discard_fraction"]),
            "mmgbsa_sample_window_ns": (
                float(task["sample_window_ns"])
                if task["sample_window_ns"] is not None and float(task["sample_window_ns"]) > 0.0
                else float("nan")
            ),
            "mmgbsa_total_time_ns": (
                float(task["total_time_ns"])
                if task["total_time_ns"] is not None and float(task["total_time_ns"]) > 0.0
                else float("nan")
            ),
            "mmgbsa_time_source": task["time_source"],
            **{
                f"abs_{key}_{stat}": value
                for key, stats in getattr(mm, "absolute_terms", ())
                for stat, value in zip(("mean", "std", "sem"), stats)
            },
            "mmgbsa_snapshot_relaxation": task["snapshot_relaxation"],
            "mmgbsa_relaxation_iterations": int(task["relaxation_iterations"]),
            "mmgbsa_contact_screened": bool(task.get("allowed_frames") is not None),
            "mmgbsa_clean_frames_available": (
                int(len(task["allowed_frames"])) if task.get("allowed_frames") is not None else float("nan")
            ),
            "mmgbsa_sample_last_frames": (
                int(task["sample_last_frames"])
                if task["sample_last_frames"] is not None and int(task["sample_last_frames"]) > 0
                else float("nan")
            ),
        }
        return True, row, ""
    except Exception as exc:
        return False, None, f"{mutation} rep{replicate}: {exc}"


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Safe incremental MM/GBSA computation")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, help="Output CSV path")
    parser.add_argument("--snapshots", type=int, default=100)
    parser.add_argument("--discard-fraction", type=float, default=0.25)
    parser.add_argument(
        "--sample-window-ns",
        type=float,
        default=0.0,
        help="If > 0, sample only the last N ns. Default 0 uses the full post-discard region.",
    )
    parser.add_argument(
        "--sample-last-frames",
        type=int,
        default=0,
        help="If > 0, sample the last N saved trajectory frames. Overrides --sample-window-ns.",
    )
    parser.add_argument(
        "--snapshot-relaxation",
        choices=("unrestrained", "h_relax", "none"),
        default="unrestrained",
        help=(
            "How each snapshot is relaxed before scoring. 'unrestrained' (default) minimises with all "
            "atoms free, which relieves both hydrogen and heavy-atom overlaps; 'h_relax' is the legacy "
            "hydrogen-only relaxation with heavy atoms restrained."
        ),
    )
    parser.add_argument(
        "--relaxation-iterations",
        type=int,
        default=100,
        help="Minimisation iteration cap. Deliberately not run to convergence; see the method note.",
    )
    parser.add_argument(
        "--contact-screen-csv",
        type=Path,
        default=None,
        help=(
            "Per-frame contact screen from screen_ligand_contact_artifacts. When given, frames "
            "flagged as artifacts are excluded and the most recent surviving frames are scored."
        ),
    )
    parser.add_argument(
        "--min-clean-frames",
        type=int,
        default=5,
        help="Fail a replicate rather than score it if the screen leaves fewer clean frames than this.",
    )
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers for per-replicate MM/GBSA")
    parser.add_argument("--force", action="store_true", help="Recompute even if checkpoint exists")
    args = parser.parse_args()

    sample_last_frames = int(args.sample_last_frames) if int(args.sample_last_frames) > 0 else None
    sample_window_ns = None if sample_last_frames is not None else _normalize_sample_window_ns(args.sample_window_ns)

    ckpt_dir = args.results_dir / ".checkpoints"
    output_path = args.output or (ckpt_dir / ".checkpoint_mmgbsa_replicate_metrics.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.force and output_path.exists():
        logging.info(f"--force active: removing previous MM/GBSA checkpoint at {output_path}")
        output_path.unlink()

    # Load MD metadata
    from src.analysis.result_collector import collect_md_results

    logging.info(f"Loading MD metadata from {args.manifest}")
    md_df = collect_md_results(args.manifest, args.results_dir)

    if md_df.empty:
        logging.error("No MD results found")
        return 1

    # Load existing checkpoint if available
    existing_results = []
    if output_path.exists() and not args.force:
        try:
            existing_df = pd.read_csv(output_path)
            existing_df = existing_df.drop_duplicates(subset=["mutation", "replicate"], keep="last")
            existing_df = existing_df.sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)
            existing_df.to_csv(output_path, index=False)
            existing_results = set(
                zip(
                    existing_df["mutation"],
                    existing_df["replicate"],
                )
            )
            logging.info(f"Loaded {len(existing_results)} existing MM/GBSA results from {output_path}")
        except Exception as exc:
            logging.warning(f"Could not load existing results: {exc}")

    clean_frames: dict[tuple[str, int], list[int]] = {}
    if args.contact_screen_csv is not None:
        if not args.contact_screen_csv.exists():
            raise FileNotFoundError(args.contact_screen_csv)
        screen = pd.read_csv(args.contact_screen_csv)
        screen = screen[(screen["status"] == "ok") & screen["is_clean"].astype(bool)]
        for (mut, rep), grp in screen.groupby(["mutation", "replicate"]):
            clean_frames[(str(mut), int(rep))] = sorted(int(f) for f in grp["frame"])
        logging.info(
            "Loaded contact screen for %d runs from %s", len(clean_frames), args.contact_screen_csv
        )

    # Build runnable task list
    tasks: list[dict] = []
    total = len(md_df)

    for idx, (_, row) in enumerate(md_df.iterrows(), 1):
        mutation = str(row["mutation"])
        replicate = int(row["replicate"])

        # Skip if already computed
        if (mutation, replicate) in existing_results:
            logging.info(f"[{idx}/{total}] Skipping {mutation} rep{replicate} (already computed)")
            continue

        logging.info(f"[{idx}/{total}] Processing {mutation} rep{replicate}")

        # Get file paths
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])

        min_pdb = _resolve_local_path(
            _nonempty_path(row.get("minimized_pdb")),
            rep_dir / f"{safe}_minimized_rep{rep:02d}.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )
        analysis_topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        ligand_sdf = _resolve_local_path(_nonempty_path(row.get("ligand_sdf")))

        # Validate paths
        if None in {min_pdb, dcd, analysis_topo, ligand_sdf}:
            logging.warning(f"  Missing inputs for {mutation} rep{replicate}")
            continue
        if not min_pdb.exists() or not dcd.exists() or not analysis_topo.exists() or not ligand_sdf.exists():
            logging.warning(f"  Unavailable paths for {mutation} rep{replicate}")
            continue

        # Derive fallback discard fraction from trajectory span when possible.
        # This is used only if DCD timing metadata is invalid inside MM/GBSA code.
        discard_fraction = float(args.discard_fraction)
        total_time_ns = None
        time_source = "trajectory_dt"
        if sample_window_ns is not None:
            total_steps = _infer_total_steps(row, rep_dir, safe, rep)
            if total_steps and total_steps > 0:
                total_time_ns = float(total_steps) * float(args.timestep_fs) / 1_000_000.0
                time_source = "run_artifacts"
                window_steps = _steps_per_ns(args.timestep_fs) * float(sample_window_ns)
                keep_fraction = min(1.0, float(window_steps) / float(total_steps))
                discard_fraction = max(0.0, min(0.95, 1.0 - keep_fraction))
                logging.info(
                    f"  {mutation} rep{replicate}: fallback discard_fraction={discard_fraction:.4f} "
                    f"(window_ns={sample_window_ns}, total_steps={total_steps}, total_ns={total_time_ns:.3f})"
                )
            else:
                logging.warning(
                    f"  {mutation} rep{replicate}: could not infer total steps; "
                    f"using discard_fraction={discard_fraction:.4f} fallback"
                )

        allowed = None
        if args.contact_screen_csv is not None:
            allowed = clean_frames.get((mutation, replicate))
            if allowed is None:
                logging.warning("  No screen rows for %s rep%d; skipping", mutation, replicate)
                continue
            if len(allowed) < int(args.min_clean_frames):
                logging.error(
                    "  %s rep%d has only %d clean frames (min %d); skipping",
                    mutation, replicate, len(allowed), args.min_clean_frames,
                )
                continue

        tasks.append(
            {
                "allowed_frames": allowed,
                "snapshot_relaxation": args.snapshot_relaxation,
                "relaxation_iterations": args.relaxation_iterations,
                "structure": row["structure"],
                "mutation": mutation,
                "safe_label": safe,
                "replicate": replicate,
                "fold_reduction": row.get("fold_reduction"),
                "min_pdb": str(min_pdb),
                "dcd": str(dcd),
                "analysis_topo": str(analysis_topo),
                "ligand_sdf": str(ligand_sdf),
                "ligand_resname": args.ligand_resname,
                "snapshots": int(args.snapshots),
                "discard_fraction": discard_fraction,
                "sample_window_ns": sample_window_ns,
                "total_time_ns": total_time_ns,
                "time_source": time_source,
                "sample_last_frames": sample_last_frames,
            }
        )

    if not tasks:
        logging.info("No MM/GBSA tasks to run.")
        return 0

    workers = max(1, int(args.workers))
    logging.info(f"Running {len(tasks)} MM/GBSA tasks with workers={workers}")

    if workers == 1:
        for i, task in enumerate(tasks, 1):
            ok, row, err = _compute_one_task(task)
            if ok and row is not None:
                _append_dedup_rows(output_path, [row])
                logging.info(
                    f"[{i}/{len(tasks)}] ✓ {task['mutation']} rep{task['replicate']} "
                    f"(ΔG={row['binding_dg']:.2f} ± {row['binding_dg_std']:.2f})"
                )
            else:
                logging.error(f"[{i}/{len(tasks)}] ✗ {err}")
            gc.collect()
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_compute_one_task, task): task for task in tasks}
            done = 0
            for fut in cf.as_completed(futures):
                done += 1
                task = futures[fut]
                try:
                    ok, row, err = fut.result()
                except Exception as exc:
                    logging.error(f"[{done}/{len(tasks)}] ✗ {task['mutation']} rep{task['replicate']}: {exc}")
                    continue
                if ok and row is not None:
                    _append_dedup_rows(output_path, [row])
                    logging.info(
                        f"[{done}/{len(tasks)}] ✓ {task['mutation']} rep{task['replicate']} "
                        f"(ΔG={row['binding_dg']:.2f} ± {row['binding_dg_std']:.2f})"
                    )
                else:
                    logging.error(f"[{done}/{len(tasks)}] ✗ {err}")

    logging.info(f"MM/GBSA computation complete. Results saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
