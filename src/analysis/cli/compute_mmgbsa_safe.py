#!/usr/bin/env python3
"""Safe, incremental MM/GBSA computation with per-replicate checkpointing.

This script processes MM/GBSA calculations one replicate at a time,
saving results after each successful computation to prevent data loss.
"""
from __future__ import annotations

import argparse
import gc
import logging
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parents[3]


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Safe incremental MM/GBSA computation")
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, help="Output CSV path")
    parser.add_argument("--snapshots", type=int, default=100)
    parser.add_argument("--discard-fraction", type=float, default=0.25)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--force", action="store_true", help="Recompute even if checkpoint exists")
    args = parser.parse_args()

    ckpt_dir = args.results_dir / ".checkpoints"
    output_path = args.output or (ckpt_dir / ".checkpoint_mmgbsa_replicate_metrics.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

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
            existing_results = set(
                zip(
                    existing_df["mutation"],
                    existing_df["replicate"],
                )
            )
            logging.info(f"Loaded {len(existing_results)} existing MM/GBSA results from {output_path}")
        except Exception as exc:
            logging.warning(f"Could not load existing results: {exc}")

    # Helper functions
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

    # Process each replicate
    from src.md.openmm.mmgbsa import compute_mmgbsa_from_trajectory

    results = []
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

        # Compute MM/GBSA
        try:
            mm = compute_mmgbsa_from_trajectory(
                minimized_pdb_path=min_pdb,
                trajectory_dcd_path=dcd,
                ligand_resname=args.ligand_resname,
                ligand_sdf=ligand_sdf,
                n_snapshots=args.snapshots,
                discard_fraction=args.discard_fraction,
                analysis_topology_pdb_path=analysis_topo,
            )

            result_row = {
                "structure": row["structure"],
                "mutation": mutation,
                "safe_label": safe,
                "replicate": replicate,
                "fold_reduction": row.get("fold_reduction"),
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
            }
            results.append(result_row)
            logging.info(f"  ✓ Success: ΔG = {mm.binding_dg_mean:.2f} ± {mm.binding_dg_std:.2f} kJ/mol")

            # Save checkpoint after each successful computation
            if results:
                new_df = pd.DataFrame(results)
                if output_path.exists():
                    existing_df = pd.read_csv(output_path)
                    combined_df = pd.concat([existing_df, new_df], ignore_index=True)
                else:
                    combined_df = new_df
                combined_df.to_csv(output_path, index=False)
                logging.info(f"  Checkpoint saved: {output_path}")
                results.clear()  # Clear to avoid duplicates

        except Exception as exc:
            logging.error(f"  ✗ Failed: {exc}")

        # Force garbage collection
        gc.collect()

    logging.info(f"MM/GBSA computation complete. Results saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
