#!/usr/bin/env python
"""Prepare apo (ligand-free) MD assets from existing holo minimized PDBs.

Scans results/md_runs/ for completed holo minimized PDBs, strips the DOR
ligand (resname 2KW), builds amber-only OpenMM systems, and writes an apo
manifest CSV.  The resulting assets use the same run_prepared_md / worker
pathway as the holo simulations.

Usage
-----
    python -m src.dor_md_pipeline_apo \\
        --mutations WT F227C V106A "V106A+P225H" "K103N+M230L" "A98G+F227C" "V106I+F227C" \\
        --holo-runs results/md_runs \\
        --apo-runs results/apo_md_runs \\
        --manifest results/apo_md_manifest.csv

All three arguments have the above defaults; the mutation list defaults to the
seven priority apo systems identified in the hypothesis analysis.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------
# Priority apo mutations (hypothesis-driven)
# -----------------------------------------------------------------------
_DEFAULT_MUTATIONS = [
    "WT",
    "F227C",
    "V106A",
    "V106A+P225H",
    "K103N+M230L",
    "A98G+F227C",
    "V106I+F227C",
]

# Ligand residue name in the 4NCG crystal / holo PDBs.
_LIGAND_RESNAME = "2KW"


def _safe_label(mutation: str) -> str:
    """Convert a mutation string to a filesystem-safe label (matching holo convention)."""
    return mutation.lower().replace("+", "_").replace(" ", "")


def prep_apo_systems(
    mutations: list[str],
    holo_runs_root: Path,
    apo_runs_root: Path,
    manifest_path: Path,
) -> None:
    from .md.manifest import MDTask, save_manifest
    from .md.openmm.md_protocol import MDProtocolConfig, prepare_apo_md_assets

    cfg = MDProtocolConfig()
    tasks: list[MDTask] = []
    task_id = 0

    for mutation in mutations:
        safe = _safe_label(mutation)
        holo_mut_dir = holo_runs_root / safe

        if not holo_mut_dir.is_dir():
            logger.warning("Holo dir not found for %s: %s — skipping", mutation, holo_mut_dir)
            continue

        rep_dirs = sorted(holo_mut_dir.glob("rep_*"))
        if not rep_dirs:
            logger.warning("No replicate dirs under %s — skipping", holo_mut_dir)
            continue

        for rep_dir in rep_dirs:
            rep_str = rep_dir.name  # e.g. "rep_01"
            rep_int = int(rep_str.split("_")[1])

            # Find the holo minimized PDB for this replicate.
            holo_min_pdb = rep_dir / f"{safe}_minimized_rep{rep_int:02d}.pdb"
            if not holo_min_pdb.exists():
                logger.warning("Missing holo minimized PDB: %s — skipping replicate", holo_min_pdb)
                continue

            # Apo output directories mirror the holo layout under apo_runs_root.
            apo_rep_dir = apo_runs_root / safe / rep_str
            apo_assets_dir = apo_rep_dir / "assets"
            apo_rep_dir.mkdir(parents=True, exist_ok=True)
            apo_assets_dir.mkdir(parents=True, exist_ok=True)

            topology_pdb = apo_assets_dir / f"{safe}_apo_md_rep{rep_int:02d}_start.pdb"
            system_xml = apo_assets_dir / f"{safe}_apo_md_rep{rep_int:02d}_system.xml"
            output_json = apo_rep_dir / f"{safe}_apo_rep{rep_int:02d}.json"

            if topology_pdb.exists() and system_xml.exists():
                logger.info("Apo assets already exist for %s rep%02d — skipping prep", mutation, rep_int)
            else:
                logger.info("Preparing apo assets for %s rep%02d …", mutation, rep_int)
                prepare_apo_md_assets(
                    minimized_pdb_path=holo_min_pdb,
                    ligand_resname=_LIGAND_RESNAME,
                    topology_pdb_path=topology_pdb,
                    system_xml_path=system_xml,
                    config=cfg,
                )
                logger.info("  → wrote %s + %s", topology_pdb, system_xml)

            tasks.append(
                MDTask(
                    task_id=task_id,
                    structure="DOR",
                    mutation=mutation,
                    safe_label=safe,
                    replicate=rep_int,
                    minimized_pdb=str(holo_min_pdb.resolve()),
                    ligand_sdf="",           # no ligand in apo
                    ligand_resname="",
                    fold_reduction=None,
                    output_json=str(output_json.resolve()),
                    leg="apo",
                    prepared_topology_pdb=str(topology_pdb.resolve()),
                    prepared_system_xml=str(system_xml.resolve()),
                )
            )
            task_id += 1

    if not tasks:
        logger.error("No apo tasks were prepared — check mutation names and holo-runs path.")
        sys.exit(1)

    save_manifest(tasks, manifest_path)
    logger.info("Wrote apo manifest with %d tasks to %s", len(tasks), manifest_path)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(
        description="Prepare apo MD assets from existing holo minimized PDBs.",
    )
    parser.add_argument(
        "--mutations",
        nargs="+",
        default=_DEFAULT_MUTATIONS,
        metavar="MUT",
        help="Mutation labels to prepare (default: 7 priority systems).",
    )
    parser.add_argument(
        "--holo-runs",
        type=Path,
        default=Path("results/md_runs"),
        help="Root directory of holo MD run results (default: results/md_runs).",
    )
    parser.add_argument(
        "--apo-runs",
        type=Path,
        default=Path("results/apo_md_runs"),
        help="Root directory for apo MD run results (default: results/apo_md_runs).",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("results/apo_md_manifest.csv"),
        help="Output apo manifest CSV path.",
    )
    args = parser.parse_args(argv)

    prep_apo_systems(
        mutations=args.mutations,
        holo_runs_root=args.holo_runs,
        apo_runs_root=args.apo_runs,
        manifest_path=args.manifest,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
