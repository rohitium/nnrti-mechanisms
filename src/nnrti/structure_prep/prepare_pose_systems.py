#!/usr/bin/env python3
"""Prepare MD assets for the perturbed DOR starting poses.

Companion to prepare_local_openmm_only_for_cluster, which is driven by genotype. Here
the protein is wild type throughout and the variable is the starting pose, so each task
is one pose CIF from pose_perturbation rather than one mutation.

Jitter is off by default. The panel runs used a 0.1 A perturbation to decorrelate
replicates seeded from the same pose; here the poses are already far apart, and a jitter
would blur the starting coordinate the whole experiment is defined by.

Usage
-----
    PYTHONPATH=src python -m nnrti.structure_prep.prepare_pose_systems \
        --poses alt_01 alt_02 alt_05 tilt_02 --production-ns 2
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.md.manifest import MDTask, save_manifest
from nnrti.md.openmm.md_protocol import MDProtocolConfig, prepare_md_assets
from nnrti.structure_prep.config import dor_4ncg_spec
from nnrti.utils import ensure_dirs


def prepare_pose_systems(
    pose_cifs: list[Path],
    *,
    results_dir: Path,
    manifest_path: Path,
    replicates: int = 1,
    production_ns: float = 2.0,
    jitter_angstrom: float = 0.0,
) -> list[MDTask]:
    from nnrti.md.openmm.minimizer import minimize_system
    from nnrti.md.openmm.require import require_module
    from nnrti.md.openmm.structure import minimize_with_restraints

    app = require_module("openmm.app")
    spec = dor_4ncg_spec(REPO_ROOT)
    cfg = MDProtocolConfig(production_ns=production_ns)
    ensure_dirs([results_dir, manifest_path.parent])

    tasks: list[MDTask] = []
    for task_id, cif_path in enumerate(pose_cifs):
        pose = cif_path.stem
        safe_label = f"pose_{pose}"
        for replicate in range(1, replicates + 1):
            run_dir = results_dir / safe_label / f"rep_{replicate:02d}"
            assets_dir = run_dir / "assets"
            ensure_dirs([run_dir, assets_dir])

            min_pdb = run_dir / f"{safe_label}_minimized_rep{replicate:02d}.pdb"
            if not min_pdb.exists():
                topology, positions, forcefield = minimize_with_restraints(
                    cif_path=cif_path,
                    ligand_resname=spec.structure.ligand_resname,
                    ligand_sdf=spec.structure.ligand_sdf,
                    restraint_radius_angstrom=spec.restraint_radius_angstrom,
                    restraint_k_kj_mol_nm2=spec.restraint_k_kj_mol_nm2,
                    output_path=min_pdb,
                    jitter_seed=None,
                    jitter_angstrom=jitter_angstrom,
                )
                _, positions = minimize_system(
                    topology, positions, forcefield,
                    restraint_indices=[], restraint_k_kj_mol_nm2=0.0,
                )
                with open(min_pdb, "w") as handle:
                    app.PDBFile.writeFile(topology, positions, handle)

            system_xml = assets_dir / f"{safe_label}_md_rep{replicate:02d}_system.xml"
            topology_pdb = assets_dir / f"{safe_label}_md_rep{replicate:02d}_start.pdb"
            if not (system_xml.exists() and topology_pdb.exists()):
                prepare_md_assets(
                    minimized_pdb_path=min_pdb,
                    ligand_resname=spec.structure.ligand_resname,
                    ligand_sdf=spec.structure.ligand_sdf,
                    topology_pdb_path=topology_pdb,
                    system_xml_path=system_xml,
                    config=cfg,
                )

            tasks.append(
                MDTask(
                    task_id=len(tasks),
                    structure="DOR",
                    mutation=f"WT/{pose}",
                    safe_label=safe_label,
                    replicate=replicate,
                    minimized_pdb=str(min_pdb.resolve()),
                    ligand_sdf=str(spec.structure.ligand_sdf.resolve()),
                    ligand_resname=spec.structure.ligand_resname,
                    fold_reduction=None,
                    output_json=str(run_dir / f"{safe_label}_rep{replicate:02d}.json"),
                    input_cif=str(cif_path.resolve()),
                    jitter_seed=None,
                    jitter_angstrom=jitter_angstrom,
                    restraint_radius=spec.restraint_radius_angstrom,
                    restraint_k=spec.restraint_k_kj_mol_nm2,
                    prepared_topology_pdb=str(topology_pdb.resolve()),
                    prepared_system_xml=str(system_xml.resolve()),
                )
            )
            logging.info("%s rep %02d ready", safe_label, replicate)

    save_manifest(tasks, manifest_path)
    logging.info("Wrote %s with %d tasks", manifest_path, len(tasks))
    return tasks


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--poses", nargs="+", required=True)
    parser.add_argument(
        "--pose-dir", type=Path, default=REPO_ROOT / "data" / "prepared" / "dor_poses"
    )
    parser.add_argument(
        "--results-dir", type=Path, default=REPO_ROOT / "results" / "md_runs" / "pose_test"
    )
    parser.add_argument(
        "--manifest", type=Path, default=REPO_ROOT / "manifests" / "pose_test_manifest.csv"
    )
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--production-ns", type=float, default=2.0)
    parser.add_argument("--jitter-angstrom", type=float, default=0.0)
    args = parser.parse_args(argv)

    cifs = []
    for pose in args.poses:
        path = args.pose_dir / f"{pose}.cif"
        if not path.exists():
            raise SystemExit(f"pose CIF not found: {path}")
        cifs.append(path)

    prepare_pose_systems(
        cifs,
        results_dir=args.results_dir,
        manifest_path=args.manifest,
        replicates=args.replicates,
        production_ns=args.production_ns,
        jitter_angstrom=args.jitter_angstrom,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
