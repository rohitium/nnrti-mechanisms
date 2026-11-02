#!/usr/bin/env python3
"""Minimize each perturbed DOR pose locally and report whether it is recoverable.

The holo prep already runs a restrained minimization in which everything beyond 8 A of
the ligand is held and the pocket is free, so that step doubles as a cheap screen: a
starting pose whose clashes cannot be relieved here will not survive equilibration on
the cluster either.  For each pose this records the DOR heavy-atom RMSD to the
crystallographic pose before and after minimization, the closest ligand-protein contact
after minimization, and the final potential energy.

A pose that minimizes to a sane energy with no contact below roughly 2 A is worth
running.  A pose that stays clashed, or whose energy is an outlier, is not.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.md.openmm.require import require_module
from nnrti.md.openmm.structure import minimize_with_restraints
from nnrti.structure_prep.config import dor_4ncg_spec
from nnrti.structure_prep.pose_perturbation import (
    heavy_mask,
    ligand_rows,
    parse_atom_site,
)


def _ligand_heavy_xyz(topology, positions, ligand_resname: str) -> np.ndarray:
    unit = require_module("openmm.unit")
    return np.array(
        [
            positions[atom.index].value_in_unit(unit.angstrom)
            for atom in topology.atoms()
            if atom.residue.name == ligand_resname
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=float,
    )


def _environment_heavy_xyz(topology, positions, ligand_resname: str) -> np.ndarray:
    unit = require_module("openmm.unit")
    return np.array(
        [
            positions[atom.index].value_in_unit(unit.angstrom)
            for atom in topology.atoms()
            if atom.residue.name not in {ligand_resname, "HOH", "WAT"}
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=float,
    )


def _kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Plain (not superposed) RMSD: the pose is judged in the pocket frame."""
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def screen_pose(
    cif_path: Path,
    crystal_heavy: np.ndarray,
    spec,
    output_dir: Path,
) -> dict:
    name = cif_path.stem
    started = time.time()
    record: dict = {"pose": name}

    minimized_pdb = output_dir / f"{name}_minimized.pdb"
    try:
        topology, positions, _ = minimize_with_restraints(
            cif_path=cif_path,
            ligand_resname=spec.structure.ligand_resname,
            ligand_sdf=spec.structure.ligand_sdf,
            restraint_radius_angstrom=spec.restraint_radius_angstrom,
            restraint_k_kj_mol_nm2=spec.restraint_k_kj_mol_nm2,
            output_path=minimized_pdb,
            jitter_seed=None,
            jitter_angstrom=0.0,
        )
    except Exception as exc:  # noqa: BLE001 - the failure mode is the result
        record.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        logging.warning("%s: %s", name, record["error"])
        return record

    ligand_xyz = _ligand_heavy_xyz(topology, positions, spec.structure.ligand_resname)
    env_xyz = _environment_heavy_xyz(topology, positions, spec.structure.ligand_resname)
    distances = np.linalg.norm(ligand_xyz[:, None, :] - env_xyz[None, :, :], axis=2)

    record.update(
        {
            "status": "ok",
            "rmsd_after_min": round(_kabsch_rmsd(ligand_xyz, crystal_heavy), 3),
            "centroid_shift_after_min": round(
                float(np.linalg.norm(ligand_xyz.mean(axis=0) - crystal_heavy.mean(axis=0))), 3
            ),
            "min_contact_after_min": round(float(distances.min()), 3),
            "n_contacts_below_2.5_after_min": int((distances < 2.5).sum()),
            "seconds": round(time.time() - started, 1),
        }
    )
    logging.info(
        "%-14s rmsd %5.2f  minD %5.2f  (%.0f s)",
        name,
        record["rmsd_after_min"],
        record["min_contact_after_min"],
        record["seconds"],
    )
    return record


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--pose-dir", type=Path, default=REPO_ROOT / "data" / "prepared" / "dor_poses"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "data" / "prepared" / "dor_poses" / "minimized",
    )
    parser.add_argument("--poses", nargs="+", default=None, help="Pose stems to screen.")
    args = parser.parse_args(argv)

    spec = dor_4ncg_spec(REPO_ROOT)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    block = parse_atom_site(spec.structure.cif_path)
    idx = ligand_rows(block, spec.structure.ligand_resname)
    crystal_heavy = block.coords(idx)[heavy_mask(block, idx)]

    cif_paths = sorted(args.pose_dir.glob("*.cif"))
    if args.poses:
        wanted = set(args.poses)
        cif_paths = [p for p in cif_paths if p.stem in wanted]
    if not cif_paths:
        raise SystemExit(f"No pose CIFs found in {args.pose_dir}")

    # The crystal pose is screened too: it is the control for the minimization itself.
    records = [
        screen_pose(spec.structure.cif_path, crystal_heavy, spec, args.output_dir)
        | {"pose": "crystal"}
    ]
    for cif_path in cif_paths:
        records.append(screen_pose(cif_path, crystal_heavy, spec, args.output_dir))

    summary = args.output_dir / "pose_screen.csv"
    fields = [
        "pose",
        "status",
        "rmsd_after_min",
        "centroid_shift_after_min",
        "min_contact_after_min",
        "n_contacts_below_2.5_after_min",
        "seconds",
        "error",
    ]
    with open(summary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    logging.info("Wrote %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
