#!/usr/bin/env python3
"""Pocket-frame pose coordinates for DOR, and where the perturbed starts sit.

The published coordinates cannot answer whether a starting pose is inside the
sampled region.  ``dor_rmsd_angstrom`` superposes on the ligand itself, so rigid-body
motion is removed by construction, and ``com_distance_angstrom`` measures the DOR
centroid against the centroid of the whole protein, which barely moves with pose.
A rotation about the ligand centroid leaves both unchanged.

This computes two coordinates that do resolve pose, both after superposing each frame
on the NNIBP-lining Ca atoms so that motion is measured relative to the pocket rather
than to the laboratory frame:

``centroid_displacement_angstrom``
    Distance from the crystallographic DOR centroid.  Pure translation.

``long_axis_tilt_deg``
    Angle between the DOR long axis and its crystallographic orientation, where the
    long axis runs from the chlorine of the chlorocyanophenyl ring to the centroid of
    the triazolinone N/O atoms.  A defined molecular vector rather than an inertial
    eigenvector, so an end-for-end flip reads as 180 degrees instead of folding back
    to 0.  Pure rotation.

``pose_rmsd_angstrom`` is also reported: heavy-atom RMSD to the crystallographic pose
with no ligand superposition, which is the scalar the minimization screen used.

Reference frame: the NNIBP Ca atoms of the prepared crystal structure.  Reference pose:
the 4NCG ligand coordinates themselves.  The two differ by the 0.39 A of Ca movement
during preparation, which is small against everything measured here.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.analysis.pbc import apply_mdtraj_pbc_correction, load_mdtraj_trajectory
from nnrti.structure_prep.pose_perturbation import (
    heavy_mask,
    ligand_rows,
    parse_atom_site,
)

LIG_RESNAME = "2KW"

# NNIBP-lining p66 residues, topology resids (canonical - 3).  Source: Cilento, Kirby
# & Sarafianos.  MDTraj `resid` is a global residue index and chain 0 is p66, so these
# index directly; the p51 contribution (E138) is left out because it falls in another
# chain and 15 Ca atoms already define the frame.
_NNBP_CA_RESIDS_P66: tuple[int, ...] = (
    97, 98, 100, 103, 104, 105,
    176, 178, 185, 186, 187,
    224, 226, 231, 315,
)

# Indices into the 29-atom heavy-atom list, which follows 4NCG CIF order.
_CL_INDEX = 20
_TRIAZOLINONE_INDICES = (25, 26, 27, 28)


def _nnibp_ca_selection(topology) -> np.ndarray:
    query = "chainid 0 and name CA and (" + " or ".join(
        f"resid {r}" for r in _NNBP_CA_RESIDS_P66
    ) + ")"
    idx = topology.select(query)
    if len(idx) != len(_NNBP_CA_RESIDS_P66):
        raise ValueError(
            f"NNIBP Ca selection returned {len(idx)}, expected {len(_NNBP_CA_RESIDS_P66)}"
        )
    return idx


def _ligand_heavy_selection(topology) -> np.ndarray:
    idx = topology.select(f"resname '{LIG_RESNAME}'")
    heavy = np.array(
        [i for i in idx if topology.atom(i).element.symbol != "H"], dtype=int
    )
    if len(heavy) != 29:
        raise ValueError(f"expected 29 DOR heavy atoms, found {len(heavy)}")
    return heavy


def _long_axis(xyz: np.ndarray) -> np.ndarray:
    vector = xyz[_CL_INDEX] - xyz[list(_TRIAZOLINONE_INDICES)].mean(axis=0)
    return vector / np.linalg.norm(vector)


def _kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    return rot, mc, tc


def pose_coordinates(ligand_xyz: np.ndarray, reference_xyz: np.ndarray) -> dict:
    """Both inputs already in the pocket frame, angstrom, 29 heavy atoms."""
    centroid = ligand_xyz.mean(axis=0)
    ref_centroid = reference_xyz.mean(axis=0)
    cosine = float(np.clip(np.dot(_long_axis(ligand_xyz), _long_axis(reference_xyz)), -1.0, 1.0))
    return {
        "centroid_displacement_angstrom": float(np.linalg.norm(centroid - ref_centroid)),
        "long_axis_tilt_deg": float(np.degrees(np.arccos(cosine))),
        "pose_rmsd_angstrom": float(
            np.sqrt(np.mean(np.sum((ligand_xyz - reference_xyz) ** 2, axis=1)))
        ),
    }


def _reference(prepared_crystal_pdb: Path, crystal_cif: Path) -> tuple[np.ndarray, np.ndarray]:
    """NNIBP Ca frame from the prepared structure, ligand pose from the crystal CIF."""
    import mdtraj as md

    prepared = md.load(str(prepared_crystal_pdb))
    frame_xyz = prepared.xyz[0, _nnibp_ca_selection(prepared.topology), :] * 10.0

    block = parse_atom_site(crystal_cif)
    idx = ligand_rows(block, LIG_RESNAME)
    ligand_xyz = block.coords(idx)[heavy_mask(block, idx)]
    return frame_xyz, ligand_xyz


def starting_pose_rows(pose_dir: Path, crystal_cif: Path, reference_ligand: np.ndarray) -> list[dict]:
    """The 12 perturbed starts, before minimization.

    Every pose CIF carries the unmodified 4NCG protein, so all of them and the crystal
    share one pocket frame and no superposition is needed.
    """
    rows = [{"pose": "crystal", "family": "reference", **pose_coordinates(reference_ligand, reference_ligand)}]
    for cif_path in sorted(pose_dir.glob("*.cif")):
        block = parse_atom_site(cif_path)
        idx = ligand_rows(block, LIG_RESNAME)
        ligand_xyz = block.coords(idx)[heavy_mask(block, idx)]
        family = "rotation" if cif_path.stem.startswith("rot_") else "alternative"
        rows.append(
            {"pose": cif_path.stem, "family": family, **pose_coordinates(ligand_xyz, reference_ligand)}
        )
    return rows


def trajectory_rows(
    run_dir: Path,
    label: str,
    reference_frame: np.ndarray,
    reference_ligand: np.ndarray,
    stride: int,
) -> list[dict]:
    import mdtraj as md

    rows: list[dict] = []
    for dcd_path in sorted(run_dir.glob("rep_*/*_analysis.dcd")):
        replicate = int(dcd_path.parent.name.split("_")[-1])
        topo_path = dcd_path.with_name(dcd_path.stem + "_topology.pdb")
        if not topo_path.exists():
            logging.warning("no topology beside %s, skipping", dcd_path.name)
            continue

        traj = load_mdtraj_trajectory(dcd_path, topo_path)
        apply_mdtraj_pbc_correction(traj, anchor_selection="protein", ligand_resname=LIG_RESNAME)
        ca_idx = _nnibp_ca_selection(traj.topology)
        lig_idx = _ligand_heavy_selection(traj.topology)
        prot_idx = traj.topology.select("protein and not element H")

        from scipy.spatial import cKDTree

        for frame in range(0, traj.n_frames, max(1, stride)):
            ca_xyz = traj.xyz[frame, ca_idx, :] * 10.0
            lig_xyz = traj.xyz[frame, lig_idx, :] * 10.0
            rot, mobile_c, target_c = _kabsch(ca_xyz, reference_frame)
            aligned = (lig_xyz - mobile_c) @ rot.T + target_c
            # The closest ligand-protein heavy-atom approach. Frames below roughly 2.2 A
            # are not physically equilibrated structures; boundness_qc only tests the
            # 6 A upper bound, so nothing upstream rejects them.
            min_contact = float(
                cKDTree(traj.xyz[frame, prot_idx, :] * 10.0).query(lig_xyz)[0].min()
            )
            rows.append(
                {
                    "mutation": label,
                    "replicate": replicate,
                    "frame_index": frame,
                    "min_contact_angstrom": min_contact,
                    **pose_coordinates(aligned, reference_ligand),
                }
            )
        logging.info("%s rep %02d: %d frames", label, replicate, traj.n_frames)
    return rows


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-dir", type=Path, default=REPO_ROOT / "results" / "md_runs" / "wt")
    parser.add_argument("--label", default="WT")
    parser.add_argument("--pose-dir", type=Path, default=REPO_ROOT / "data" / "prepared" / "dor_poses")
    parser.add_argument(
        "--prepared-crystal",
        type=Path,
        default=REPO_ROOT / "data" / "prepared" / "dor_poses" / "minimized" / "4NCG_minimized.pdb",
    )
    parser.add_argument("--crystal-cif", type=Path, default=REPO_ROOT / "data" / "structures" / "4NCG.cif")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results" / "analysis" / "pose_coordinates" / "tables",
    )
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference_frame, reference_ligand = _reference(args.prepared_crystal, args.crystal_cif)

    starts = pd.DataFrame(starting_pose_rows(args.pose_dir, args.crystal_cif, reference_ligand))
    starts_path = args.output_dir / "starting_poses.csv"
    starts.to_csv(starts_path, index=False)
    logging.info("Wrote %s (%d poses)", starts_path, len(starts))

    frames = pd.DataFrame(
        trajectory_rows(args.run_dir, args.label, reference_frame, reference_ligand, args.stride)
    )
    frames_path = args.output_dir / "sampled_frames.csv"
    frames.to_csv(frames_path, index=False)
    logging.info("Wrote %s (%d frames)", frames_path, len(frames))

    if not frames.empty:
        for column in ("centroid_displacement_angstrom", "long_axis_tilt_deg", "pose_rmsd_angstrom"):
            series = frames[column]
            logging.info(
                "%s sampled range: %.2f - %.2f (median %.2f, 99th %.2f)",
                column,
                series.min(),
                series.max(),
                series.median(),
                series.quantile(0.99),
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
