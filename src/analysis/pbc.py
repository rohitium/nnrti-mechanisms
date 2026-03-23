from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class PBCAuditSummary:
    dcd_path: Path
    topology_path: Path
    n_frames: int
    n_atoms: int
    n_molecules: int
    sampled_frames: int
    has_unitcell: bool
    max_bond_length_angstrom: float
    max_anchor_center_offset_angstrom: float
    max_anchor_internal_jump_angstrom: float
    max_ligand_anchor_gap_angstrom: float
    max_ligand_anchor_jump_angstrom: float
    motion_outlier: bool
    passed: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "trajectory_dcd": str(self.dcd_path),
            "topology_pdb": str(self.topology_path),
            "n_frames": int(self.n_frames),
            "n_atoms": int(self.n_atoms),
            "n_molecules": int(self.n_molecules),
            "sampled_frames": int(self.sampled_frames),
            "has_unitcell": bool(self.has_unitcell),
            "max_bond_length_angstrom": float(self.max_bond_length_angstrom),
            "max_anchor_center_offset_angstrom": float(self.max_anchor_center_offset_angstrom),
            "max_anchor_internal_jump_angstrom": float(self.max_anchor_internal_jump_angstrom),
            "max_ligand_anchor_gap_angstrom": float(self.max_ligand_anchor_gap_angstrom),
            "max_ligand_anchor_jump_angstrom": float(self.max_ligand_anchor_jump_angstrom),
            "motion_outlier": bool(self.motion_outlier),
            "passed": bool(self.passed),
        }


def topology_for_analysis_dcd(dcd_path: Path) -> Path:
    name = dcd_path.name
    if not name.endswith("_analysis.dcd"):
        raise ValueError(f"Unexpected DCD name (expected *_analysis.dcd): {dcd_path}")
    topo_name = name.replace("_analysis.dcd", "_analysis_topology.pdb")
    return dcd_path.with_name(topo_name)


def load_mdtraj_trajectory(dcd_path: Path, topo_path: Path):
    import mdtraj as md
    import mdtraj.formats

    if str(dcd_path).endswith(".bak"):
        ref = md.load(str(topo_path))
        with mdtraj.formats.DCDTrajectoryFile(str(dcd_path), "r") as handle:
            xyz, lengths, angles = handle.read()
        return md.Trajectory(
            xyz / 10.0,
            ref.topology,
            unitcell_lengths=(lengths / 10.0) if lengths is not None else None,
            unitcell_angles=angles,
        )
    return md.load(str(dcd_path), top=str(topo_path))


def _atom_indices_from_molecule(molecule) -> np.ndarray:
    return np.fromiter((atom.index for atom in molecule), dtype=int)


def _safe_box_lengths(traj, frame_i: int) -> np.ndarray | None:
    if getattr(traj, "unitcell_lengths", None) is None:
        return None
    box = traj.unitcell_lengths[frame_i]
    if box is None or not np.all(np.isfinite(box)) or not np.all(box > 0):
        return None
    return np.asarray(box, dtype=float)


def _anchor_molecules(topology, anchor_selection: str = "protein") -> list[set]:
    anchor_atom_idx = set(int(i) for i in topology.select(str(anchor_selection)))
    molecules = list(topology.find_molecules())
    anchors = [mol for mol in molecules if any(atom.index in anchor_atom_idx for atom in mol)]
    if anchors:
        return anchors
    if molecules:
        return [max(molecules, key=len)]
    return []


def _anchor_atom_indices(topology, anchor_selection: str = "protein") -> np.ndarray:
    anchor_idx = topology.select(str(anchor_selection))
    if len(anchor_idx) > 0:
        return np.asarray(anchor_idx, dtype=int)
    return np.arange(topology.n_atoms, dtype=int)


def _wrap_selection_to_anchor(traj, selection_idx: np.ndarray, anchor_idx: np.ndarray) -> None:
    if selection_idx.size == 0 or anchor_idx.size == 0:
        return
    for frame_i in range(traj.n_frames):
        box = _safe_box_lengths(traj, frame_i)
        if box is None:
            continue
        anchor_com = traj.xyz[frame_i, anchor_idx].mean(axis=0)
        sel_com = traj.xyz[frame_i, selection_idx].mean(axis=0)
        traj.xyz[frame_i, selection_idx] += -box * np.round((sel_com - anchor_com) / box)


def apply_mdtraj_pbc_correction(
    traj,
    *,
    anchor_selection: str = "protein",
    ligand_resname: str | None = "2KW",
) -> None:
    if traj.n_atoms < 1:
        raise ValueError("Cannot PBC-correct an empty trajectory.")

    topology = traj.topology
    traj.make_molecules_whole(inplace=True)

    molecules = list(topology.find_molecules())
    molecule_indices = [_atom_indices_from_molecule(mol) for mol in molecules]

    for frame_i in range(1, traj.n_frames):
        box = _safe_box_lengths(traj, frame_i)
        if box is None:
            continue
        for atom_idx in molecule_indices:
            if atom_idx.size == 0:
                continue
            prev_com = traj.xyz[frame_i - 1, atom_idx].mean(axis=0)
            curr_com = traj.xyz[frame_i, atom_idx].mean(axis=0)
            traj.xyz[frame_i, atom_idx] += -box * np.round((curr_com - prev_com) / box)

    anchors = _anchor_molecules(topology, anchor_selection=anchor_selection)
    anchor_idx = _anchor_atom_indices(topology, anchor_selection=anchor_selection)

    if anchors and getattr(traj, "unitcell_lengths", None) is not None:
        traj.image_molecules(inplace=True, anchor_molecules=anchors, make_whole=False)

    for frame_i in range(traj.n_frames):
        box = _safe_box_lengths(traj, frame_i)
        if box is None:
            continue
        anchor_com = traj.xyz[frame_i, anchor_idx].mean(axis=0)
        traj.xyz[frame_i, :, :] += (0.5 * box) - anchor_com

    if ligand_resname:
        ligand_idx = np.asarray(topology.select(f"resname '{ligand_resname}'"), dtype=int)
        _wrap_selection_to_anchor(traj, ligand_idx, anchor_idx)


def audit_mdtraj_trajectory(
    traj,
    *,
    dcd_path: Path,
    topo_path: Path,
    anchor_selection: str = "protein",
    ligand_resname: str | None = "2KW",
    frame_stride: int = 10,
    max_bond_length_angstrom: float = 3.0,
    max_anchor_center_offset_angstrom: float = 3.0,
    max_ligand_anchor_gap_angstrom: float = 1.0,
    max_anchor_internal_jump_angstrom: float = 15.0,
    max_ligand_anchor_jump_angstrom: float = 8.0,
) -> PBCAuditSummary:
    topology = traj.topology
    frame_idx = np.arange(0, traj.n_frames, max(1, int(frame_stride)), dtype=int)
    if frame_idx.size == 0 and traj.n_frames > 0:
        frame_idx = np.array([0], dtype=int)

    anchor_idx = _anchor_atom_indices(topology, anchor_selection=anchor_selection)
    ligand_idx = (
        np.asarray(topology.select(f"resname '{ligand_resname}'"), dtype=int)
        if ligand_resname
        else np.array([], dtype=int)
    )

    molecules = list(topology.find_molecules())
    anchor_molecules = _anchor_molecules(topology, anchor_selection=anchor_selection)
    anchor_mol_idx = [_atom_indices_from_molecule(mol) for mol in anchor_molecules]

    bond_pairs: np.ndarray
    bonds = [(bond.atom1.index, bond.atom2.index) for bond in topology.bonds]
    if bonds:
        bond_pairs = np.asarray(bonds, dtype=int)
    else:
        bond_pairs = np.empty((0, 2), dtype=int)

    max_bond = 0.0
    max_center_offset = 0.0
    max_anchor_jump = 0.0
    max_gap = 0.0
    max_ligand_jump = 0.0
    prev_anchor_rel = None
    prev_ligand_rel = None

    for fi in frame_idx:
        xyz = traj.xyz[int(fi)]
        box = _safe_box_lengths(traj, int(fi))
        anchor_com = xyz[anchor_idx].mean(axis=0)

        if bond_pairs.size > 0:
            delta = xyz[bond_pairs[:, 0]] - xyz[bond_pairs[:, 1]]
            if box is not None:
                delta = delta - box * np.round(delta / box)
            dist = np.linalg.norm(delta, axis=1) * 10.0
            if dist.size > 0:
                max_bond = max(max_bond, float(np.nanmax(dist)))

        if box is not None:
            max_center_offset = max(
                max_center_offset,
                float(np.linalg.norm(anchor_com - (0.5 * box)) * 10.0),
            )

        if anchor_mol_idx:
            anchor_rel = np.vstack([xyz[idx].mean(axis=0) - anchor_com for idx in anchor_mol_idx])
            if prev_anchor_rel is not None and prev_anchor_rel.shape == anchor_rel.shape:
                jump = np.linalg.norm(anchor_rel - prev_anchor_rel, axis=1) * 10.0
                if jump.size > 0:
                    max_anchor_jump = max(max_anchor_jump, float(np.nanmax(jump)))
            prev_anchor_rel = anchor_rel

        if ligand_idx.size > 0:
            ligand_com = xyz[ligand_idx].mean(axis=0)
            delta = ligand_com - anchor_com
            direct_dist = float(np.linalg.norm(delta) * 10.0)
            if box is not None:
                min_image = delta - box * np.round(delta / box)
                min_image_dist = float(np.linalg.norm(min_image) * 10.0)
            else:
                min_image_dist = direct_dist
            max_gap = max(max_gap, abs(direct_dist - min_image_dist))

            ligand_rel = delta
            if prev_ligand_rel is not None:
                jump_delta = ligand_rel - prev_ligand_rel
                if box is not None:
                    jump_delta = jump_delta - box * np.round(jump_delta / box)
                max_ligand_jump = max(max_ligand_jump, float(np.linalg.norm(jump_delta) * 10.0))
            prev_ligand_rel = ligand_rel

    hard_passed = (
        max_bond <= float(max_bond_length_angstrom)
        and max_center_offset <= float(max_anchor_center_offset_angstrom)
        and max_gap <= float(max_ligand_anchor_gap_angstrom)
    )
    motion_outlier = (
        max_anchor_jump > float(max_anchor_internal_jump_angstrom)
        or max_ligand_jump > float(max_ligand_anchor_jump_angstrom)
    )

    return PBCAuditSummary(
        dcd_path=dcd_path,
        topology_path=topo_path,
        n_frames=int(traj.n_frames),
        n_atoms=int(traj.n_atoms),
        n_molecules=int(len(molecules)),
        sampled_frames=int(frame_idx.size),
        has_unitcell=getattr(traj, "unitcell_lengths", None) is not None,
        max_bond_length_angstrom=float(max_bond),
        max_anchor_center_offset_angstrom=float(max_center_offset),
        max_anchor_internal_jump_angstrom=float(max_anchor_jump),
        max_ligand_anchor_gap_angstrom=float(max_gap),
        max_ligand_anchor_jump_angstrom=float(max_ligand_jump),
        motion_outlier=bool(motion_outlier),
        passed=bool(hard_passed),
    )
