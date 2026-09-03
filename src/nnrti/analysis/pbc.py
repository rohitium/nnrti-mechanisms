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


def pbcfix_dcd_for(analysis_dcd: Path) -> Path:
    name = analysis_dcd.name
    if name.endswith("_analysis_pbcfix.dcd"):
        return analysis_dcd
    if name.endswith("_analysis.dcd"):
        return analysis_dcd.with_name(name.replace("_analysis.dcd", "_analysis_pbcfix.dcd"))
    return analysis_dcd.with_name(f"{analysis_dcd.stem}_pbcfix.dcd")


def raw_analysis_dcd_for(dcd_path: Path) -> Path:
    name = dcd_path.name
    if name.endswith("_analysis_pbcfix.dcd"):
        return dcd_path.with_name(name.replace("_analysis_pbcfix.dcd", "_analysis.dcd"))
    return dcd_path


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


def topology_for_analysis_dcd(dcd_path: Path) -> Path:
    name = dcd_path.name
    if name.endswith("_analysis_pbcfix.dcd"):
        return dcd_path.with_name(name.replace("_analysis_pbcfix.dcd", "_analysis_topology.pdb"))
    if name.endswith("_analysis.dcd"):
        return dcd_path.with_name(name.replace("_analysis.dcd", "_analysis_topology.pdb"))
    raise ValueError(f"Unexpected DCD name (expected *_analysis.dcd or *_analysis_pbcfix.dcd): {dcd_path}")

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
