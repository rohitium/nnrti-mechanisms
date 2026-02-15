from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ContactMetrics:
    contact_count: int
    hbond_count: int | None


@dataclass(frozen=True)
class EnsembleMetrics:
    contact_count_mean: float
    contact_count_std: float
    hbond_count_mean: float
    hbond_count_std: float
    pocket_volume_proxy_mean: float
    pocket_volume_proxy_std: float
    n_frames: int


def compute_contacts(
    pdbx_path: Path, ligand_resname: str, cutoff_angstrom: float = 4.0
) -> ContactMetrics:
    import MDAnalysis
    from MDAnalysis.lib.distances import capped_distance
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    u = MDAnalysis.Universe(str(pdbx_path))
    ligand = u.select_atoms(f"resname {ligand_resname} and not name H*")
    protein = u.select_atoms(f"protein and not resname {ligand_resname} and not name H*")
    if ligand.n_atoms == 0:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {pdbx_path}.")
    if protein.n_atoms == 0:
        raise ValueError("Protein selection is empty.")

    pairs = capped_distance(
        ligand.positions,
        protein.positions,
        max_cutoff=cutoff_angstrom,
        box=u.dimensions,
        return_distances=False,
    )
    if len(pairs) == 0:
        contact_count = 0
    else:
        protein_idx = np.asarray(pairs)[:, 1]
        contact_resids = {int(protein.atoms[i].resindex) for i in protein_idx}
        contact_count = int(len(contact_resids))

    hbond = HydrogenBondAnalysis(
        u,
        donors_sel=f"resname {ligand_resname} or protein",
        hydrogens_sel="name H*",
        acceptors_sel=f"resname {ligand_resname} or protein",
        between=[f"resname {ligand_resname}", "protein"],
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=135.0,
    )
    hbond.run(start=0, stop=1, step=1)
    hb = hbond.results.hbonds
    if hb is None or len(hb) == 0:
        hbond_count = 0
    else:
        unique_pairs = {(int(row[1]), int(row[3])) for row in hb}
        hbond_count = int(len(unique_pairs))

    return ContactMetrics(contact_count=contact_count, hbond_count=hbond_count)


def _pocket_volume_proxy_universe(
    universe,
    ligand_resname: str,
    grid_spacing: float = 0.5,
    radius_angstrom: float = 8.0,
) -> float:
    ligand = universe.select_atoms(f"resname {ligand_resname}")
    receptor = universe.select_atoms(f"not resname {ligand_resname} and not name H*")
    if ligand.n_atoms == 0:
        raise ValueError("Ligand not found for pocket-volume calculation.")

    center = ligand.positions.mean(axis=0)
    spacing = grid_spacing
    radius = radius_angstrom

    mins = center - radius
    maxs = center + radius
    xs = np.arange(mins[0], maxs[0] + spacing, spacing)
    ys = np.arange(mins[1], maxs[1] + spacing, spacing)
    zs = np.arange(mins[2], maxs[2] + spacing, spacing)

    grid = np.array(np.meshgrid(xs, ys, zs, indexing="ij")).reshape(3, -1).T
    d_center = np.linalg.norm(grid - center, axis=1)
    grid = grid[d_center <= radius]

    vdw = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8}
    receptor_pos = receptor.positions
    receptor_elements: list[str] = []
    for atom in receptor.atoms:
        elem = atom.element
        if elem is None:
            elem = atom.name[0].upper()
        receptor_elements.append(elem)

    free_mask = np.ones(len(grid), dtype=bool)
    for pos, elem in zip(receptor_pos, receptor_elements):
        rad = vdw.get(elem, 1.7)
        d = np.linalg.norm(grid - pos, axis=1)
        free_mask &= d > rad

    voxel_volume = spacing**3
    return float(np.sum(free_mask) * voxel_volume)


def pocket_volume_proxy(
    pdbx_path: Path,
    ligand_resname: str,
    grid_spacing: float = 0.5,
    radius_angstrom: float = 8.0,
) -> float:
    import MDAnalysis

    u = MDAnalysis.Universe(str(pdbx_path))
    return _pocket_volume_proxy_universe(
        u,
        ligand_resname=ligand_resname,
        grid_spacing=grid_spacing,
        radius_angstrom=radius_angstrom,
    )


def compute_ensemble_metrics(
    topology_pdb_path: Path,
    trajectory_dcd_path: Path,
    ligand_resname: str,
    frame_stride: int = 1,
    max_frames: int | None = 200,
    sample_window_ns: float | None = 1.0,
    contact_cutoff_angstrom: float = 4.0,
    grid_spacing: float = 0.75,
    pocket_radius_angstrom: float = 8.0,
) -> EnsembleMetrics:
    """Compute ensemble-averaged ligand/protein metrics from a trajectory."""
    import MDAnalysis
    from MDAnalysis.lib.distances import capped_distance
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    u = MDAnalysis.Universe(str(topology_pdb_path), str(trajectory_dcd_path))
    ligand = u.select_atoms(f"resname {ligand_resname} and not name H*")
    protein = u.select_atoms(f"protein and not resname {ligand_resname} and not name H*")
    if ligand.n_atoms == 0:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {topology_pdb_path}.")
    if protein.n_atoms == 0:
        raise ValueError("Protein selection is empty for ensemble metric calculation.")

    frame_step = max(1, frame_stride)
    sampled_frames = list(range(0, len(u.trajectory), frame_step))
    if not sampled_frames:
        raise ValueError("No frames available for ensemble metric calculation.")

    frame_indices: list[int]
    if sample_window_ns is not None and sample_window_ns > 0.0 and len(sampled_frames) > 1:
        dt_ps = getattr(u.trajectory, "dt", None)
        corrected_dt_ps = dt_ps
        if corrected_dt_ps is not None and np.isfinite(corrected_dt_ps) and corrected_dt_ps > 1000.0:
            candidate = float(corrected_dt_ps) / 1000.0
            corrected_dt_ps = candidate if candidate <= 1000.0 else None
        if corrected_dt_ps is not None and np.isfinite(corrected_dt_ps) and corrected_dt_ps > 0:
            total_time_ps = (len(u.trajectory) - 1) * float(corrected_dt_ps)
            window_ps = float(sample_window_ns) * 1000.0
            start_time_ps = max(0.0, total_time_ps - window_ps)
            sampled_in_window = [
                frame_id
                for frame_id in sampled_frames
                if (frame_id * float(corrected_dt_ps)) >= start_time_ps
            ]
            frame_indices = sampled_in_window if sampled_in_window else sampled_frames
        else:
            frame_indices = sampled_frames
    else:
        frame_indices = sampled_frames

    if max_frames is not None and len(frame_indices) > max_frames:
        select_idx = np.linspace(0, len(frame_indices) - 1, num=max_frames, dtype=int)
        frame_indices = [frame_indices[i] for i in select_idx.tolist()]

    contact_values: list[float] = []
    pocket_values: list[float] = []
    for frame_id in frame_indices:
        u.trajectory[frame_id]
        pairs = capped_distance(
            ligand.positions,
            protein.positions,
            max_cutoff=contact_cutoff_angstrom,
            box=u.dimensions,
            return_distances=False,
        )
        if len(pairs) == 0:
            contact_values.append(0.0)
        else:
            protein_idx = np.asarray(pairs)[:, 1]
            contact_resids = {int(protein.atoms[i].resindex) for i in protein_idx}
            contact_values.append(float(len(contact_resids)))
        pocket_values.append(
            _pocket_volume_proxy_universe(
                u,
                ligand_resname=ligand_resname,
                grid_spacing=grid_spacing,
                radius_angstrom=pocket_radius_angstrom,
            )
        )

    if not frame_indices:
        raise ValueError("No frames available for ensemble metric calculation.")

    frame_set = set(frame_indices)
    hbond = HydrogenBondAnalysis(
        u,
        donors_sel=f"resname {ligand_resname} or protein",
        hydrogens_sel="name H*",
        acceptors_sel=f"resname {ligand_resname} or protein",
        between=[f"resname {ligand_resname}", "protein"],
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=135.0,
    )
    hbond.run(
        start=min(frame_indices),
        stop=max(frame_indices) + 1,
        step=frame_step,
    )
    frame_to_hbond: dict[int, int] = {idx: 0 for idx in frame_indices}
    hb = hbond.results.hbonds
    if hb is not None and len(hb) > 0:
        frame_pair_counts: dict[int, set[tuple[int, int]]] = {
            idx: set() for idx in frame_indices
        }
        for row in hb:
            frame_id = int(row[0])
            if frame_id in frame_set:
                frame_pair_counts[frame_id].add((int(row[1]), int(row[3])))
        for frame_id, pair_set in frame_pair_counts.items():
            frame_to_hbond[frame_id] = len(pair_set)
    hbond_values = [float(frame_to_hbond[idx]) for idx in frame_indices]

    return EnsembleMetrics(
        contact_count_mean=float(np.mean(contact_values)),
        contact_count_std=float(np.std(contact_values)),
        hbond_count_mean=float(np.mean(hbond_values)),
        hbond_count_std=float(np.std(hbond_values)),
        pocket_volume_proxy_mean=float(np.mean(pocket_values)),
        pocket_volume_proxy_std=float(np.std(pocket_values)),
        n_frames=len(frame_indices),
    )
