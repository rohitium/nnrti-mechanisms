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


# NNBP-lining residue topology resids (canonical HIV-1 RT numbering − 3, resid_offset=−3).
# p66 subunit (15 residues): L100, K101, K103, V106, T107, V108, V179, Y181, Y188,
#                             V189, G190, F227, W229, L234, Y318
# p51 subunit (1 residue):   E138
# Source: Cilento, Kirby & Sarafianos, "Avoiding Drug Resistance in HIV Reverse Transcriptase"
_NNBP_CA_RESIDS_P66: tuple[int, ...] = (
    97, 98, 100, 103, 104, 105,   # L100, K101, K103, V106, T107, V108
    176, 178, 185, 186, 187,      # V179, Y181, Y188, V189, G190
    224, 226, 231, 315,           # F227, W229, L234, Y318
)
_NNBP_CA_RESID_P51: int = 135   # E138, topology resid = 138 − 3
_VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8}
_PROBE_RADIUS = 1.4  # Å — solvent probe


def _p66_segid(universe) -> str:
    """Return the segid of the p66 subunit (largest protein segment = p66 in HIV-1 RT)."""
    from collections import Counter
    prot_ca = universe.select_atoms("protein and name CA")
    if prot_ca.n_atoms == 0:
        return ""
    cnt = Counter(prot_ca.segids.tolist())
    return max(cnt, key=cnt.get)


def _p51_segid(universe) -> str:
    """Return the segid of the p51 subunit (second largest protein segment in HIV-1 RT)."""
    from collections import Counter
    prot_ca = universe.select_atoms("protein and name CA")
    if prot_ca.n_atoms == 0:
        return ""
    cnt = Counter(prot_ca.segids.tolist())
    sorted_segs = sorted(cnt, key=cnt.get, reverse=True)
    return sorted_segs[1] if len(sorted_segs) > 1 else ""


def _pocket_volume_proxy_universe(
    universe,
    ligand_resname: str = "",  # kept for API compat, no longer used
    grid_spacing: float = 0.75,
    radius_angstrom: float = 10.0,
) -> float:
    """Compute NNBP pocket volume using NNBP residue Cα centroid as center.

    Works for both apo and holo trajectories (does not require a drug).
    Center = per-frame centroid of Cα atoms from NNBP-lining residues:
      p66 (15 residues): L100, K101, K103, V106, T107, V108, V179, Y181, Y188,
                         V189, G190, F227, W229, L234, Y318
      p51 (1 residue):   E138
    Source: Cilento, Kirby & Sarafianos, "Avoiding Drug Resistance in HIV RT"
    """
    from scipy.spatial import cKDTree

    # p66 Cα selection
    seg = _p66_segid(universe)
    seg_filter = f" and segid {seg}" if seg else ""
    ca_sel_str = ("protein and name CA" + seg_filter + " and (" +
                  " or ".join(f"resid {r}" for r in _NNBP_CA_RESIDS_P66) + ")")
    ca_sel = universe.select_atoms(ca_sel_str)

    # p51 Cα selection (E138)
    p51_seg = _p51_segid(universe)
    ca_p51_pos = np.empty((0, 3), dtype=float)
    if p51_seg:
        p51_ca = universe.select_atoms(
            f"protein and name CA and segid {p51_seg} and resid {_NNBP_CA_RESID_P51}"
        )
        if p51_ca.n_atoms > 0:
            ca_p51_pos = p51_ca.positions

    total_ca = ca_sel.n_atoms + ca_p51_pos.shape[0]
    if total_ca < 3:
        raise ValueError(f"Too few NNBP Cα found ({total_ca}); check resid offset.")

    ca_positions = (np.vstack([ca_sel.positions, ca_p51_pos])
                    if ca_p51_pos.shape[0] > 0 else ca_sel.positions)
    receptor = universe.select_atoms("protein and not name H*")

    center = ca_positions.mean(axis=0)
    s = grid_spacing
    r = radius_angstrom
    axes = [np.arange(center[i] - r, center[i] + r + s, s) for i in range(3)]
    grid = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    grid = grid[np.linalg.norm(grid - center, axis=1) <= r]
    if len(grid) == 0:
        return 0.0

    receptor_pos = receptor.positions
    rec_vdw = np.array([
        _VDW.get((a.element or a.name[0]).upper()[:1], 1.7)
        for a in receptor.atoms
    ], dtype=float)

    tree = cKDTree(receptor_pos)
    max_excl = float(rec_vdw.max()) + _PROBE_RADIUS
    candidate_idx = tree.query_ball_point(grid, r=max_excl)

    free_mask = np.ones(len(grid), dtype=bool)
    for i, neighbours in enumerate(candidate_idx):
        if not neighbours:
            continue
        nb = np.array(neighbours)
        d = np.linalg.norm(grid[i] - receptor_pos[nb], axis=1)
        if np.any(d < rec_vdw[nb] + _PROBE_RADIUS):
            free_mask[i] = False

    return float(np.sum(free_mask) * s**3)


def pocket_volume_proxy_from_universe(
    universe,
    ligand_resname: str = "",
    grid_spacing: float = 0.75,
    radius_angstrom: float = 10.0,
) -> float:
    """Compute NNBP pocket volume for the current frame of an MDAnalysis Universe."""
    return _pocket_volume_proxy_universe(
        universe,
        grid_spacing=grid_spacing,
        radius_angstrom=radius_angstrom,
    )


def pocket_volume_proxy(
    pdbx_path: Path,
    ligand_resname: str = "",
    grid_spacing: float = 0.75,
    radius_angstrom: float = 10.0,
) -> float:
    import MDAnalysis

    u = MDAnalysis.Universe(str(pdbx_path))
    return _pocket_volume_proxy_universe(
        u,
        grid_spacing=grid_spacing,
        radius_angstrom=radius_angstrom,
    )


def compute_ensemble_metrics(
    topology_pdb_path: Path,
    trajectory_dcd_path: Path,
    ligand_resname: str,
    frame_stride: int = 1,
    max_frames: int | None = 200,
    sample_window_ns: float | None = None,
    total_time_ns: float | None = None,
    contact_cutoff_angstrom: float = 4.0,
    grid_spacing: float = 0.75,
    pocket_radius_angstrom: float = 10.0,
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
        # Prefer a frame-fraction mapping when the total simulated time is known, because
        # DCD dt metadata is not always present or reliable across all trajectories.
        if total_time_ns is not None and np.isfinite(total_time_ns) and float(total_time_ns) > 0.0:
            total_ns = float(total_time_ns)
            window_ns = float(sample_window_ns)
            if window_ns >= total_ns:
                start_frame = 0
            else:
                # Map time->frame assuming frames are uniformly spaced over production.
                frac_start = max(0.0, (total_ns - window_ns) / total_ns)
                start_frame = int(np.ceil(frac_start * float(len(u.trajectory) - 1)))
            frame_indices = [frame_id for frame_id in sampled_frames if frame_id >= start_frame]
            if not frame_indices:
                frame_indices = sampled_frames
        else:
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
