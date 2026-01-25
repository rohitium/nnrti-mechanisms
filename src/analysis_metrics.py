from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


@dataclass(frozen=True)
class ContactMetrics:
    contact_count: int
    hbond_count: int | None


def compute_contacts(
    pdbx_path: Path, ligand_resname: str, cutoff_angstrom: float = 4.0
) -> ContactMetrics:
    u = MDAnalysis.Universe(str(pdbx_path))
    ligand = u.select_atoms(f"resname {ligand_resname}")
    protein = u.select_atoms(f"protein and not resname {ligand_resname}")
    if ligand.n_atoms == 0:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {pdbx_path}.")

    distances = MDAnalysis.lib.distances.distance_array(
        ligand.positions, protein.positions
    )
    contact_count = int(np.sum(distances < cutoff_angstrom))

    hbond = HydrogenBondAnalysis(
        u,
        donors_sel=f"resname {ligand_resname} or protein",
        hydrogens_sel="name H*",
        acceptors_sel=f"resname {ligand_resname} or protein",
        between=[f"resname {ligand_resname}", "protein"],
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=135.0,
    )
    hbond.run()
    hbond_count = int(len(hbond.results.hbonds))

    return ContactMetrics(contact_count=contact_count, hbond_count=hbond_count)


def pocket_volume_proxy(
    pdbx_path: Path,
    ligand_resname: str,
    grid_spacing: float = 0.5,
    radius_angstrom: float = 8.0,
) -> float:
    u = MDAnalysis.Universe(str(pdbx_path))
    ligand = u.select_atoms(f"resname {ligand_resname}")
    receptor = u.select_atoms(f"not resname {ligand_resname} and not name H*")
    if ligand.n_atoms == 0:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {pdbx_path}.")

    center = ligand.positions.mean(axis=0)
    radius = radius_angstrom
    spacing = grid_spacing

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
    receptor_elements = []
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
