from __future__ import annotations

from typing import Sequence

import numpy as np

from .require import require_module


def heavy_atom_indices(topology, exclude_resname: str) -> list[int]:
    indices = []
    for atom in topology.atoms():
        if atom.element is None:
            continue
        if atom.residue.name == exclude_resname:
            continue
        if atom.element.symbol != "H":
            indices.append(atom.index)
    return indices


def restrained_indices(
    positions,
    ligand_indices: Sequence[int],
    candidate_indices: Sequence[int],
    radius_angstrom: float,
) -> list[int]:
    unit = require_module("openmm.unit")
    radius_nm = radius_angstrom / 10.0
    ligand_pos = np.array(
        [positions[i].value_in_unit(unit.nanometer) for i in ligand_indices]
    )
    restrained = []
    for idx in candidate_indices:
        pos = positions[idx].value_in_unit(unit.nanometer)
        d = np.min(np.linalg.norm(ligand_pos - pos, axis=1))
        if d > radius_nm:
            restrained.append(idx)
    return restrained
