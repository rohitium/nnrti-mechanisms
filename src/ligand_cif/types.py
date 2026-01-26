from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LigandAtom:
    name: str
    element: str
    x: float
    y: float
    z: float
    aromatic: bool


@dataclass(frozen=True)
class LigandBond:
    atom1: str
    atom2: str
    order: int
    aromatic: bool
