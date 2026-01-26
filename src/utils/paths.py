from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Paths:
    root: Path
    data: Path
    structures: Path
    ligands: Path
    generated: Path
    results: Path
    plots: Path


def project_paths(root: Path) -> Paths:
    data = root / "data"
    return Paths(
        root=root,
        data=data,
        structures=data / "structures",
        ligands=data / "ligands",
        generated=data / "generated",
        results=root / "results",
        plots=root / "results" / "plots",
    )


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)
