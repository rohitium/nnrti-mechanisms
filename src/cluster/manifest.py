from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class FEPTask:
    """Represents a single FEP computation task for cluster execution.

    The cluster worker will:
    1. Minimize the structure from input_cif (with jitter)
    2. Run the specified FEP leg (complex or solvent)
    3. Save results to output_json
    """

    task_id: int
    structure: str
    mutation: str
    safe_label: str
    replicate: int
    leg: str  # "complex" or "solvent"
    minimized_pdb: str  # Output path for minimized structure
    ligand_sdf: str
    ligand_resname: str
    fold_reduction: float | None
    output_json: str
    # Fields for cluster minimization
    input_cif: str = ""  # Input CIF to minimize
    jitter_seed: int | None = None
    jitter_angstrom: float = 0.1
    restraint_radius: float = 8.0
    restraint_k: float = 500.0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> FEPTask:
        fold = data.get("fold_reduction")
        if fold is not None and pd.isna(fold):
            fold = None

        jitter_seed = data.get("jitter_seed")
        if jitter_seed is not None and pd.isna(jitter_seed):
            jitter_seed = None
        elif jitter_seed is not None:
            jitter_seed = int(jitter_seed)

        return cls(
            task_id=int(data["task_id"]),
            structure=str(data["structure"]),
            mutation=str(data["mutation"]),
            safe_label=str(data["safe_label"]),
            replicate=int(data["replicate"]),
            leg=str(data["leg"]),
            minimized_pdb=str(data.get("minimized_pdb", "")),
            ligand_sdf=str(data["ligand_sdf"]),
            ligand_resname=str(data["ligand_resname"]),
            fold_reduction=float(fold) if fold is not None else None,
            output_json=str(data["output_json"]),
            input_cif=str(data.get("input_cif", "")),
            jitter_seed=jitter_seed,
            jitter_angstrom=float(data.get("jitter_angstrom", 0.1)),
            restraint_radius=float(data.get("restraint_radius", 8.0)),
            restraint_k=float(data.get("restraint_k", 500.0)),
        )


def save_manifest(tasks: list[FEPTask], output_path: Path) -> None:
    """Save a list of FEPTasks to a CSV manifest file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([t.to_dict() for t in tasks])
    df.to_csv(output_path, index=False)


def load_manifest(manifest_path: Path) -> list[FEPTask]:
    """Load FEPTasks from a CSV manifest file."""
    df = pd.read_csv(manifest_path)
    return [FEPTask.from_dict(row) for _, row in df.iterrows()]


def get_task_by_id(manifest_path: Path, task_id: int) -> FEPTask:
    """Load a single FEPTask by its task_id from the manifest."""
    df = pd.read_csv(manifest_path)
    row = df[df["task_id"] == task_id]
    if row.empty:
        raise ValueError(f"Task ID {task_id} not found in manifest {manifest_path}")
    return FEPTask.from_dict(row.iloc[0])
