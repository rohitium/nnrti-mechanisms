from __future__ import annotations

import csv
from dataclasses import dataclass, asdict, fields
from pathlib import Path


@dataclass(frozen=True)
class MDTask:
    """Represents a single MD simulation task for cluster execution.

    The cluster worker will:
    1. Minimize the structure from input_cif (with jitter)
    2. Run explicit-solvent MD
    3. Save results to output_json
    """

    task_id: int
    structure: str
    mutation: str
    safe_label: str
    replicate: int
    minimized_pdb: str  # Output path for minimized structure
    ligand_sdf: str
    ligand_resname: str
    fold_reduction: float | None
    output_json: str
    leg: str = "complex"
    # Fields for cluster minimization
    input_cif: str = ""  # Input CIF to minimize
    jitter_seed: int | None = None
    jitter_angstrom: float = 0.1
    restraint_radius: float = 8.0
    restraint_k: float = 500.0
    # Optional OpenMM-only execution assets (prepared locally)
    prepared_topology_pdb: str = ""
    prepared_system_xml: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> MDTask:
        def _norm(value):
            if value is None:
                return None
            text = str(value).strip()
            if not text or text.lower() in {"nan", "none"}:
                return None
            return text

        fold_txt = _norm(data.get("fold_reduction"))
        fold = float(fold_txt) if fold_txt is not None else None

        jitter_txt = _norm(data.get("jitter_seed"))
        jitter_seed = int(jitter_txt) if jitter_txt is not None else None

        return cls(
            task_id=int(data["task_id"]),
            structure=str(data["structure"]),
            mutation=str(data["mutation"]),
            safe_label=str(data["safe_label"]),
            replicate=int(data["replicate"]),
            minimized_pdb=str(data.get("minimized_pdb", "")),
            ligand_sdf=str(data["ligand_sdf"]),
            ligand_resname=str(data["ligand_resname"]),
            fold_reduction=float(fold) if fold is not None else None,
            output_json=str(data["output_json"]),
            leg=str(data.get("leg", "complex")),
            input_cif=str(data.get("input_cif", "")),
            jitter_seed=jitter_seed,
            jitter_angstrom=float(data.get("jitter_angstrom", 0.1)),
            restraint_radius=float(data.get("restraint_radius", 8.0)),
            restraint_k=float(data.get("restraint_k", 500.0)),
            prepared_topology_pdb=str(data.get("prepared_topology_pdb", "")),
            prepared_system_xml=str(data.get("prepared_system_xml", "")),
        )


# Backward-compatibility alias.
FEPTask = MDTask


def save_manifest(tasks: list[MDTask], output_path: Path) -> None:
    """Save a list of MDTasks to a CSV manifest file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [f.name for f in fields(MDTask)]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for task in tasks:
            writer.writerow(task.to_dict())


def load_manifest(manifest_path: Path) -> list[MDTask]:
    """Load MDTasks from a CSV manifest file."""
    with manifest_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        return [MDTask.from_dict(row) for row in reader]


def get_task_by_id(manifest_path: Path, task_id: int) -> MDTask:
    """Load a single MDTask by its task_id from the manifest."""
    for task in load_manifest(manifest_path):
        if task.task_id == task_id:
            return task
    raise ValueError(f"Task ID {task_id} not found in manifest {manifest_path}")
