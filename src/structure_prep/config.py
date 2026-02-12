from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class StructureSpec:
    name: str
    cif_path: Path
    ligand_resname: str
    ligand_chain: str
    protein_chains: tuple[str, ...]
    dna_chains: tuple[str, ...]
    ligand_sdf: Path


@dataclass(frozen=True)
class RunSpec:
    structure: StructureSpec
    restraint_radius_angstrom: float = 8.0
    restraint_k_kj_mol_nm2: float = 500.0


def dor_4ncg_spec(root: Path) -> RunSpec:
    return RunSpec(
        structure=StructureSpec(
            name="DOR",
            cif_path=root / "data" / "structures" / "4NCG.cif",
            ligand_resname="2KW",
            ligand_chain="C",
            protein_chains=("A", "B"),
            dna_chains=(),
            ligand_sdf=root / "data" / "ligands" / "dor.sdf",
        ),
    )
