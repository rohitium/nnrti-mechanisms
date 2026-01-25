from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MutationSpec:
    chain_id: str
    residue_id: str
    new_residue: str
    label: str


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
    mutation: MutationSpec
    restraint_radius_angstrom: float = 8.0
    restraint_k_kj_mol_nm2: float = 500.0


def rpv_spec(root: Path) -> RunSpec:
    return RunSpec(
        structure=StructureSpec(
            name="RPV",
            cif_path=root / "data" / "structures" / "7Z2D.cif",
            ligand_resname="T27",
            ligand_chain="D",
            protein_chains=("A", "B"),
            dna_chains=("C",),
            ligand_sdf=root / "data" / "ligands" / "rpv.sdf",
        ),
        mutation=MutationSpec(
            chain_id="A", residue_id="138", new_residue="LYS", label="E138K"
        ),
    )


def dor_spec(root: Path) -> RunSpec:
    return RunSpec(
        structure=StructureSpec(
            name="DOR",
            cif_path=root / "data" / "structures" / "7Z2G.cif",
            ligand_resname="2KW",
            ligand_chain="D",
            protein_chains=("A", "B"),
            dna_chains=("C",),
            ligand_sdf=root / "data" / "ligands" / "dor.sdf",
        ),
        mutation=MutationSpec(
            chain_id="A", residue_id="106", new_residue="ALA", label="V106A"
        ),
    )
