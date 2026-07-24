"""Structure extraction helpers for Perses FEP inputs."""

from __future__ import annotations

from pathlib import Path

SOLVENT_RESNAMES = {"HOH", "WAT", "SOL"}
ION_RESNAMES = {"NA", "CL", "K", "Na+", "Cl-"}


def extract_protein_only(
    source_pdb: Path,
    output_dir: Path,
    *,
    ligand_resname: str = "2KW",
    output_name: str = "protein_no_ligand.pdb",
) -> Path:
    """Write an unsolvated protein PDB from a solvated holo or apo MD structure."""
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_pdb = output_dir / output_name
    protein_lines: list[str] = []
    for line in source_pdb.read_text().splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue
        resname = line[17:20].strip()
        if resname == ligand_resname or resname in SOLVENT_RESNAMES or resname in ION_RESNAMES:
            continue
        protein_lines.append(line)
    if not protein_lines:
        raise ValueError(f"No protein atoms found in {source_pdb}")
    protein_pdb.write_text("\n".join(protein_lines + ["END", ""]))
    return protein_pdb
