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


def normalize_openmm_for_pmx(protein_pdb: Path) -> None:
    """Rename OpenMM/CHARMM atom labels to names pmx amber14sbmut expects.

    Amber names the second/third methylene hydrogens HB2/HB3 (β) and HA2/HA3
    (glycine α); pmx wants HB1/HB2 and HA1/HA2. The rename is **idempotent** and
    **per-residue**: a residue that already carries HB1 (or HA1) is skipped, so
    re-normalizing an already-converted structure cannot collide HB3→HB2→HB1
    onto an existing HB1 (which produced duplicate HB1/HB1 atoms pmx could not
    resolve). HA2/HA3 only occur on glycine, so that rename is glycine-specific.
    """
    lines = protein_pdb.read_text().splitlines()

    def _residues_with(atom: str) -> set[tuple[str, str]]:
        return {
            (line[21], line[22:26].strip())
            for line in lines
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == atom
        }

    have_hb1 = _residues_with("HB1")
    have_ha1 = _residues_with("HA1")
    rename = {"HB2": " HB1", "HB3": " HB2", "HA2": " HA1", "HA3": " HA2"}

    normalized: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:26].strip())
            name = line[12:16].strip()
            already = have_hb1 if name in ("HB2", "HB3") else have_ha1
            if name in rename and key not in already:
                line = line[:12] + rename[name] + line[16:]
        normalized.append(line)
    protein_pdb.write_text("\n".join(normalized + ["END", ""]))
