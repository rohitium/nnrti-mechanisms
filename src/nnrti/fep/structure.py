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


# γ/δ/ε methylene hydrogens follow amber's "2/3" convention (HG2/HG3, HD2/HD3,
# HE2/HE3); pmx's hybrid residues expect "1/2" (HG1/HG2, ...). These renames are
# RESIDUE-SCOPED: the same labels denote non-methylene atoms elsewhere (His ring
# HD1/HD2/HE1, Asn/Gln amide HD2/HE2, aromatic HD/HE, and the branched methyls of
# Ile/Leu/Thr/Val), which must not be touched. The residue sets below are exactly
# those whose γ/δ/ε positions are CH2 groups.
_HG_RESIDUES = frozenset({"ARG", "GLN", "GLU", "LYS", "LYP", "LYN", "MET", "PRO"})
_HD_RESIDUES = frozenset({"ARG", "LYS", "LYP", "LYN", "PRO"})
_HE_RESIDUES = frozenset({"LYS", "LYP", "LYN"})
# atom -> (allowed residue names, new name, first-atom-of-pair for idempotency)
_METHYLENE_RENAME = {
    "HG2": (_HG_RESIDUES, " HG1", "HG1"), "HG3": (_HG_RESIDUES, " HG2", "HG1"),
    "HD2": (_HD_RESIDUES, " HD1", "HD1"), "HD3": (_HD_RESIDUES, " HD2", "HD1"),
    "HE2": (_HE_RESIDUES, " HE1", "HE1"), "HE3": (_HE_RESIDUES, " HE2", "HE1"),
}


def normalize_openmm_for_pmx(protein_pdb: Path) -> None:
    """Rename OpenMM/CHARMM atom labels to names pmx amber14sbmut expects.

    Amber names methylene hydrogens with a "2/3" convention (β: HB2/HB3; glycine
    α: HA2/HA3; and the γ/δ/ε methylenes HG2/HG3, HD2/HD3, HE2/HE3); pmx's hybrid
    residues expect "1/2" (HB1/HB2, HA1/HA2, HG1/HG2, ...). Without this, mutating
    a residue that has an unconverted methylene crashes pmx ``_set_conformation``
    (``old_res[name]`` IndexError) while copying the A-state coordinates — seen
    for proline (P225H) and lysine (K103N).

    HB is renamed on every residue that has it; HA2/HA3 occur only on glycine.
    The γ/δ/ε renames are **residue-scoped** (see ``_METHYLENE_RENAME``) because
    HG/HD/HE "2/3" labels denote different atoms in His (ring), Asn/Gln (amide),
    aromatics, and branched-methyl residues — those must be left alone. Every
    rename is **idempotent per-residue**: a residue already carrying the "1"
    hydrogen (HB1/HG1/HD1/HE1) is skipped, so re-normalizing cannot collide
    HB3→HB2→HB1 onto an existing HB1.
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
    hb_ha_rename = {"HB2": " HB1", "HB3": " HB2", "HA2": " HA1", "HA3": " HA2"}
    have_first = {name: _residues_with(name) for name in ("HG1", "HD1", "HE1")}

    normalized: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:26].strip())
            name = line[12:16].strip()
            resname = line[17:20].strip()
            if name in hb_ha_rename:
                already = have_hb1 if name in ("HB2", "HB3") else have_ha1
                if key not in already:
                    line = line[:12] + hb_ha_rename[name] + line[16:]
            elif name in _METHYLENE_RENAME:
                res_ok, new_name, first_atom = _METHYLENE_RENAME[name]
                if resname in res_ok and key not in have_first[first_atom]:
                    line = line[:12] + new_name + line[16:]
        normalized.append(line)
    protein_pdb.write_text("\n".join(normalized + ["END", ""]))
