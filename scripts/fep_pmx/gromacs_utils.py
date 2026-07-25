"""GROMACS helpers for pmx NEQ system builds."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


class GromacsError(RuntimeError):
    pass


def find_gmx() -> str:
    gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
    if gmx is None:
        raise GromacsError(
            "gmx not found on PATH. On Sherlock: source scripts/sherlock/load_gromacs_module.sh"
        )
    return gmx


def run_gmx(
    gmx: str,
    args: list[str],
    *,
    cwd: Path,
    input_text: str = "",
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    cmd = [gmx, *args]
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if check and proc.returncode != 0:
        raise GromacsError(
            f"gmx {' '.join(args)} failed ({proc.returncode})\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


def parse_gro_atom_count(gro_path: Path) -> int:
    lines = gro_path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid gro file: {gro_path}")
    return int(lines[1].strip())


def read_gro_coords(gro_path: Path) -> tuple[str, list[str], list[tuple[float, float, float]], str]:
    lines = gro_path.read_text().splitlines()
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms] if len(lines) > 2 + n_atoms else "0 0 0"
    coords: list[tuple[float, float, float]] = []
    for line in atom_lines:
        coords.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    return title, atom_lines, coords, box_line


def write_merged_gro(
    *,
    protein_gro: Path,
    ligand_gro: Path | None,
    output_gro: Path,
    title: str = "Merged solute",
) -> None:
    p_title, p_lines, _, p_box = read_gro_coords(protein_gro)
    out_lines = list(p_lines)
    if ligand_gro is not None:
        _, l_lines, _, _ = read_gro_coords(ligand_gro)
        out_lines.extend(l_lines)
    n_atoms = len(out_lines)
    output_gro.write_text(
        title + "\n"
        + f"{n_atoms}\n"
        + "\n".join(out_lines)
        + "\n"
        + p_box
        + "\n"
    )


def extract_ligand_coords_nm(
    source_pdb: Path,
    *,
    resname: str,
    atom_names: list[str],
) -> list[tuple[float, float, float]]:
    """Extract ligand coordinates in GROMACS nm order (matches dor.top atom list)."""
    raw_atoms: list[tuple[str, float, float, float]] = []
    for line in source_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != resname:
            continue
        name = line[12:16].strip()
        x = float(line[30:38]) / 10.0
        y = float(line[38:46]) / 10.0
        z = float(line[46:54]) / 10.0
        raw_atoms.append((name, x, y, z))

    if len(raw_atoms) != len(atom_names):
        raise ValueError(
            f"Ligand atom count mismatch in {source_pdb}: "
            f"found {len(raw_atoms)}, expected {len(atom_names)} for {resname}"
        )

    # OpenMM exports ligand atoms in the same order as OpenFF / dor.top.
    return [(x, y, z) for _, x, y, z in raw_atoms]


def parse_dor_atom_names(dor_top: Path, resname: str) -> list[str]:
    names: list[str] = []
    in_atoms = False
    for line in dor_top.read_text().splitlines():
        if line.strip().startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms and line.strip().startswith("["):
            break
        if not in_atoms or not line.strip() or line.strip().startswith(";"):
            continue
        parts = line.split()
        if len(parts) >= 5 and parts[3] == resname:
            names.append(parts[4])
    if not names:
        raise ValueError(f"No atoms parsed for {resname} in {dor_top}")
    return names


def write_ligand_gro(
    *,
    coords_nm: list[tuple[float, float, float]],
    atom_names: list[str],
    resname: str,
    output_gro: Path,
    residue_number: int = 1,
) -> None:
    if len(coords_nm) != len(atom_names):
        raise ValueError("coords and atom_names length mismatch")
    lines: list[str] = ["Ligand\n", f"{len(atom_names)}\n"]
    for idx, (name, (x, y, z)) in enumerate(zip(atom_names, coords_nm), start=1):
        lines.append(
            f"{residue_number:5d}{resname:5s}{name:>5s}{idx:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
        )
    lines.append("   0.000   0.000   0.000\n")
    output_gro.write_text("".join(lines))


def write_dor_molecule_itp(dor_top: Path, output_itp: Path) -> None:
    """Strip [ system ] / [ molecules ] from exported OpenFF top → molecule .itp."""
    lines = dor_top.read_text().splitlines(keepends=True)
    out: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[ system ]") or stripped.startswith("[ molecules ]"):
            break
        out.append(line)
    output_itp.write_text("".join(out))


def append_ligand_to_topology(
    protein_top: Path,
    *,
    ligand_itp: Path,
    ligand_name: str,
    output_top: Path,
) -> None:
    text = protein_top.read_text().splitlines()
    out: list[str] = []
    molecules_idx: int | None = None

    for i, line in enumerate(text):
        if line.strip().startswith("[ molecules ]"):
            molecules_idx = i
        out.append(line)

    if molecules_idx is None:
        raise ValueError(f"No [ molecules ] section in {protein_top}")

    rel_itp = ligand_itp.name
    out.insert(molecules_idx, f'#include "{rel_itp}"\n')
    molecules_idx += 1

    # Append ligand after last molecule line
    insert_at = len(out)
    for j in range(molecules_idx + 1, len(out)):
        if out[j].strip().startswith("["):
            insert_at = j
            break
    out.insert(insert_at, f"{ligand_name:<15} 1\n")
    output_top.write_text("\n".join(out) + "\n")


def count_net_charge_from_top(top_path: Path) -> float:
    charge = 0.0
    in_atoms = False
    for line in top_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms and stripped.startswith("["):
            break
        if not in_atoms or not stripped or stripped.startswith(";"):
            continue
        parts = stripped.split()
        if len(parts) >= 7:
            try:
                charge += float(parts[6])
            except ValueError:
                continue
    return charge
