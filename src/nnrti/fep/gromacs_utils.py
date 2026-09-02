"""GROMACS helpers for pmx NEQ system builds."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


class GromacsError(RuntimeError):
    pass


def find_gmx() -> str:
    gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
    if gmx is None:
        raise GromacsError(
            "gmx not found on PATH. On Sherlock: source ops/slurm/cluster/load_gromacs_module.sh"
        )
    return gmx


def resolve_gmxlib(env: dict[str, str] | None = None) -> str:
    """Return pmx mutff path for hybrid topology includes."""
    env = env or os.environ
    existing = env.get("GMXLIB", "").strip()
    if existing:
        return existing

    try:
        import pmx  # type: ignore import-not-found
    except ImportError as exc:
        raise GromacsError(
            "GMXLIB not set and pmx is not importable. "
            "On Sherlock: module load python/3.9.0 && source ~/.venvs/pmx/bin/activate "
            "before submitting, or export GMXLIB to pmx/data/mutff."
        ) from exc

    gmxlib = os.path.join(os.path.dirname(pmx.__file__), "data", "mutff")
    if not os.path.isdir(gmxlib):
        raise GromacsError(f"pmx mutff directory not found: {gmxlib}")
    return gmxlib


def gromacs_env(base: dict[str, str] | None = None) -> dict[str, str]:
    """Copy env and ensure GMXLIB is set for hybrid topologies."""
    env = dict(base or os.environ)
    env["GMXLIB"] = resolve_gmxlib(env)
    return env


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


_GENERIC_HIS_RESNAMES = frozenset({"HIS", "HSD", "HSE", "HSP", "HISD", "HISE", "HISH"})


def _amber_his_resname(atom_names: set[str]) -> str:
    """Infer Amber histidine protonation from input hydrogen names."""
    has_hd1 = "HD1" in atom_names
    has_he2 = "HE2" in atom_names
    if has_hd1 and has_he2:
        return "HIP"
    if has_hd1:
        return "HID"
    if has_he2:
        return "HIE"
    return "HIE"


def normalize_hybrid_his_for_pdb2gmx(pdb_path: Path) -> int:
    """Rename generic HIS residues so pdb2gmx uses input protonation.

    OpenMM MD structures export HIS with HD1/HE2 atom names but generic HIS
    resnames. pdb2gmx then reassigns protonation via its H-bond network and
    fails when input atoms (e.g. HD1) do not match the chosen state (HIE).
    """
    lines = pdb_path.read_text().splitlines()
    residue_atoms: dict[tuple[str, str, str], set[str]] = {}
    atom_line_indices: dict[tuple[str, str, str], list[int]] = {}

    for idx, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname not in _GENERIC_HIS_RESNAMES and resname not in {"HID", "HIE", "HIP"}:
            continue
        key = (line[21], line[22:26].strip(), resname)
        atom_name = line[12:16].strip()
        residue_atoms.setdefault(key, set()).add(atom_name)
        atom_line_indices.setdefault(key, []).append(idx)

    renames: dict[tuple[str, str, str], str] = {}
    for key, atoms in residue_atoms.items():
        _, _, resname = key
        if resname in {"HID", "HIE", "HIP"}:
            continue
        renames[key] = _amber_his_resname(atoms)

    if not renames:
        return 0

    for key, new_resname in renames.items():
        for idx in atom_line_indices[key]:
            line = lines[idx]
            lines[idx] = f"{line[:17]}{new_resname:>3s}{line[20:]}"

    pdb_path.write_text("\n".join(lines) + ("\n" if lines else ""))
    return len(renames)


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


def write_dor_ligand_itps(
    dor_top: Path,
    *,
    atomtypes_itp: Path,
    molecule_itp: Path,
) -> None:
    """Split OpenFF dor.top into atomtypes + molecule itps for protein merge."""
    lines = dor_top.read_text().splitlines(keepends=True)
    atomtypes_lines: list[str] = []
    molecule_lines: list[str] = []
    section: str | None = None

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            if section == "defaults":
                continue
            if section in ("system", "molecules"):
                break
            if section == "atomtypes":
                atomtypes_lines.append(line)
                continue
            if section == "moleculetype" or molecule_lines:
                molecule_lines.append(line)
                continue
            continue
        if section == "atomtypes":
            atomtypes_lines.append(line)
        elif section not in (None, "defaults") and (section == "moleculetype" or molecule_lines):
            molecule_lines.append(line)

    if not atomtypes_lines or not molecule_lines:
        raise ValueError(f"Failed to split ligand topology from {dor_top}")

    atomtypes_itp.write_text("".join(atomtypes_lines))
    molecule_itp.write_text("".join(molecule_lines))


def _ligand_atomtypes_insert_index(lines: list[str]) -> int:
    """Return index to insert ligand atomtypes after forcefield includes."""
    last_ff_include = 0
    first_molecule_include = len(lines)
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("#include") and ".ff/" in stripped:
            last_ff_include = i + 1
        elif stripped.startswith("#include"):
            first_molecule_include = i
            break
    return first_molecule_include if first_molecule_include < len(lines) else last_ff_include


def append_ligand_to_topology(
    protein_top: Path,
    *,
    ligand_atomtypes_itp: Path,
    ligand_itp: Path,
    ligand_name: str,
    output_top: Path,
) -> None:
    out = protein_top.read_text().splitlines()

    def _molecules_header_index() -> int:
        for i, line in enumerate(out):
            if line.strip().startswith("[ molecules ]"):
                return i
        raise ValueError(f"No [ molecules ] section in {protein_top}")

    # Sanity-check the section exists before mutating the buffer.
    _molecules_header_index()

    # 1) Ligand [ atomtypes ] include: after the forcefield includes, before the
    #    first molecule-type include.
    atomtypes_idx = _ligand_atomtypes_insert_index(out)
    out.insert(atomtypes_idx, f'#include "{ligand_atomtypes_itp.name}"')

    # 2) Ligand moleculetype include: immediately before [ molecules ]. Re-find the
    #    header rather than reusing a pre-insert index (which the atomtypes insert
    #    above shifts) so this stays correct regardless of how many lines moved.
    molecules_idx = _molecules_header_index()
    out.insert(molecules_idx, f'#include "{ligand_itp.name}"')

    # 3) Append the ligand to the [ molecules ] list, after the last molecule row
    #    (i.e. before the next section header, or at EOF if it is the last section).
    molecules_idx = _molecules_header_index()
    insert_at = len(out)
    for j in range(molecules_idx + 1, len(out)):
        if out[j].strip().startswith("["):
            insert_at = j
            break
    out.insert(insert_at, f"{ligand_name:<15} 1")
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
