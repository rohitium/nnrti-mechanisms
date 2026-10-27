#!/usr/bin/env python3
"""Generate perturbed DOR starting poses as modified mmCIF files.

The holo prep pipeline reads the bound pose from the CIF: ``minimize_with_restraints``
deletes the CIF ligand and rigidly aligns the SDF template onto its heavy atoms, so a
CIF with rewritten ligand coordinates is all that is needed to seed MD from a different
starting pose.  Only the ``Cartn_x/y/z`` columns of the ligand rows are touched; every
other line is passed through byte-identical.

Two pose families are produced:

``rot_pa{1,2,3}_{angle}``
    Rigid rotation about each principal axis of the ligand heavy atoms, through the
    heavy-atom centroid.  Axis 1 is the long axis, so that rotation is close to
    isosteric; axes 2 and 3 swap the ends of the molecule and clash hard.

``alt_{nn}``
    Random rigid-body placements screened for pocket occupancy and sterics: heavy-atom
    RMSD above ``--min-rmsd`` from the crystal pose, centroid still within
    ``--max-centroid-shift`` of the crystal centroid, and no protein heavy atom closer
    than ``--min-contact``.  This is pocket-constrained random placement rather than
    docking, which keeps the test free of an external docking dependency.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# _atom_site column order in the loop header, used to index whitespace-split rows.
_COL_GROUP = 0
_COL_TYPE_SYMBOL = 2
_COL_LABEL_COMP = 5
_COL_X = 10
_COL_Y = 11
_COL_Z = 12
_COL_AUTH_COMP = 17

_SOLVENT_COMPS = {"HOH", "WAT", "DOD"}

# Heavy-atom indices in 4NCG CIF order, defining the ligand long axis: the
# chlorocyanophenyl chlorine against the triazolinone N/O centroid.  A fixed molecular
# vector rather than an inertial eigenvector, so an end-for-end flip reads as 180
# degrees instead of folding back to 0.
_CL_INDEX = 20
_TRIAZOLINONE_INDICES = [25, 26, 27, 28]


@dataclass(frozen=True)
class AtomSiteBlock:
    """The parsed ``_atom_site`` loop of an mmCIF file."""

    lines: list[str]
    first_row: int
    last_row: int
    rows: list[list[str]]
    row_line_numbers: list[int]

    def coords(self, indices) -> np.ndarray:
        return np.array(
            [
                [
                    float(self.rows[i][_COL_X]),
                    float(self.rows[i][_COL_Y]),
                    float(self.rows[i][_COL_Z]),
                ]
                for i in indices
            ],
            dtype=float,
        )


def parse_atom_site(cif_path: Path) -> AtomSiteBlock:
    lines = cif_path.read_text().splitlines(keepends=True)
    header_seen = False
    first_row = None
    rows: list[list[str]] = []
    row_line_numbers: list[int] = []

    for idx, line in enumerate(lines):
        if line.startswith("_atom_site."):
            header_seen = True
            continue
        if not header_seen:
            continue
        if line.startswith(("ATOM ", "HETATM")):
            if first_row is None:
                first_row = idx
            rows.append(line.split())
            row_line_numbers.append(idx)
        elif first_row is not None:
            break

    if first_row is None:
        raise ValueError(f"No _atom_site rows found in {cif_path}")

    return AtomSiteBlock(
        lines=lines,
        first_row=first_row,
        last_row=row_line_numbers[-1],
        rows=rows,
        row_line_numbers=row_line_numbers,
    )


def ligand_rows(block: AtomSiteBlock, resname: str) -> list[int]:
    hits = [
        i
        for i, row in enumerate(block.rows)
        if row[_COL_LABEL_COMP] == resname or row[_COL_AUTH_COMP] == resname
    ]
    if not hits:
        raise ValueError(f"Ligand residue '{resname}' not present in _atom_site rows")
    return hits


def heavy_mask(block: AtomSiteBlock, indices) -> np.ndarray:
    return np.array(
        [block.rows[i][_COL_TYPE_SYMBOL].upper() != "H" for i in indices],
        dtype=bool,
    )


def environment_rows(block: AtomSiteBlock, ligand_idx: set[int]) -> list[int]:
    """Non-ligand, non-solvent heavy atoms, used for the steric screen."""
    return [
        i
        for i, row in enumerate(block.rows)
        if i not in ligand_idx
        and row[_COL_TYPE_SYMBOL].upper() != "H"
        and row[_COL_AUTH_COMP] not in _SOLVENT_COMPS
        and row[_COL_LABEL_COMP] not in _SOLVENT_COMPS
    ]


def principal_axes(xyz: np.ndarray) -> np.ndarray:
    """Principal axes of the point set, ordered longest extent first."""
    centred = xyz - xyz.mean(axis=0)
    cov = centred.T @ centred
    values, vectors = np.linalg.eigh(cov)
    return vectors[:, np.argsort(values)[::-1]].T


def rotation_about_axis(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    theta = np.deg2rad(angle_deg)
    k = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    return np.eye(3) + np.sin(theta) * k + (1.0 - np.cos(theta)) * (k @ k)


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniform random rotation matrix via a random unit quaternion."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def apply_transform(xyz: np.ndarray, rot: np.ndarray, centroid: np.ndarray, shift: np.ndarray) -> np.ndarray:
    return ((xyz - centroid) @ rot.T) + centroid + shift


def long_axis(xyz: np.ndarray) -> np.ndarray:
    vector = xyz[_CL_INDEX] - xyz[_TRIAZOLINONE_INDICES].mean(axis=0)
    return vector / np.linalg.norm(vector)


def pose_diagnostics(
    moved_all: np.ndarray,
    heavy: np.ndarray,
    crystal_heavy: np.ndarray,
    env_xyz: np.ndarray,
) -> dict:
    moved_heavy = moved_all[heavy]
    rmsd = float(np.sqrt(np.mean(np.sum((moved_heavy - crystal_heavy) ** 2, axis=1))))
    distances = np.linalg.norm(moved_heavy[:, None, :] - env_xyz[None, :, :], axis=2)
    nearest = distances.min(axis=1)
    cosine = float(np.clip(np.dot(long_axis(moved_heavy), long_axis(crystal_heavy)), -1.0, 1.0))
    return {
        "rmsd_angstrom": round(rmsd, 3),
        "long_axis_tilt_deg": round(float(np.degrees(np.arccos(cosine))), 2),
        "centroid_shift_angstrom": round(
            float(np.linalg.norm(moved_heavy.mean(axis=0) - crystal_heavy.mean(axis=0))), 3
        ),
        "min_contact_angstrom": round(float(nearest.min()), 3),
        "n_contacts_below_2.0": int((distances < 2.0).sum()),
        "n_contacts_below_2.5": int((distances < 2.5).sum()),
        "n_contacts_below_3.0": int((distances < 3.0).sum()),
    }


def write_perturbed_cif(
    block: AtomSiteBlock,
    ligand_idx: list[int],
    new_xyz: np.ndarray,
    output_path: Path,
) -> None:
    lines = list(block.lines)
    for slot, row_index in enumerate(ligand_idx):
        row = list(block.rows[row_index])
        row[_COL_X] = f"{new_xyz[slot, 0]:.3f}"
        row[_COL_Y] = f"{new_xyz[slot, 1]:.3f}"
        row[_COL_Z] = f"{new_xyz[slot, 2]:.3f}"
        lines[block.row_line_numbers[row_index]] = " ".join(row) + "\n"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(lines))


def generate_poses(
    cif_path: Path,
    output_dir: Path,
    *,
    ligand_resname: str,
    rotation_angles: list[float],
    n_alternative: int,
    min_rmsd: float,
    max_rmsd: float | None,
    max_centroid_shift: float,
    min_contact: float,
    seed: int,
    max_attempts: int,
    prefix: str = "alt",
    min_tilt: float = 0.0,
    max_displacement: float | None = None,
) -> list[dict]:
    block = parse_atom_site(cif_path)
    ligand_idx = ligand_rows(block, ligand_resname)
    heavy = heavy_mask(block, ligand_idx)
    ligand_xyz = block.coords(ligand_idx)
    crystal_heavy = ligand_xyz[heavy]
    centroid = crystal_heavy.mean(axis=0)
    axes = principal_axes(crystal_heavy)

    # Prune the environment to a shell around the pocket. Nothing beyond the ligand's own
    # reach plus the allowed centroid shift plus the widest contact cutoff can affect any
    # reported distance, and the screen runs on a few hundred atoms instead of thousands.
    env_all = block.coords(environment_rows(block, set(ligand_idx)))
    reach = (
        float(np.linalg.norm(crystal_heavy - centroid, axis=1).max())
        + max_centroid_shift
        + 3.0
    )
    env_xyz = env_all[np.linalg.norm(env_all - centroid, axis=1) <= reach]

    logging.info(
        "%s: %d ligand atoms (%d heavy), %d environment heavy atoms",
        ligand_resname,
        len(ligand_idx),
        int(heavy.sum()),
        len(env_xyz),
    )

    records: list[dict] = []

    crystal = pose_diagnostics(ligand_xyz, heavy, crystal_heavy, env_xyz)
    crystal.update({"pose": "crystal", "family": "reference", "accepted": True})
    records.append(crystal)

    for axis_number, axis in enumerate(axes, start=1):
        for angle in rotation_angles:
            name = f"rot_pa{axis_number}_{int(angle)}"
            rot = rotation_about_axis(axis, angle)
            moved = apply_transform(ligand_xyz, rot, centroid, np.zeros(3))
            record = pose_diagnostics(moved, heavy, crystal_heavy, env_xyz)
            record.update({"pose": name, "family": "rotation", "accepted": True})
            records.append(record)
            write_perturbed_cif(block, ligand_idx, moved, output_dir / f"{name}.cif")

    rng = np.random.default_rng(seed)
    accepted: list[np.ndarray] = []
    attempts = 0
    while len(accepted) < n_alternative and attempts < max_attempts:
        attempts += 1
        rot = random_rotation(rng)
        shift = rng.normal(scale=max_centroid_shift / 2.0, size=3)
        if np.linalg.norm(shift) > max_centroid_shift:
            continue
        moved = apply_transform(ligand_xyz, rot, centroid, shift)
        diag = pose_diagnostics(moved, heavy, crystal_heavy, env_xyz)
        if diag["rmsd_angstrom"] < min_rmsd:
            continue
        if max_rmsd is not None and diag["rmsd_angstrom"] > max_rmsd:
            continue
        if diag["min_contact_angstrom"] < min_contact:
            continue
        if diag["long_axis_tilt_deg"] < min_tilt:
            continue
        if max_displacement is not None and diag["centroid_shift_angstrom"] > max_displacement:
            continue
        moved_heavy = moved[heavy]
        # keep accepted poses distinct from one another, but do not demand a
        # separation as wide as the band itself -- in a narrow band that is
        # unsatisfiable
        separation = min_rmsd if max_rmsd is None else max(0.5, (max_rmsd - min_rmsd) / 2.0)
        if any(
            np.sqrt(np.mean(np.sum((moved_heavy - prior[heavy]) ** 2, axis=1))) < separation
            for prior in accepted
        ):
            continue
        accepted.append(moved)
        name = f"{prefix}_{len(accepted):02d}"
        diag.update({"pose": name, "family": "alternative", "accepted": True})
        records.append(diag)
        write_perturbed_cif(block, ligand_idx, moved, output_dir / f"{name}.cif")

    logging.info(
        "Accepted %d/%d alternative poses after %d attempts",
        len(accepted),
        n_alternative,
        attempts,
    )
    if len(accepted) < n_alternative:
        logging.warning(
            "Steric screen is tight: loosen --min-contact or --min-rmsd, "
            "or raise --max-attempts"
        )
    return records


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    root = Path(__file__).resolve().parents[3]
    parser.add_argument("--cif", type=Path, default=root / "data" / "structures" / "4NCG.cif")
    parser.add_argument(
        "--output-dir", type=Path, default=root / "data" / "prepared" / "dor_poses"
    )
    parser.add_argument("--ligand-resname", default="2KW")
    parser.add_argument(
        "--rotation-angle",
        type=float,
        nargs="*",
        default=[45.0, 60.0, 90.0],
        dest="rotation_angles",
    )
    parser.add_argument("--n-alternative", type=int, default=3)
    parser.add_argument("--min-rmsd", type=float, default=2.0)
    parser.add_argument(
        "--max-rmsd", type=float, default=None,
        help="Only keep placements whose heavy-atom RMSD is below this (A). "
             "With --min-rmsd this selects a band, which is what a "
             "radius-of-convergence ladder needs.",
    )
    parser.add_argument("--max-centroid-shift", type=float, default=2.5)
    parser.add_argument("--min-contact", type=float, default=1.9)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--max-attempts", type=int, default=200000)
    parser.add_argument("--prefix", default="alt", help="Stem for the generated placements.")
    parser.add_argument(
        "--min-tilt", type=float, default=0.0,
        help="Only keep placements whose long-axis tilt exceeds this (degrees).",
    )
    parser.add_argument(
        "--max-displacement", type=float, default=None,
        help="Only keep placements whose centroid moves less than this (A).",
    )
    args = parser.parse_args(argv)

    records = generate_poses(
        args.cif,
        args.output_dir,
        ligand_resname=args.ligand_resname,
        rotation_angles=args.rotation_angles,
        n_alternative=args.n_alternative,
        min_rmsd=args.min_rmsd,
        max_rmsd=args.max_rmsd,
        max_centroid_shift=args.max_centroid_shift,
        min_contact=args.min_contact,
        seed=args.seed,
        max_attempts=args.max_attempts,
        prefix=args.prefix,
        min_tilt=args.min_tilt,
        max_displacement=args.max_displacement,
    )

    summary = args.output_dir / "pose_diagnostics.csv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "pose",
        "family",
        "rmsd_angstrom",
        "long_axis_tilt_deg",
        "centroid_shift_angstrom",
        "min_contact_angstrom",
        "n_contacts_below_2.0",
        "n_contacts_below_2.5",
        "n_contacts_below_3.0",
    ]
    with open(summary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    logging.info("Wrote %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
