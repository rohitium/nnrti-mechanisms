#!/usr/bin/env python3
"""Export DOR (OpenFF 2.0.0) to GROMACS .top/.gro for holo system builds."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import DOR_ITP_DIR, LIGAND_SDF, OPENFF_FORCE_FIELD


def export_dor_gromacs(
    ligand_sdf: Path,
    output_dir: Path,
    *,
    prefix: str = "dor",
    box_nm: float = 5.0,
) -> tuple[Path, Path]:
    from openff.interchange import Interchange
    from openff.toolkit import ForceField, Molecule
    from openff.units import unit

    if not ligand_sdf.is_file():
        raise FileNotFoundError(ligand_sdf)

    output_dir.mkdir(parents=True, exist_ok=True)
    mol = Molecule.from_file(ligand_sdf)
    mol.name = "2KW"
    if mol.n_conformers == 0:
        mol.generate_conformers(n_conformers=1)

    force_field = ForceField(OPENFF_FORCE_FIELD)
    box = box_nm * np.eye(3) * unit.nanometer
    interchange = Interchange.from_smirnoff(force_field=force_field, topology=[mol], box=box)
    out_prefix = str(output_dir / prefix)
    interchange.to_gromacs(prefix=out_prefix)

    top_path = output_dir / f"{prefix}.top"
    gro_path = output_dir / f"{prefix}.gro"
    if not top_path.is_file() or not gro_path.is_file():
        raise RuntimeError(f"Expected GROMACS files under {output_dir}")
    return gro_path, top_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export DOR OpenFF parameters to GROMACS.")
    parser.add_argument("--ligand-sdf", type=Path, default=LIGAND_SDF)
    parser.add_argument("--output-dir", type=Path, default=DOR_ITP_DIR)
    parser.add_argument("--prefix", default="dor")
    parser.add_argument(
        "--conda-env",
        default="nnrti-prep",
        help="Hint only: run with `conda activate nnrti-prep` (needs openff-interchange).",
    )
    args = parser.parse_args(argv)

    try:
        gro, top = export_dor_gromacs(args.ligand_sdf, args.output_dir, prefix=args.prefix)
    except ModuleNotFoundError as exc:
        print(
            f"Missing dependency ({exc}). Activate conda env: conda activate {args.conda_env}",
            file=sys.stderr,
        )
        return 1

    print(f"Wrote {gro}")
    print(f"Wrote {top}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
