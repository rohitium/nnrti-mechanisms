#!/usr/bin/env python3
"""Write DOR with AM1-BCC partial charges baked into an SDF.

DOR is the same molecule in every system, so its AM1-BCC charges need to be
derived once. Writing them into an SDF lets the scoring job load the ligand
with charges already assigned, so the cluster does not need AmberTools --
only a toolkit that can read the file.

Usage
-----
    PYTHONPATH=src python ops/slurm/apo_collapse/export_am1bcc_charges.py \
        data/ligands/dor.sdf data/ligands/dor_am1bcc.sdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

import numpy as np
from nnrti.md.openmm.ligand import load_ligand_molecule


def main() -> int:
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    mol = load_ligand_molecule(src, "am1bcc")
    q = np.array([float(x.m) for x in mol.partial_charges])
    print(f"{mol.n_atoms} atoms, net charge {q.sum():+.6f}, mean |q| {np.abs(q).mean():.4f}")
    if abs(q.sum()) > 1e-4:
        raise SystemExit(f"net charge {q.sum()} is not neutral; refusing to write")
    mol.to_file(str(dst), file_format="sdf")

    # read back and confirm the charges survived the round trip
    from openff.toolkit.topology import Molecule
    back = Molecule.from_file(str(dst), allow_undefined_stereo=True)
    qb = np.array([float(x.m) for x in back.partial_charges])
    if not np.allclose(q, qb, atol=1e-6):
        raise SystemExit("charges did not survive the SDF round trip")
    print(f"wrote {dst}; charges round-trip to {np.abs(q - qb).max():.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
