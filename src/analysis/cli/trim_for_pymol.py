#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write protein+ligand-only PDB/DCD files for easier PyMOL loading."
    )
    parser.add_argument("--topology", type=Path, required=True, help="Input topology PDB")
    parser.add_argument("--trajectory", type=Path, required=True, help="Input trajectory DCD")
    parser.add_argument("--ligand-resname", type=str, required=True, help="Ligand residue name (e.g. 2KW)")
    parser.add_argument("--out-pdb", type=Path, required=True, help="Output trimmed PDB")
    parser.add_argument("--out-dcd", type=Path, required=True, help="Output trimmed DCD")
    parser.add_argument(
        "--no-pbc-correct",
        action="store_true",
        help="Disable NoJump + center_in_box transform before writing.",
    )
    args = parser.parse_args()

    import MDAnalysis as mda

    u = mda.Universe(str(args.topology), str(args.trajectory))
    if not args.no_pbc_correct:
        try:
            from MDAnalysis import transformations as trans

            prot = u.select_atoms("protein")
            anchor = prot if prot.n_atoms > 0 else u.atoms
            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(anchor, center="geometry", wrap=False),
            )
        except Exception as exc:
            print(f"Warning: failed to apply PBC correction transforms ({exc})")

    sel = u.select_atoms(f"protein or resname {args.ligand_resname}")
    if sel.n_atoms == 0:
        raise SystemExit(f"Selection returned 0 atoms: protein or resname {args.ligand_resname}")

    args.out_pdb.parent.mkdir(parents=True, exist_ok=True)
    args.out_dcd.parent.mkdir(parents=True, exist_ok=True)

    u.trajectory[0]
    sel.write(str(args.out_pdb))

    with mda.Writer(str(args.out_dcd), n_atoms=sel.n_atoms) as writer:
        for _ in u.trajectory:
            writer.write(sel)

    print(f"Selection atoms: {sel.n_atoms}")
    print(f"Frames written: {len(u.trajectory)}")
    print(f"Wrote PDB: {args.out_pdb}")
    print(f"Wrote DCD: {args.out_dcd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
