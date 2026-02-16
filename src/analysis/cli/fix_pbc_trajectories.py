#!/usr/bin/env python3
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path


def _topology_for_dcd(dcd_path: Path) -> Path:
    name = dcd_path.name
    if not name.endswith("_analysis.dcd"):
        raise ValueError(f"Unexpected DCD name (expected *_analysis.dcd): {dcd_path}")
    topo_name = name.replace("_analysis.dcd", "_analysis_topology.pdb")
    return dcd_path.with_name(topo_name)


def _correct_one(
    dcd_path: Path,
    topo_path: Path,
    out_path: Path,
) -> tuple[int, int]:
    import MDAnalysis as mda
    from MDAnalysis import transformations as trans

    u = mda.Universe(str(topo_path), str(dcd_path))
    if u.atoms.n_atoms == 0:
        raise ValueError(f"No atoms in topology for {dcd_path}")

    protein = u.select_atoms("protein")
    anchor = protein if protein.n_atoms > 0 else u.atoms
    u.trajectory.add_transformations(
        trans.NoJump(check_continuity=False),
        trans.center_in_box(anchor, center="geometry", wrap=False),
    )

    n_frames = 0
    with mda.Writer(str(out_path), n_atoms=u.atoms.n_atoms) as writer:
        for _ in u.trajectory:
            writer.write(u.atoms)
            n_frames += 1
    return n_frames, u.atoms.n_atoms


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Apply PBC correction (NoJump + center_in_box) to analysis trajectories."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("results/md_runs"),
        help="Root folder to scan for *_analysis.dcd",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*_analysis.dcd",
        help="Glob pattern under --root",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite original DCDs (default writes *_pbcfix.dcd)",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="_pbcfix",
        help="Suffix added before .dcd when not using --in-place",
    )
    parser.add_argument(
        "--backup-ext",
        type=str,
        default=".bak",
        help="Backup extension used only with --in-place",
    )
    args = parser.parse_args()

    root = args.root
    dcd_paths = sorted(root.rglob(args.pattern))
    if not dcd_paths:
        print(f"No trajectories found under {root} matching {args.pattern}")
        return 0

    ok = 0
    failed = 0
    for dcd_path in dcd_paths:
        topo_path = _topology_for_dcd(dcd_path)
        if not topo_path.exists():
            print(f"[skip] missing topology: {topo_path}")
            failed += 1
            continue

        try:
            if args.in_place:
                backup_path = dcd_path.with_name(dcd_path.name + args.backup_ext)
                with tempfile.NamedTemporaryFile(
                    prefix=dcd_path.stem + ".",
                    suffix=".tmp.dcd",
                    dir=str(dcd_path.parent),
                    delete=False,
                ) as tmp:
                    tmp_path = Path(tmp.name)
                n_frames, n_atoms = _correct_one(dcd_path, topo_path, tmp_path)
                if backup_path.exists():
                    backup_path.unlink()
                dcd_path.replace(backup_path)
                tmp_path.replace(dcd_path)
                print(f"[ok] {dcd_path} frames={n_frames} atoms={n_atoms} backup={backup_path.name}")
            else:
                out_path = dcd_path.with_name(f"{dcd_path.stem}{args.suffix}.dcd")
                n_frames, n_atoms = _correct_one(dcd_path, topo_path, out_path)
                print(f"[ok] {out_path} frames={n_frames} atoms={n_atoms}")
            ok += 1
        except Exception as exc:
            print(f"[fail] {dcd_path}: {exc}")
            failed += 1

    print(f"Completed: ok={ok} failed={failed} total={len(dcd_paths)}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
