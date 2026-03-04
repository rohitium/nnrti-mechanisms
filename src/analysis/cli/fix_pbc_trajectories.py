#!/usr/bin/env python3
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np


def _topology_for_dcd(dcd_path: Path) -> Path:
    name = dcd_path.name
    if not name.endswith("_analysis.dcd"):
        raise ValueError(f"Unexpected DCD name (expected *_analysis.dcd): {dcd_path}")
    topo_name = name.replace("_analysis.dcd", "_analysis_topology.pdb")
    return dcd_path.with_name(topo_name)


def _load_mdtraj_trajectory(dcd_path: Path, topo_path: Path):
    import mdtraj as md
    import mdtraj.formats

    if str(dcd_path).endswith(".bak"):
        ref = md.load(str(topo_path))
        with mdtraj.formats.DCDTrajectoryFile(str(dcd_path), "r") as handle:
            xyz, lengths, angles = handle.read()
        return md.Trajectory(
            xyz / 10.0,  # A -> nm
            ref.topology,
            unitcell_lengths=(lengths / 10.0) if lengths is not None else None,
            unitcell_angles=angles,
        )
    return md.load(str(dcd_path), top=str(topo_path))


def _correct_one_mdtraj(
    dcd_path: Path,
    topo_path: Path,
    out_path: Path,
) -> tuple[int, int]:
    traj = _load_mdtraj_trajectory(dcd_path=dcd_path, topo_path=topo_path)
    if traj.n_atoms < 1:
        raise ValueError(f"No atoms in topology for {dcd_path}")

    # 1) Rebuild each molecule from its bond graph.
    traj.make_molecules_whole(inplace=True)

    # 2) Apply molecule-wise no-jump unwrapping across frames. This removes
    # frame-to-frame image hops that still remain after "whole molecule"
    # reconstruction (observed as abrupt RMSD/COM spikes).
    mol_indices = [
        np.asarray([atom.index for atom in mol], dtype=int)
        for mol in traj.topology.find_molecules()
    ]
    for frame_i in range(1, traj.n_frames):
        if traj.unitcell_lengths is None:
            break
        box = traj.unitcell_lengths[frame_i]
        if box is None or not np.all(np.isfinite(box)) or not np.all(box > 0):
            continue
        for mol_idx in mol_indices:
            if mol_idx.size == 0:
                continue
            prev_com = traj.xyz[frame_i - 1, mol_idx].mean(axis=0)
            curr_com = traj.xyz[frame_i, mol_idx].mean(axis=0)
            shift = -box * np.round((curr_com - prev_com) / box)
            traj.xyz[frame_i, mol_idx] += shift

    traj.save_dcd(str(out_path))
    return int(traj.n_frames), int(traj.n_atoms)


def _correct_one_mdanalysis(
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

    def _make_whole(ag):
        def _apply(ts):
            try:
                ag.unwrap(compound="segments", inplace=True)
            except Exception:
                try:
                    ag.unwrap(compound="residues", inplace=True)
                except Exception:
                    pass
            return ts
        return _apply

    u.trajectory.add_transformations(
        _make_whole(anchor),
        trans.NoJump(check_continuity=False),
        trans.center_in_box(anchor, center="geometry", wrap=False),
    )

    n_frames = 0
    with mda.Writer(str(out_path), n_atoms=u.atoms.n_atoms) as writer:
        for _ in u.trajectory:
            writer.write(u.atoms)
            n_frames += 1
    return n_frames, u.atoms.n_atoms


def _correct_one(
    dcd_path: Path,
    topo_path: Path,
    out_path: Path,
    *,
    backend: str,
) -> tuple[int, int]:
    if backend == "mdtraj":
        return _correct_one_mdtraj(dcd_path=dcd_path, topo_path=topo_path, out_path=out_path)
    if backend == "mdanalysis":
        return _correct_one_mdanalysis(dcd_path=dcd_path, topo_path=topo_path, out_path=out_path)
    raise ValueError(f"Unsupported backend: {backend}")


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
    parser.add_argument(
        "--backend",
        choices=["mdtraj", "mdanalysis"],
        default="mdtraj",
        help="PBC correction backend (default: mdtraj).",
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
                n_frames, n_atoms = _correct_one(dcd_path, topo_path, tmp_path, backend=args.backend)
                if backup_path.exists():
                    backup_path.unlink()
                dcd_path.replace(backup_path)
                tmp_path.replace(dcd_path)
                print(
                    f"[ok] {dcd_path} frames={n_frames} atoms={n_atoms} "
                    f"backup={backup_path.name} backend={args.backend}"
                )
            else:
                out_path = dcd_path.with_name(f"{dcd_path.stem}{args.suffix}.dcd")
                n_frames, n_atoms = _correct_one(dcd_path, topo_path, out_path, backend=args.backend)
                print(f"[ok] {out_path} frames={n_frames} atoms={n_atoms} backend={args.backend}")
            ok += 1
        except Exception as exc:
            print(f"[fail] {dcd_path}: {exc}")
            failed += 1

    print(f"Completed: ok={ok} failed={failed} total={len(dcd_paths)}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
