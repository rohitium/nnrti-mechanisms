#!/usr/bin/env python3
"""Write PBC-imaged analysis DCDs once, then reuse them for metrics.

Default: results/md_runs/{mut}/rep_XX/{safe}_repXX_analysis_pbcfix.dcd
Skips _archive. Copies NSTEP/DELTA/NSAVC from the source so timing survives mdtraj.save_dcd.
"""
from __future__ import annotations

import argparse
import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from ..pbc import apply_mdtraj_pbc_correction, load_mdtraj_trajectory, pbcfix_dcd_for, topology_for_analysis_dcd

SKIP_PARTS = {"_archive", "apo", "sherlock_rerun"}


def _copy_timing_header(src: Path, dest: Path) -> None:
    raw = src.read_bytes()[:92]
    if raw[4:8] != b"CORD":
        return
    with dest.open("r+b") as handle:
        handle.seek(16)
        handle.write(raw[16:20])  # NSAVC
        handle.seek(20)
        handle.write(raw[20:24])  # NSTEP
        handle.seek(44)
        handle.write(raw[44:48])  # DELTA


def _correct_one(dcd_path: str, topo_path: str, out_path: str) -> tuple[int, int]:
    import tempfile

    src = Path(dcd_path)
    dest = Path(out_path)
    traj = load_mdtraj_trajectory(dcd_path=src, topo_path=Path(topo_path))
    if traj.n_atoms < 1:
        raise ValueError(f"No atoms in topology for {src}")
    apply_mdtraj_pbc_correction(traj, anchor_selection="protein", ligand_resname="2KW")
    write_to = dest
    tmp_path = None
    if dest.resolve() == src.resolve():
        with tempfile.NamedTemporaryFile(
            prefix=src.stem + ".",
            suffix=".tmp.dcd",
            dir=str(src.parent),
            delete=False,
        ) as tmp:
            tmp_path = Path(tmp.name)
        write_to = tmp_path
    traj.save_dcd(str(write_to))
    _copy_timing_header(src, write_to)
    if tmp_path is not None:
        backup = src.with_name(src.name + ".bak")
        if backup.exists():
            backup.unlink()
        src.replace(backup)
        tmp_path.replace(src)
    return int(traj.n_frames), int(traj.n_atoms)


def _keep(path: Path) -> bool:
    return not any(part in SKIP_PARTS for part in path.parts)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write PBC-imaged *_analysis_pbcfix.dcd files.")
    parser.add_argument("--root", type=Path, default=Path("results/md_runs"))
    parser.add_argument("--pattern", type=str, default="*_analysis.dcd")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true", help="Overwrite existing pbcfix DCDs.")
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the source DCD after writing a .bak (legacy). Prefer sidecar pbcfix files.",
    )
    args = parser.parse_args()

    dcd_paths = [p for p in sorted(args.root.rglob(args.pattern)) if _keep(p)]
    jobs = []
    skipped = 0
    for dcd_path in dcd_paths:
        topo_path = topology_for_analysis_dcd(dcd_path)
        if not topo_path.exists():
            print(f"[skip] missing topology: {topo_path}")
            continue
        if args.in_place:
            out_path = dcd_path
            if dcd_path.with_name(dcd_path.name + ".bak").exists() and not args.force:
                skipped += 1
                continue
        else:
            out_path = pbcfix_dcd_for(dcd_path)
            if out_path.exists() and not args.force:
                skipped += 1
                continue
        jobs.append((str(dcd_path), str(topo_path), str(out_path)))

    print(f"imaging {len(jobs)} DCDs ({skipped} already present) with {args.workers} workers")
    ok = 0
    failed = 0
    if jobs:
        with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
            futs = {pool.submit(_correct_one, *job): job for job in jobs}
            done = 0
            for fut in as_completed(futs):
                dcd_path, _topo, out_path = futs[fut]
                done += 1
                try:
                    n_frames, n_atoms = fut.result()
                    ok += 1
                    print(f"[ok {done}/{len(jobs)}] {out_path} frames={n_frames} atoms={n_atoms}")
                except Exception as exc:
                    failed += 1
                    print(f"[fail {done}/{len(jobs)}] {dcd_path}: {exc}")

    print(f"Completed: ok={ok} failed={failed} skipped={skipped} scanned={len(dcd_paths)}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
