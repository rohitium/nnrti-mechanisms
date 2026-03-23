#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ..pbc import (
    apply_mdtraj_pbc_correction,
    audit_mdtraj_trajectory,
    load_mdtraj_trajectory,
    topology_for_analysis_dcd,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit MDTraj PBC-corrected analysis trajectories and flag residual artifacts.",
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
        "--output-csv",
        type=Path,
        default=Path("results/tables/pbc_audit.csv"),
        help="Where to write the audit table.",
    )
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=10,
        help="Audit every Nth frame (default: 10).",
    )
    parser.add_argument(
        "--apply-correction-in-memory",
        action="store_true",
        help="Apply the MDTraj PBC fix before auditing, without modifying files on disk.",
    )
    parser.add_argument(
        "--ligand-resname",
        type=str,
        default="2KW",
        help="Ligand residue name used for ligand/protein audit checks.",
    )
    parser.add_argument(
        "--max-bond-length-angstrom",
        type=float,
        default=3.0,
        help="Fail if any sampled bonded pair exceeds this distance.",
    )
    parser.add_argument(
        "--max-anchor-center-offset-angstrom",
        type=float,
        default=3.0,
        help="Fail if the protein COM is farther than this from the box center.",
    )
    parser.add_argument(
        "--max-ligand-anchor-gap-angstrom",
        type=float,
        default=1.0,
        help="Fail if ligand/protein direct vs min-image COM distances differ by more than this.",
    )
    parser.add_argument(
        "--max-anchor-internal-jump-angstrom",
        type=float,
        default=15.0,
        help="Fail if protein-anchor molecule relative COM jumps exceed this between sampled frames.",
    )
    parser.add_argument(
        "--max-ligand-anchor-jump-angstrom",
        type=float,
        default=8.0,
        help="Fail if ligand COM jumps relative to the protein exceed this between sampled frames.",
    )
    args = parser.parse_args()

    dcd_paths = sorted(args.root.rglob(args.pattern))
    if not dcd_paths:
        print(f"No trajectories found under {args.root} matching {args.pattern}")
        return 0

    rows: list[dict[str, object]] = []
    failures = 0
    for dcd_path in dcd_paths:
        topo_path = topology_for_analysis_dcd(dcd_path)
        if not topo_path.exists():
            print(f"[fail] missing topology: {topo_path}")
            failures += 1
            continue

        try:
            traj = load_mdtraj_trajectory(dcd_path=dcd_path, topo_path=topo_path)
            if args.apply_correction_in_memory:
                apply_mdtraj_pbc_correction(
                    traj,
                    anchor_selection="protein",
                    ligand_resname=args.ligand_resname,
                )
            summary = audit_mdtraj_trajectory(
                traj,
                dcd_path=dcd_path,
                topo_path=topo_path,
                anchor_selection="protein",
                ligand_resname=args.ligand_resname,
                frame_stride=args.frame_stride,
                max_bond_length_angstrom=args.max_bond_length_angstrom,
                max_anchor_center_offset_angstrom=args.max_anchor_center_offset_angstrom,
                max_ligand_anchor_gap_angstrom=args.max_ligand_anchor_gap_angstrom,
                max_anchor_internal_jump_angstrom=args.max_anchor_internal_jump_angstrom,
                max_ligand_anchor_jump_angstrom=args.max_ligand_anchor_jump_angstrom,
            )
            rows.append(summary.to_dict())
            status = "ok" if summary.passed else "fail"
            motion_note = " motion_outlier=1" if summary.motion_outlier else ""
            print(
                f"[{status}] {dcd_path} "
                f"bond={summary.max_bond_length_angstrom:.2f}A "
                f"center={summary.max_anchor_center_offset_angstrom:.2f}A "
                f"gap={summary.max_ligand_anchor_gap_angstrom:.2f}A "
                f"anchor_jump={summary.max_anchor_internal_jump_angstrom:.2f}A "
                f"lig_jump={summary.max_ligand_anchor_jump_angstrom:.2f}A"
                f"{motion_note}"
            )
            if not summary.passed:
                failures += 1
        except Exception as exc:
            print(f"[fail] {dcd_path}: {exc}")
            failures += 1

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["passed", "motion_outlier", "trajectory_dcd"], ascending=[True, False, True])
    df.to_csv(args.output_csv, index=False)
    passed = int(df["passed"].fillna(False).sum()) if "passed" in df.columns else 0
    print(
        f"Completed: passed={passed} failed={failures} total={len(dcd_paths)} "
        f"report={args.output_csv}"
    )
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
