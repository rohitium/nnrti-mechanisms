#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..pbc import load_mdtraj_trajectory, topology_for_analysis_dcd


@dataclass(frozen=True)
class AlignmentResult:
    input_dcd: Path
    topology_pdb: Path
    output_dcd: Path
    n_frames: int
    n_atoms: int
    n_alignment_atoms: int
    rmsd_first_frame_before_angstrom: float
    rmsd_first_frame_after_angstrom: float

    def to_dict(self) -> dict[str, object]:
        return {
            "input_dcd": str(self.input_dcd),
            "topology_pdb": str(self.topology_pdb),
            "output_dcd": str(self.output_dcd),
            "n_frames": int(self.n_frames),
            "n_atoms": int(self.n_atoms),
            "n_alignment_atoms": int(self.n_alignment_atoms),
            "rmsd_first_frame_before_angstrom": float(self.rmsd_first_frame_before_angstrom),
            "rmsd_first_frame_after_angstrom": float(self.rmsd_first_frame_after_angstrom),
        }


def _protein_ca_indices(topology) -> list[int]:
    return [int(i) for i in topology.select("protein and name CA")]


def _chain_protein_ca_signature(topology) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for chain in topology.chains:
        count = sum(1 for atom in chain.atoms if atom.residue.is_protein and atom.name == "CA")
        if count > 0:
            out.append((int(chain.index), int(count)))
    return out


def _paired_rmsd_angstrom(a_xyz, b_xyz) -> float:
    import numpy as np

    delta = a_xyz - b_xyz
    return float((np.mean(np.sum(delta * delta, axis=1)) ** 0.5) * 10.0)


def _align_one(
    dcd_path: Path,
    topo_path: Path,
    reference,
    ref_ca_idx: list[int],
    out_path: Path,
) -> AlignmentResult:
    traj = load_mdtraj_trajectory(dcd_path=dcd_path, topo_path=topo_path)
    traj_ca_idx = _protein_ca_indices(traj.topology)
    if len(traj_ca_idx) != len(ref_ca_idx):
        raise ValueError(
            f"protein CA count mismatch: traj={len(traj_ca_idx)} ref={len(ref_ca_idx)} "
            f"traj_signature={_chain_protein_ca_signature(traj.topology)} "
            f"ref_signature={_chain_protein_ca_signature(reference.topology)}"
        )

    before = _paired_rmsd_angstrom(traj.xyz[0, traj_ca_idx], reference.xyz[0, ref_ca_idx])
    traj.superpose(reference, atom_indices=traj_ca_idx, ref_atom_indices=ref_ca_idx)
    after = _paired_rmsd_angstrom(traj.xyz[0, traj_ca_idx], reference.xyz[0, ref_ca_idx])
    traj.save_dcd(str(out_path))
    return AlignmentResult(
        input_dcd=dcd_path,
        topology_pdb=topo_path,
        output_dcd=out_path,
        n_frames=int(traj.n_frames),
        n_atoms=int(traj.n_atoms),
        n_alignment_atoms=int(len(traj_ca_idx)),
        rmsd_first_frame_before_angstrom=float(before),
        rmsd_first_frame_after_angstrom=float(after),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Align corrected analysis trajectories to 4NCG using protein C-alpha atoms.",
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
        "--reference-cif",
        type=Path,
        default=Path("data/structures/4NCG.cif"),
        help="Reference structure for alignment.",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="_aligned_4ncg_ca",
        help="Suffix added before .dcd for aligned outputs.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite aligned outputs if they already exist.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("results/tables/aligned_4ncg_ca_summary.csv"),
        help="Where to write the alignment summary table.",
    )
    args = parser.parse_args()

    if not args.reference_cif.exists():
        raise FileNotFoundError(f"Missing reference CIF: {args.reference_cif}")

    import mdtraj as md

    reference = md.load(str(args.reference_cif))
    ref_ca_idx = _protein_ca_indices(reference.topology)
    if len(ref_ca_idx) < 3:
        raise ValueError(f"Reference has too few protein CA atoms: {len(ref_ca_idx)}")

    dcd_paths = sorted(args.root.rglob(args.pattern))
    if not dcd_paths:
        print(f"No trajectories found under {args.root} matching {args.pattern}")
        return 0

    rows: list[dict[str, object]] = []
    failed = 0
    for dcd_path in dcd_paths:
        topo_path = topology_for_analysis_dcd(dcd_path)
        if not topo_path.exists():
            print(f"[fail] missing topology: {topo_path}")
            failed += 1
            continue
        out_path = dcd_path.with_name(f"{dcd_path.stem}{args.suffix}.dcd")
        if out_path.exists() and not args.overwrite:
            print(f"[skip] {out_path}")
            continue
        try:
            result = _align_one(
                dcd_path=dcd_path,
                topo_path=topo_path,
                reference=reference,
                ref_ca_idx=ref_ca_idx,
                out_path=out_path,
            )
            rows.append(result.to_dict())
            print(
                f"[ok] {out_path} frames={result.n_frames} atoms={result.n_atoms} "
                f"align_atoms={result.n_alignment_atoms} "
                f"first_frame_rmsd={result.rmsd_first_frame_before_angstrom:.2f}->{result.rmsd_first_frame_after_angstrom:.2f}A"
            )
        except Exception as exc:
            print(f"[fail] {dcd_path}: {exc}")
            failed += 1

    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.summary_csv, index=False)
    print(f"Completed: ok={len(rows)} failed={failed} total={len(dcd_paths)} summary={args.summary_csv}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
