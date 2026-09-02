#!/usr/bin/env python3
"""Summarize WT replicate stability and bound-pose QC.

Designed to run on Sherlock inside the repo root:

    module load chemistry py-openmm/8.1.1_py312
    PYTHONPATH=. python ops/slurm/cluster/check_wt_replicates.py \
      --root /scratch/users/rsatija/nnrti-mechanisms \
      --replicates 1-6
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd


def _parse_replicates(spec: Optional[str]) -> Optional[Set[int]]:
    if not spec or spec.lower() == "all":
        return None
    out: Set[int] = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(token))
    return out


def _resolve(root: Path, value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    p = Path(value)
    if p.is_absolute():
        return p if p.exists() else None
    q = root / p
    return q if q.exists() else None


def _rep_from_path(path: Path) -> Optional[int]:
    for part in path.parts:
        m = re.fullmatch(r"rep_(\d+)", part)
        if m:
            return int(m.group(1))
    m = re.search(r"rep(\d+)", path.name)
    return int(m.group(1)) if m else None


def _discover_jsons(root: Path, replicates: Optional[Set[int]]) -> List[Path]:
    candidates: Dict[str, Path] = {}
    for base in [root / "results/md_runs/wt", root / "results/md_runs/WT"]:
        if base.exists():
            for p in base.glob("rep_*/*.json"):
                rp = p.resolve()
                candidates.setdefault(str(rp).lower(), rp)
    filtered = []
    for path in candidates.values():
        if "quick" in path.name.lower():
            continue
        rep = _rep_from_path(path)
        if rep is None:
            continue
        if replicates is not None and rep not in replicates:
            continue
        filtered.append(path)
    return sorted(filtered, key=lambda p: (_rep_from_path(p) or 0, str(p)))


def _state_summary(path: Optional[Path], discard_fraction: float) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {"state_csv_exists": False}
    df = pd.read_csv(path)
    start = int(round(len(df) * discard_fraction))
    tail = df.iloc[start:].copy() if len(df) else df

    def stats(col: str, prefix: str) -> Dict[str, Any]:
        if col not in tail.columns or tail.empty:
            return {}
        s = pd.to_numeric(tail[col], errors="coerce").dropna()
        if s.empty:
            return {}
        return {
            f"{prefix}_mean": float(s.mean()),
            f"{prefix}_std": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
            f"{prefix}_min": float(s.min()),
            f"{prefix}_max": float(s.max()),
            f"{prefix}_last": float(s.iloc[-1]),
        }

    out: Dict[str, Any] = {
        "state_csv_exists": True,
        "state_rows": int(len(df)),
        "state_rows_analyzed": int(len(tail)),
    }
    out.update(stats("Temperature (K)", "temperature_k"))
    out.update(stats("Density (g/mL)", "density_g_ml"))
    out.update(stats("Potential Energy (kJ/mole)", "potential_kj_mol"))
    out.update(stats("Speed (ns/day)", "speed_ns_day"))
    return out


def _trajectory_summary(
    topo_path: Optional[Path],
    dcd_path: Optional[Path],
    ligand_resname: str,
    stride: int,
    max_frames: int,
) -> Dict[str, Any]:
    if topo_path is None or dcd_path is None or not topo_path.exists() or not dcd_path.exists():
        return {"trajectory_exists": False}

    import mdtraj as md
    from scipy.spatial import cKDTree

    traj = md.load(str(dcd_path), top=str(topo_path), stride=max(1, stride))
    if max_frames > 0 and traj.n_frames > max_frames:
        idx = np.linspace(0, traj.n_frames - 1, max_frames, dtype=int)
        traj = traj[idx]

    ca = traj.topology.select("protein and name CA")
    lig = traj.topology.select(f"resname '{ligand_resname}' and not element H")
    prot_heavy = traj.topology.select(f"protein and not resname '{ligand_resname}' and not element H")
    if len(ca) == 0:
        return {"trajectory_exists": True, "trajectory_error": "No protein CA atoms selected"}
    if len(lig) == 0:
        return {"trajectory_exists": True, "trajectory_error": f"No ligand atoms selected for {ligand_resname}"}

    ref = traj[0]
    traj.superpose(ref, atom_indices=ca, ref_atom_indices=ca)
    ca_rmsd = md.rmsd(traj, ref, atom_indices=ca) * 10.0

    lig_ref = ref.xyz[0, lig, :]
    lig_xyz = traj.xyz[:, lig, :]
    lig_rmsd = np.sqrt(np.mean((lig_xyz - lig_ref[None, :, :]) ** 2, axis=(1, 2))) * 10.0

    prot_resindex = np.array([traj.topology.atom(int(i)).residue.index for i in prot_heavy], dtype=int)
    min_dist = []
    contact_count = []
    for frame in range(traj.n_frames):
        lig_pos = traj.xyz[frame, lig, :]
        prot_pos = traj.xyz[frame, prot_heavy, :]
        tree = cKDTree(prot_pos)
        nearest, _ = tree.query(lig_pos, k=1)
        min_dist.append(float(np.min(nearest) * 10.0))
        close = tree.query_ball_point(lig_pos, r=0.4)
        residues = set()
        for hits in close:
            residues.update(int(prot_resindex[h]) for h in hits)
        contact_count.append(len(residues))

    def series_stats(values: np.ndarray, prefix: str) -> Dict[str, float]:
        values = np.asarray(values, dtype=float)
        return {
            f"{prefix}_mean": float(values.mean()),
            f"{prefix}_std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            f"{prefix}_min": float(values.min()),
            f"{prefix}_max": float(values.max()),
            f"{prefix}_last": float(values[-1]),
        }

    out: Dict[str, Any] = {
        "trajectory_exists": True,
        "trajectory_frames_sampled": int(traj.n_frames),
        "trajectory_atoms": int(traj.n_atoms),
        "protein_ca_atoms": int(len(ca)),
        "ligand_heavy_atoms": int(len(lig)),
    }
    out.update(series_stats(ca_rmsd, "ca_rmsd_a"))
    out.update(series_stats(lig_rmsd, "ligand_pose_rmsd_a"))
    out.update(series_stats(np.asarray(min_dist), "ligand_protein_min_distance_a"))
    out.update(series_stats(np.asarray(contact_count), "ligand_contact_residue_count"))
    return out


def _qc_flags(row: Dict[str, Any]) -> List[str]:
    flags = []
    if row.get("status") != "ok":
        flags.append("status_not_ok")
    completed = row.get("md_production_steps_completed")
    target = row.get("md_production_steps")
    if completed is not None and target is not None and completed < target:
        flags.append("incomplete_steps")
    if not row.get("state_csv_exists"):
        flags.append("missing_state_csv")
    if not row.get("trajectory_exists"):
        flags.append("missing_trajectory")
    if row.get("temperature_k_mean") is not None and not (295.0 <= row["temperature_k_mean"] <= 305.0):
        flags.append("temperature_mean_outside_300K_pm5")
    if row.get("density_g_ml_mean") is not None and not (0.95 <= row["density_g_ml_mean"] <= 1.08):
        flags.append("density_mean_outside_expected")
    if row.get("ligand_protein_min_distance_a_max") is not None and row["ligand_protein_min_distance_a_max"] > 6.0:
        flags.append("ligand_min_distance_large")
    if row.get("ligand_pose_rmsd_a_max") is not None and row["ligand_pose_rmsd_a_max"] > 8.0:
        flags.append("ligand_pose_rmsd_large")
    return flags


def main() -> int:
    parser = argparse.ArgumentParser(description="QC summary for WT holo MD replicates.")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--replicates", default="all", help="all, comma list, or range like 1-6")
    parser.add_argument("--ligand-resname", default="2KW")
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--max-frames", type=int, default=200)
    parser.add_argument("--discard-fraction", type=float, default=0.1)
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/wt_replicate_qc"))
    args = parser.parse_args()

    root = args.root.resolve()
    reps = _parse_replicates(args.replicates)
    jsons = _discover_jsons(root, reps)
    if not jsons:
        raise SystemExit(f"No WT replicate JSON files found under {root}")

    rows = []
    for json_path in jsons:
        data = json.loads(json_path.read_text())
        row: Dict[str, Any] = {
            "replicate": _rep_from_path(json_path),
            "json_path": str(json_path),
            "status": data.get("status"),
            "md_production_steps": data.get("md_production_steps"),
            "md_production_steps_completed": data.get("md_production_steps_completed"),
            "elapsed_seconds": data.get("elapsed_seconds"),
        }
        state_csv = _resolve(root, data.get("state_csv"))
        topo = _resolve(root, data.get("analysis_topology_pdb"))
        dcd = _resolve(root, data.get("analysis_dcd"))
        row.update(_state_summary(state_csv, args.discard_fraction))
        row.update(_trajectory_summary(topo, dcd, args.ligand_resname, args.stride, args.max_frames))
        flags = _qc_flags(row)
        row["qc_flags"] = ";".join(flags)
        row["qc_pass"] = len(flags) == 0
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("replicate")
    out_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "wt_replicate_qc_summary.csv"
    json_path = out_dir / "wt_replicate_qc_summary.json"
    df.to_csv(csv_path, index=False)
    json_path.write_text(json.dumps(rows, indent=2))

    display_cols = [
        "replicate",
        "qc_pass",
        "qc_flags",
        "md_production_steps_completed",
        "temperature_k_mean",
        "density_g_ml_mean",
        "ca_rmsd_a_mean",
        "ca_rmsd_a_max",
        "ligand_pose_rmsd_a_mean",
        "ligand_pose_rmsd_a_max",
        "ligand_protein_min_distance_a_max",
        "ligand_contact_residue_count_mean",
    ]
    existing = [c for c in display_cols if c in df.columns]
    print(df[existing].to_string(index=False))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
