#!/usr/bin/env python3
"""Plot a 100 ns triplet story for any frame-level metric in frame_features.csv."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..pbc import load_mdtraj_trajectory
from ..susceptibility import load_dor_susceptibilities
from .compute_ligand_pocket_states import (
    _aligned_dcd_path,
    _find_reference_ligand_atom_index,
    _infer_total_ns,
    _resolve_ligand_atom_index,
)


MUTATION_COLORS = {
    "WT": "#333333",
    "Y181C": "#4c78a8",
    "Y188L": "#e45756",
    "V106I": "#4c78a8",
    "V106A": "#e45756",
    "G190A": "#4c78a8",
    "G190E": "#e45756",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a WT/comparator/DRM frame-metric story over 100 ns.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument("--reference-cif", type=Path, default=Path("data/structures/4NCG.cif"))
    parser.add_argument(
        "--contact-defs",
        type=Path,
        default=Path("results/tables/holo/dor_key_contact_definitions_4ncg.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/geometry_wt_v106a_f227l_v106i_f227c_pose_rmsd"),
    )
    parser.add_argument("--triplet", type=str, default="WT,V106A+F227L,V106I+F227C")
    parser.add_argument("--metric", type=str, default="ligand_pose_rmsd_angstrom")
    parser.add_argument("--ylabel", type=str, default="Ligand Pose RMSD (A)")
    parser.add_argument("--title", type=str, default="Ligand Pose RMSD")
    parser.add_argument("--output-prefix", type=str, default="triplet_story_100ns_WT_V106A_F227L_V106I_F227C_POSE_RMSD")
    parser.add_argument("--max-time-ns", type=float, default=100.0)
    parser.add_argument("--aligned-suffix", type=str, default="_aligned_4ncg_ca")
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--triplet-colors",
        type=str,
        default="",
        help="Optional comma-separated colors to apply to the triplet in order.",
    )
    return parser.parse_args()


def _load_fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    fold_map = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    fold_map["WT"] = 1.0
    return fold_map


def _build_full_time_axis_ns(n_frames: int, total_ns: float | None) -> np.ndarray:
    if n_frames <= 1:
        return np.zeros(max(1, n_frames), dtype=float)
    if total_ns is not None and np.isfinite(total_ns) and float(total_ns) > 0.0:
        return np.linspace(0.0, float(total_ns), int(n_frames), dtype=float)
    return np.arange(int(n_frames), dtype=float)


def _load_replicate_rows(manifest_csv: Path, needed_mutations: set[str]) -> pd.DataFrame:
    mf = pd.read_csv(manifest_csv)
    needed_cols = {"mutation", "replicate", "output_json"}
    missing = needed_cols - set(mf.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    return mf[mf["mutation"].astype(str).isin(sorted(needed_mutations))].copy()


def _compute_ligand_pose_rmsd_from_aligned(
    manifest_csv: Path,
    *,
    triplet: list[str],
    reference_cif: Path,
    contact_defs_csv: Path,
    aligned_suffix: str,
    ligand_resname: str,
    max_time_ns: float,
) -> pd.DataFrame:
    import mdtraj as md

    if not manifest_csv.exists():
        raise FileNotFoundError(manifest_csv)
    if not reference_cif.exists():
        raise FileNotFoundError(reference_cif)
    if not contact_defs_csv.exists():
        raise FileNotFoundError(contact_defs_csv)

    mf = _load_replicate_rows(manifest_csv, set(triplet))
    if mf.empty:
        raise ValueError("No requested mutations found in manifest")

    reference = md.load(str(reference_cif))
    contact_defs = pd.read_csv(contact_defs_csv)
    pose_atom_names = sorted({str(x) for x in contact_defs["ligand_atom"].dropna().astype(str)})
    ref_pose_atom_indices: list[int] = []
    ref_pose_atom_xyz: dict[str, np.ndarray] = {}
    for atom_name in pose_atom_names:
        idx = _find_reference_ligand_atom_index(reference.topology, ligand_resname, atom_name)
        if idx is None:
            raise ValueError(f"Could not map reference ligand atom {atom_name} in {reference_cif}")
        ref_pose_atom_indices.append(int(idx))
        ref_pose_atom_xyz[str(atom_name)] = reference.xyz[0, int(idx)].copy()

    repo_root = Path(__file__).resolve().parents[3]
    rows: list[pd.DataFrame] = []
    for _, row in mf.sort_values(["mutation", "replicate"]).iterrows():
        mutation = str(row["mutation"])
        replicate = int(pd.to_numeric(row["replicate"], errors="coerce"))
        out_json = Path(str(row["output_json"]))
        if not out_json.exists():
            text = str(out_json)
            marker = "nnrti-mechanisms/"
            if marker in text:
                out_json = repo_root / text.split(marker, 1)[1]
        if not out_json.exists():
            continue
        data = json.loads(out_json.read_text())
        topo_path = Path(str(data["analysis_topology_pdb"]))
        analysis_dcd = Path(str(data["analysis_dcd"]))
        if not topo_path.exists():
            topo_path = repo_root / str(topo_path)
        if not analysis_dcd.exists():
            analysis_dcd = repo_root / str(analysis_dcd)
        aligned_dcd = _aligned_dcd_path(analysis_dcd, aligned_suffix)
        if not aligned_dcd.exists():
            raise FileNotFoundError(f"Missing aligned DCD: {aligned_dcd}")
        traj = load_mdtraj_trajectory(aligned_dcd, topo_path)
        if traj.n_frames < 1:
            continue

        frame0_xyz = traj.xyz[0]
        traj_pose_atom_idx: list[int] = []
        for atom_name in pose_atom_names:
            atom_idx = _resolve_ligand_atom_index(
                traj.topology,
                xyz_frame0=frame0_xyz,
                ligand_resname=ligand_resname,
                atom_name=atom_name,
                ref_atom_xyz=ref_pose_atom_xyz[atom_name],
            )
            if atom_idx is None:
                raise ValueError(f"Could not map ligand pose atom {atom_name} in {aligned_dcd}")
            traj_pose_atom_idx.append(int(atom_idx))

        rmsd = md.rmsd(
            traj,
            reference,
            frame=0,
            atom_indices=np.asarray(traj_pose_atom_idx, dtype=int),
            ref_atom_indices=np.asarray(ref_pose_atom_indices, dtype=int),
        ) * 10.0
        total_ns = _infer_total_ns(pd.Series({"analysis_dcd": str(analysis_dcd), "safe_label": data.get("safe_label", mutation), "replicate": replicate}))
        time_ns = _build_full_time_axis_ns(int(traj.n_frames), total_ns)
        keep = time_ns <= float(max_time_ns)
        if int(np.sum(keep)) < 2:
            keep = np.ones(int(traj.n_frames), dtype=bool)
        rows.append(
            pd.DataFrame(
                {
                    "mutation": mutation,
                    "replicate": replicate,
                    "time_ns": time_ns[keep].astype(float),
                    "ligand_pose_rmsd_angstrom": rmsd[keep].astype(float),
                }
            )
        )
    if not rows:
        raise ValueError("No aligned pose-RMSD traces could be computed")
    return pd.concat(rows, ignore_index=True)


def _build_common_time_grid(trace_df: pd.DataFrame, max_time_ns: float) -> np.ndarray:
    times = pd.to_numeric(trace_df["time_ns"], errors="coerce").dropna().to_numpy(dtype=float)
    times = np.unique(np.round(times, 6))
    times = times[(times >= 0.0) & (times <= float(max_time_ns))]
    if times.size == 0:
        return np.linspace(0.0, float(max_time_ns), num=101, dtype=float)
    if times[0] > 0.0:
        times = np.concatenate([np.array([0.0], dtype=float), times])
    if times[-1] < float(max_time_ns):
        times = np.concatenate([times, np.array([float(max_time_ns)], dtype=float)])
    return times


def _interpolate_replicates_to_grid(trace_df: pd.DataFrame, metric: str, time_grid: np.ndarray) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for (mutation, replicate), rep_df in trace_df.groupby(["mutation", "replicate"], sort=True):
        rep_df = rep_df.sort_values("time_ns", kind="stable").copy()
        x = pd.to_numeric(rep_df["time_ns"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(rep_df[metric], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        x = x[ok]
        y = y[ok]
        if x.size == 0:
            continue
        x_unique, unique_idx = np.unique(x, return_index=True)
        y_unique = y[unique_idx]
        y_interp = np.interp(time_grid, x_unique, y_unique)
        rows.append(
            pd.DataFrame(
                {
                    "mutation": str(mutation),
                    "replicate": int(replicate),
                    "time_ns": time_grid.astype(float),
                    metric: y_interp.astype(float),
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["mutation", "replicate", "time_ns", metric])
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    args = _parse_args()
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    triplet = [token.strip() for token in str(args.triplet).split(",") if token.strip()]
    if len(triplet) != 3:
        raise ValueError("--triplet must contain exactly three comma-separated mutations")

    triplet_colors = [token.strip() for token in str(args.triplet_colors).split(",") if token.strip()]
    if triplet_colors and len(triplet_colors) != len(triplet):
        raise ValueError("--triplet-colors must provide exactly one color per triplet entry")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    metric = str(args.metric)
    fold_map = _load_fold_map(args.susceptibility_xlsx)
    data_source: str
    if metric == "ligand_pose_rmsd_angstrom":
        raw_trace_df = _compute_ligand_pose_rmsd_from_aligned(
            args.manifest,
            triplet=triplet,
            reference_cif=args.reference_cif,
            contact_defs_csv=args.contact_defs,
            aligned_suffix=str(args.aligned_suffix),
            ligand_resname=str(args.ligand_resname),
            max_time_ns=float(args.max_time_ns),
        )
        data_source = "full_aligned_trajectory"
    else:
        if not args.frame_feature_csv.exists():
            raise FileNotFoundError(args.frame_feature_csv)
        frame_df = pd.read_csv(args.frame_feature_csv)
        if metric not in frame_df.columns:
            raise ValueError(f"Metric {metric!r} not found in {args.frame_feature_csv}")
        raw_trace_df = frame_df[
            frame_df["mutation"].astype(str).isin(triplet)
            & (pd.to_numeric(frame_df["time_ns"], errors="coerce") <= float(args.max_time_ns))
        ][["mutation", "replicate", "time_ns", metric]].copy()
        data_source = "frame_features_csv"
    raw_trace_df["time_ns"] = pd.to_numeric(raw_trace_df["time_ns"], errors="coerce").astype(float)
    raw_trace_df = raw_trace_df.sort_values(["mutation", "replicate", "time_ns"], kind="stable").reset_index(drop=True)
    raw_trace_df.to_csv(out_tables / "trace_values_raw.csv", index=False)

    time_grid = _build_common_time_grid(raw_trace_df, max_time_ns=float(args.max_time_ns))
    trace_df = _interpolate_replicates_to_grid(raw_trace_df, metric, time_grid)
    trace_df.to_csv(out_tables / "trace_values.csv", index=False)

    mean_df = (
        trace_df.groupby(["mutation", "time_ns"], as_index=False)[metric]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "metric_mean", "std": "metric_std", "count": "n_replicates"})
    )
    mean_df["metric_sem"] = (
        mean_df["metric_std"].fillna(0.0) / mean_df["n_replicates"].clip(lower=1).pow(0.5)
    ).astype(float)
    mean_df.to_csv(out_tables / "mean_traces.csv", index=False)

    fig, ax = plt.subplots(figsize=(13.8, 5.6), constrained_layout=True)
    xmax = 0.0
    for idx, mutation in enumerate(triplet):
        color = triplet_colors[idx] if triplet_colors else MUTATION_COLORS.get(mutation, "#555555")
        mut_mean = mean_df[mean_df["mutation"].astype(str) == mutation].copy()
        x = mut_mean["time_ns"].to_numpy(dtype=float)
        y = mut_mean["metric_mean"].to_numpy(dtype=float)
        sem = mut_mean["metric_sem"].to_numpy(dtype=float)
        xmax = max(xmax, float(np.nanmax(x)) if len(x) else 0.0)
        fold = fold_map.get(mutation, float("nan"))
        label = f"{mutation} ({fold:.1f}x)" if pd.notna(fold) else str(mutation)
        ax.plot(x, y, color=color, linewidth=2.1, alpha=0.95, label=label)
        lo = y - sem
        hi = y + sem
        ok = np.isfinite(lo) & np.isfinite(hi)
        if np.any(ok):
            ax.fill_between(x[ok], lo[ok], hi[ok], color=color, alpha=0.16, linewidth=0)
    ax.set_xlim(0.0, xmax if xmax > 0 else float(args.max_time_ns))
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(str(args.ylabel))
    ax.set_title(str(args.title))
    ax.grid(alpha=0.22, linestyle=":")
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    png = out_plots / f"{str(args.output_prefix)}.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "frame_feature_csv": str(args.frame_feature_csv),
                "reference_cif": str(args.reference_cif),
                "contact_defs": str(args.contact_defs),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "triplet": triplet,
                "max_time_ns": float(args.max_time_ns),
                "metric": metric,
                "ylabel": str(args.ylabel),
                "title": str(args.title),
                "output_prefix": str(args.output_prefix),
                "triplet_colors": triplet_colors,
                "aligned_suffix": str(args.aligned_suffix),
                "ligand_resname": str(args.ligand_resname),
                "data_source": data_source,
                "aggregation": "replicates interpolated onto common time grid before mean/SEM",
                "common_time_grid_size": int(len(time_grid)),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
